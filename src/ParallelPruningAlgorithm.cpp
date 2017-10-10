#include <RcppArmadillo.h>
#include <algorithm>
#include <vector>
#include <math.h>
#include <R_ext/Rdynload.h>


//[[Rcpp::plugins("cpp11")]]
#ifdef _OPENMP

// Need to decide wheter to use '#pragma omp for' or '#pragma omp for simd'
#if _OPENMP >= 201307  // OMP 4.0 or higher

#define _PRAGMA_OMP_FOR_SIMD _Pragma("omp for simd")
#define _PRAGMA_OMP_FOR _Pragma("omp for")
#define _PRAGMA_OMP_SIMD _Pragma("omp simd")

#else // #if _OPENMP >= 201307

#define _PRAGMA_OMP_FOR_SIMD _Pragma("omp for")
#define _PRAGMA_OMP_FOR _Pragma("omp for")
#define _PRAGMA_OMP_SIMD _Pragma("omp simd")

#endif // _OPENMP >= 201307

#include <omp.h>

#else // #ifdef _OPENMP

// the preprocessor directives should simply be ignored at compile-time
#define _PRAGMA_OMP_FOR_SIMD _Pragma("omp for simd")
#define _PRAGMA_OMP_FOR _Pragma("omp for")
#define _PRAGMA_OMP_SIMD _Pragma("omp simd")

#endif // #ifdef _OPENMP

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;
using namespace std;

// For each i a[i] with b[i] in x
// [[Rcpp::export]]
arma::uvec multiReplaceC(arma::uvec const& x, arma::uvec const& a, arma::uvec const& b) {
  arma::uvec x2 = x;
  arma::uvec ind = sort_index(x2);
  arma::uvec xInd = x2(ind);
  std::pair<arma::uvec::iterator, arma::uvec::iterator> bounds;
  for(int i = 0; i < a.n_elem; ++i) {
    bounds = std::equal_range(xInd.begin(), xInd.end(), a.at(i));
    if(bounds.first != bounds.second) {
      uint first = bounds.first - xInd.begin();
      uint last = bounds.second - xInd.begin() - 1;
      x2.elem(ind.subvec(first, last)).fill(b.at(i));
    }
  }
  return x2;
}

class ParallelPruningAlgorithm {
protected:
  uint nLevels;
  // number of all nodes
  uint M;
  // number of tips
  uint N;

  // edge matrix
  arma::umat edge;

  arma::vec tReord;
  arma::uvec eReord;
  arma::uvec reord;

  arma::uvec tipsVector;
  arma::uvec tipsVectorIndex;
  arma::uvec edgeVector;
  arma::uvec edgeVectorIndex;

public:
  // Default constructor;
  ParallelPruningAlgorithm() {
  };

  uint get_nLevels() const {
    return nLevels;
  }

  uint get_M() const {
    return M;
  }

  uint get_N() const {
    return N;
  }

  arma::vec get_tReord() const {
    return (this->tReord);
  }

  arma::uvec get_eReord() const {
    return (this->eReord);
  }

  arma::uvec get_tipsVector() const {
    return tipsVector;
  }

  arma::uvec get_tipsVectorIndex() const {
    return tipsVectorIndex;
  }

  arma::umat get_edge() const {
    return edge;
  }

  arma::uvec get_edgeVector() const {
    return edgeVector;
  }

  arma::uvec get_edgeVectorIndex() const {
    return edgeVectorIndex;
  }

  void createPruningOrder(Rcpp::List const& tree) {
    if(!tree.inherits("phylo")) {
      stop("Input must be a phylo object.");
    }

    // 0-based indices
    IntegerMatrix edge = as<IntegerMatrix>(tree["edge"]) - 1;

    // number of tips
    this->N = as<CharacterVector>(tree["tip.label"]).length();
    // number of all nodes
    this->M = unique(as<IntegerVector>(edge)).length();

    IntegerVector edgeEnds = edge(_, 1);
    edgeEnds.push_back(N);
    IntegerVector endingAt = match(seq(0, M - 1), edgeEnds) - 1;

    IntegerVector nonPrunedChildren(M);
    IntegerVector ee1 = edge(_, 0);

    while(ee1.length() > 0) {
      IntegerVector matchp = match(seq(N, M - 1), ee1) - 1;
      matchp = matchp[!Rcpp::is_na(matchp)];

      for(int m : matchp) {
        nonPrunedChildren[ee1[m]]++;
      }

      ee1[matchp] = -1;
      ee1 = ee1[ee1 != (-1)];
    }

    IntegerVector tipsVector;
    IntegerVector tipsVectorIndex = IntegerVector(1);

    IntegerVector edgeVector;
    IntegerVector edgeVectorIndex = IntegerVector(1);

    // start by pruning the tips
    IntegerVector tips = seq(0, N - 1);

    while(tips[0] != N) { // while the root has not become a tip itself

      tipsVectorIndex.push_back(
        tipsVectorIndex[tipsVectorIndex.length() - 1] + tips.length());

      // add the tips to be pruned to the tipsVector
      for(int n : tips) tipsVector.push_back(n);

      // indices in the edge-matrix of the to be pruned edges (pointing to tips)
      IntegerVector edgesToTips = endingAt[tips];

      // empty the tips vector so it can be filled in with new tips to be pruned
      tips = IntegerVector();

      int nEdgesDone = 0;
      IntegerVector remainingParents;

      while(nEdgesDone != edgesToTips.length() ) {
        // vectorized update of the parent nodes.
        // resolving order of sibling edges being pruned at the same time
        if(nEdgesDone == 0) {
          remainingParents = edge(_, 0);
          remainingParents = remainingParents[edgesToTips];
        } else {
          remainingParents = edge(_, 0);
          remainingParents = remainingParents[
          as<IntegerVector>(edgesToTips[! is_na(edgesToTips)])];
        }

        remainingParents = unique(remainingParents);

        IntegerVector edgeStarts;
        for(int ett : edgesToTips) {
          if(ett != NA_INTEGER) {
            edgeStarts.push_back(edge(ett, 0));
          } else {
            edgeStarts.push_back(NA_INTEGER);
          }
        }

        // sib- edges that are first in edge[edgesToTips, 0] get served
        // first
        IntegerVector edgesNext = match(remainingParents, edgeStarts);
        edgesNext = edgesNext - 1;
        edgesNext = edgesNext.sort();

        for(int u : edgesNext) edgeVector.push_back(u);

        // attach the index of the current last element of edgeVector
        // to edgeVectorIndex
        edgeVectorIndex.push_back(
          edgeVectorIndex[edgeVectorIndex.length() - 1] + edgesNext.length());

        // For the parent nodes, decrement the amount of non-visited
        // children
        for(int u : edgesNext) {
          nonPrunedChildren[edge(edgesToTips[u], 0)]--;
          if(nonPrunedChildren[edge(edgesToTips[u], 0)] == 0) {
            tips.push_back(edge(edgesToTips[u], 0));
          }
        }

        for(int u : edgesNext) edgesToTips[u] = NA_INTEGER;

        nEdgesDone += edgesNext.length();
      }
    }

    this->nLevels = tipsVectorIndex.length() - 1;
    this->tipsVector = as<arma::uvec>(tipsVector);
    this->tipsVectorIndex = as<arma::uvec>(tipsVectorIndex);
    this->edgeVector = as<arma::uvec>(edgeVector);
    this->edgeVectorIndex = as<arma::uvec>(edgeVectorIndex);
    this->edge = as<arma::umat>(edge);

    NumericVector edgeLength = tree["edge.length"];
    this->tReord = as<arma::vec>(edgeLength);

    IntegerVector ord = seq(0, M - 1);
    this->reord = as<arma::uvec>(ord);

    reorderEdges();
  };

  void reorderEdges() {
    IntegerVector edgeEnds;
    for(int i = 0; i<edge.n_rows; i++) edgeEnds.push_back(edge(i, 1));
    edgeEnds.push_back(N);
    IntegerVector endingAtRcpp = match(seq(0, M - 1), edgeEnds) - 1;
    arma::uvec endingAt = as<arma::uvec>(endingAtRcpp);

    arma::uvec ONE(1);
    ONE.fill(1);
    arma::uvec parents = edge.col(0);
    parents.replace(N, 2*M - 1);
    (this->eReord) = parents;
    // duplicate tReord and reord since it will be shuffled during the reordering
    arma::vec t = (this->tReord);
    arma::uvec ord = (this->reord);


    // edges pointing to tips
    uvec edgesToTips = endingAt.elem(tipsVector(span(tipsVectorIndex(0),
                                                     tipsVectorIndex(1) - 1)));

    uint jEVI = 0;

    uint nEdgesDone = 0;
    while(nEdgesDone != edgesToTips.n_elem) {
      uvec edgesNext = edgeVector(span(edgeVectorIndex(jEVI),
                                       edgeVectorIndex(jEVI + 1) - 1));

      uvec edgeEnds = edge(edgesToTips(edgesNext), ONE);

      uvec edgeEndsNew = arma::regspace<uvec>(edgeVectorIndex(jEVI),
                                              edgeVectorIndex(jEVI + 1) - 1);

      parents = multiReplaceC(parents, edgeEnds, M+edgeEndsNew);

      (this->eReord).subvec(edgeVectorIndex(jEVI), edgeVectorIndex(jEVI + 1) - 1) =
        parents(edgesToTips(edgesNext));
      (this->eReord) = multiReplaceC((this->eReord), edgeEnds, M+edgeEndsNew);

      (this->tReord).subvec(edgeVectorIndex(jEVI), edgeVectorIndex(jEVI + 1) - 1) =
        t.elem(edgesToTips(edgesNext));

      (this->reord).subvec(edgeVectorIndex(jEVI), edgeVectorIndex(jEVI + 1) - 1) =
        ord(edgeEnds);

      ++jEVI;
      nEdgesDone += edgesNext.n_elem;
    }

    // edges pointing to internal nodes that have become tips
    for(int i = 1; i < nLevels; ++i) {
      edgesToTips = endingAt.elem(tipsVector(span(tipsVectorIndex(i),
                                                  tipsVectorIndex(i + 1) - 1)));

      uint nEdgesDone = 0;
      while(nEdgesDone != edgesToTips.n_elem) {
        uvec edgesNext = edgeVector(span(edgeVectorIndex(jEVI),
                                         edgeVectorIndex(jEVI + 1) - 1));

        //cout<<"edgesNext"<<edgesNext<<endl;

        uvec edgeEnds = edge(edgesToTips(edgesNext), ONE);

        //cout<<"edgeEnds:"<<edgeEnds<<endl;


        uvec edgeEndsNew = arma::regspace<uvec>(edgeVectorIndex(jEVI),
                                                edgeVectorIndex(jEVI + 1) - 1);

        //cout<<"edgeEndsNew:"<<edgeEndsNew<<endl;

        parents = multiReplaceC(parents, edgeEnds, M+edgeEndsNew);

        //cout<<"parents:"<<parents<<endl;

        (this->eReord).subvec(edgeVectorIndex(jEVI), edgeVectorIndex(jEVI + 1) - 1) =
          parents(edgesToTips(edgesNext));
        (this->eReord) = multiReplaceC((this->eReord), edgeEnds, M+edgeEndsNew);

        //cout<<"eReord:"<<(this->eReord)<<endl;

        (this->tReord).subvec(edgeVectorIndex(jEVI), edgeVectorIndex(jEVI + 1) - 1) =
          t.elem(edgesToTips(edgesNext));

        //cout<<"tReord:"<<(this->tReord)<<endl;

        (this->reord).subvec(edgeVectorIndex(jEVI), edgeVectorIndex(jEVI + 1) - 1) =
          ord(edgeEnds);

        //cout<<"reord:"<<(this->reord)<<endl;

        ++jEVI;
        nEdgesDone += edgesNext.n_elem;
      }
    }

    (this->eReord) = (this->eReord) - M;
  };

  virtual void set_parameters(Rcpp::List const& par) {
    //cout<<"setParameters()"<<endl;
  };

  virtual void prepareBranch(uint i) {
    //cout<<"preapareCachForBranch("<<i<<")"<<endl;
  };

  virtual void pruneBranch(uint i) {
    //cout<<"pruneBranch("<<i<<")"<<endl;
  };
  virtual void addToParent(uint i) {
    //cout<<"addToParent("<<i<<")"<<endl;
  };

  virtual ~ParallelPruningAlgorithm() {};

  void do_pruning(Rcpp::List const& par) {
    set_parameters(par);

#pragma omp parallel
{
  _PRAGMA_OMP_FOR_SIMD
  for(uint i = 0; i < M - 1; i++) {
    // calculated and store cache information for branch i
    prepareBranch(i);
  }
  uint jEVI = 0;

  for(int j = 0; j < nLevels; j++) {
    const uint eFirst = tipsVectorIndex(j);
    const uint eLast = tipsVectorIndex(j + 1) - 1;

    _PRAGMA_OMP_FOR_SIMD
      for(uint i = eFirst; i < eLast + 1; i++) {
        // perform the main calculation for branch i based on the pruning
        // result from its daughter branches.
        pruneBranch(i);
      }

      uint nEdgesDone = 0;
    while(nEdgesDone != eLast - eFirst + 1) {
      const uint unFirst = edgeVectorIndex(jEVI);
      const uint unLast = edgeVectorIndex(jEVI + 1) - 1;
      _PRAGMA_OMP_FOR_SIMD
      for(uint i = unFirst; i < unLast + 1; i++) {
        // store or add up the result from branch i to the results from its
        // sibling branches, so that these results can be used for the
        // pruning of the parent branch.
        addToParent(i);
      }
      nEdgesDone +=  edgeVectorIndex(jEVI + 1) - edgeVectorIndex(jEVI);
      ++jEVI;
    }
  }
}
  };
};

RCPP_MODULE(ParallelPruningAlgorithm) {
  class_<ParallelPruningAlgorithm>( "ParallelPruningAlgorithm" )
  .constructor()
  .method( "createPruningOrder", &ParallelPruningAlgorithm::createPruningOrder )
  .method( "reorderEdges", &ParallelPruningAlgorithm::reorderEdges )
  .method( "do_pruning", &ParallelPruningAlgorithm::do_pruning )
  .property("nLevels", &ParallelPruningAlgorithm::get_nLevels )
  .property("M", &ParallelPruningAlgorithm::get_M )
  .property("N", &ParallelPruningAlgorithm::get_N )
  .property("tReord", &ParallelPruningAlgorithm::get_tReord )
  .property("eReord", &ParallelPruningAlgorithm::get_eReord )
  .property("edge", &ParallelPruningAlgorithm::get_edge )
  .property("tipsVector", &ParallelPruningAlgorithm::get_tipsVector )
  .property("tipsVectorIndex", &ParallelPruningAlgorithm::get_tipsVectorIndex )
  .property("edgeVector", &ParallelPruningAlgorithm::get_edgeVector )
  .property("edgeVectorIndex", &ParallelPruningAlgorithm::get_edgeVectorIndex )
  ;
}

class POUMM_Likelihood: public ParallelPruningAlgorithm {
protected:
  arma::vec zReord;
  arma::vec seReord;

  bool has_se;

  arma::vec sum_se2_sigmae2;

  arma::mat abcMat;
  arma::vec a, b, c;

  uint count_abc_calls;


  double g0, alpha, theta, sigma, sigmae;
  double sigmae2, sigma2, logsigma;

  arma::vec log_se_total;
  arma::vec talpha;
  arma::vec etalpha;
  arma::vec e2talpha;
  arma::vec fe2talpha;
  arma::vec gutalphasigma2;

  arma::vec z1, z1z1;

public:
  void do_pruning(Rcpp::List const& par) {
    ParallelPruningAlgorithm::do_pruning(par);
  };

  arma::vec get_abc() const {
    arma::vec res(3);
    res(0) = a(M - 1);
    res(1) = b(M - 1);
    res(2) = c(M - 1);
    return res;
  };

  arma::mat get_abcMat() const {
    arma::mat res(M, 3);
    res.col(0) = this->a;
    res.col(1) = this->b;
    res.col(2) = this->c;
    return res;
  };

  void set_treeAndData(Rcpp::List const& tree,
                       Rcpp::NumericVector const& z,
                       Rcpp::NumericVector const& se) {
    createPruningOrder(tree);
    if(z.length() != N) {
      stop("The trait vector z must be the same length as the number of tips.");
    } else {
      arma::vec zArma = z;
      this -> zReord = zArma((this->reord)(span(0, N - 1)));
      arma::vec seArma = se;
      this -> seReord = seArma((this->reord)(span(0, N - 1)));

      this -> a = arma::vec(M);
      this -> b = arma::vec(M);
      this -> c = arma::vec(M);
      this -> talpha = arma::vec(M - 1);
      this -> etalpha = arma::vec(M - 1);
      this -> e2talpha = arma::vec(M - 1);
      this -> fe2talpha = arma::vec(M - 1);
      this -> gutalphasigma2 = arma::vec(M - 1);
      this -> z1 = arma::vec(N);
      this -> z1z1 = arma::vec(N);

    }
  };

  void set_parameters(Rcpp::List const& par) {
    this -> g0 = par["g0"];
    this -> alpha = par["alpha"];
    this -> theta = par["theta"];
    this -> sigma = par["sigma"];
    this -> sigmae = par["sigmae"];

    this -> sigmae2 = sigmae*sigmae;
    this -> sum_se2_sigmae2 = sigmae2 + seReord % seReord;

    this -> log_se_total = log(sqrt(sum_se2_sigmae2));

    this -> sigma2 = sigma*sigma;
    this -> logsigma = log(sigma);

    this -> a.fill(0);
    this -> b.fill(0);
    this -> c.fill(0);
  };

  inline void prepareBranch(uint i) {
    ParallelPruningAlgorithm::prepareBranch(i);
    if(alpha != 0) {
      talpha[i] = (this->tReord)[i] * alpha;
      etalpha[i] = exp(talpha[i]);
      e2talpha[i] = etalpha[i] * etalpha[i];
      fe2talpha[i] = alpha / (1 - e2talpha[i]);
    } else {
      talpha[i] = (this->tReord)[i] * alpha;
      etalpha[i] = exp(talpha[i]);
      e2talpha[i] = etalpha[i] * etalpha[i];
      fe2talpha[i] = -0.5 / (this->tReord)[i];
    }
  };

  inline void pruneBranch(uint i) {
    ParallelPruningAlgorithm::pruneBranch(i);
    if(i < N) {
      // branch leading to a tip
      double z1, z1z1;
      gutalphasigma2[i] = e2talpha[i] + ((-0.5 / sum_se2_sigmae2[i]) * sigma2) / fe2talpha[i];
      z1 = zReord[i] - theta;
      z1z1 = z1 * z1;

      // integration over g1 including e1 = z1 - g1
      c[i] = -0.5 * log(gutalphasigma2[i]) -
        0.25 * sigma2 * z1z1 / (sum_se2_sigmae2[i]*sum_se2_sigmae2[i]) /
          (fe2talpha[i] - alpha + (-0.5 / sum_se2_sigmae2[i]) * sigma2) +
            talpha[i] + (-0.5 * (M_LN_2PI  + z1z1 / sum_se2_sigmae2[i]) - log_se_total[i]);
      b[i] = (etalpha[i] * (z1 / sum_se2_sigmae2[i])) / gutalphasigma2[i];
      a[i] = (-0.5 / sum_se2_sigmae2[i]) / gutalphasigma2[i];
    } else {
      gutalphasigma2[i] = e2talpha[i] + (a[i] * sigma2) / fe2talpha[i];
      c[i] = -0.5 * log(gutalphasigma2[i]) - 0.25 * sigma2 * b[i] * b[i] /
        (fe2talpha[i] - alpha + a[i] * sigma2) + talpha[i] + c[i];
      b[i] = (etalpha[i] * b[i]) / gutalphasigma2[i];
      a[i] /= gutalphasigma2[i];
    }
  };
  inline void addToParent(uint i) {
    ParallelPruningAlgorithm::addToParent(i);
    a[(this->eReord)[i]] += a[i];
    b[(this->eReord)[i]] += b[i];
    c[(this->eReord)[i]] += c[i];
  }
};

RCPP_MODULE(POUMM_Likelihood) {
  class_<POUMM_Likelihood>( "POUMM_Likelihood" )
  .constructor()
  .method( "set_treeAndData", &POUMM_Likelihood::set_treeAndData )
  .method( "abcMat", &POUMM_Likelihood::get_abcMat )
  .method( "abc", &POUMM_Likelihood::get_abc )
  .method( "do_pruning", &POUMM_Likelihood::do_pruning )
  .method( "set_parameters", &POUMM_Likelihood::set_parameters )
  ;
  class_<POUMM_Likelihood::ParallelPruningAlgorithm>( "POUMM_Likelihood::ParallelPruningAlgorithm" )
  //.constructor()
  .method( "createPruningOrder", &POUMM_Likelihood::ParallelPruningAlgorithm::createPruningOrder )
  .method( "reorderEdges", &POUMM_Likelihood::ParallelPruningAlgorithm::reorderEdges )
  .method( "do_pruning", &POUMM_Likelihood::ParallelPruningAlgorithm::do_pruning )
  .property("nLevels", &POUMM_Likelihood::ParallelPruningAlgorithm::get_nLevels )
  .property("M", &POUMM_Likelihood::ParallelPruningAlgorithm::get_M )
  .property("N", &POUMM_Likelihood::ParallelPruningAlgorithm::get_N )
  .property("tReord", &POUMM_Likelihood::ParallelPruningAlgorithm::get_tReord )
  .property("eReord", &POUMM_Likelihood::ParallelPruningAlgorithm::get_eReord )
  .property("edge", &POUMM_Likelihood::ParallelPruningAlgorithm::get_edge )
  .property("tipsVector", &POUMM_Likelihood::ParallelPruningAlgorithm::get_tipsVector )
  .property("tipsVectorIndex", &POUMM_Likelihood::ParallelPruningAlgorithm::get_tipsVectorIndex )
  .property("edgeVector", &POUMM_Likelihood::ParallelPruningAlgorithm::get_edgeVector )
  .property("edgeVectorIndex", &POUMM_Likelihood::ParallelPruningAlgorithm::get_edgeVectorIndex )
  ;
}
