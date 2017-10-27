#ifndef POUMM_lnDetV_Q_1D_H_
#define POUMM_lnDetV_Q_1D_H_

#include "ThreePointV_lnDetV_Q_1d.h"

namespace ppa {
template<class Node>
class POUMM_lnDetV_Q_1d: public ThreePointV_lnDetV_Q_1d<Node> {
  ParallelPruningAlgorithm<ParallelPruningTree<Node, double>, POUMM_lnDetV_Q_1d> ppalgorithm;
public:
  // univariate trait vector
  ppa::vec z;
  ppa::vec h, u;
  double g0, alpha, theta, sigma, sigmae, e2alphaT;
  double T; // tree height


  POUMM_lnDetV_Q_1d(
    std::vector<Node> const& brStarts, std::vector<Node> const& brEnds,
    vec const& t,
    std::vector<Node> const& keys, vec const& z):
    ThreePointV_lnDetV_Q_1d<Node>(brStarts, brEnds, t),
    ppalgorithm(this->pptree, *this), z(z) {

    if(z.size() != this->pptree.num_tips() ) {
      throw std::invalid_argument("The trait vector must have N elements (N is the number of tips).");
    } else {
      uvec ordNodes = this->pptree.order_nodes(keys);
      this->z = At(z, ordNodes);

      // vec z2 = z;
      // for(int i = 0; i < this->pptree.N; ++i) {
      //   this->z[i] = z2[this->pptree.orderNodes[i]];
      // }
      vec X(this->pptree.num_tips());
      this->set_X_and_Y(X, X);
    }

    this->h = this->pptree.get_nodeHeights();

    this->T = *std::max_element(h.begin(), h.begin()+this->pptree.num_tips());
    this->u = ppa::vec(this->pptree.num_tips());
    for(int i = 0; i < this->pptree.num_tips(); i++) u[i] = T - h[i];
  }

  void set_parameters(double g0, double alpha, double theta, double sigma, double sigmae) {
    this->g0 = g0;
    this->alpha = alpha;
    this->theta = theta;
    this->sigma = sigma;
    this->sigmae = sigmae;
    this->e2alphaT = exp(-2*alpha*T);
  }

  void prepareBranch(uint i) {
    uint iParent = this->pptree.ParentPruneIndex(i);
    // tTransf[i] =
    //   sigma*sigma/(2*alpha) * ( (1-exp(-2*alpha*h[i]))*exp(-2*alpha*(T-h[i])) -
    //     (1-exp(-2*alpha*h[iParent]))*exp(-2*alpha*(T-h[iParent])) );

    double ealphahi = exp(alpha*h[i]);

    this->tTransf[i] = sigma*sigma/(2*alpha) *
      (e2alphaT*(ealphahi*ealphahi - exp(2*alpha*h[iParent])));

    if(i < this->pptree.num_tips()) {
      //double mu = exp(-alpha*h[i])*g0 + (1-exp(-alpha*h[i]))*theta;
      double mu = g0 / ealphahi + (1 - 1/ealphahi)*theta;

      double ealphaui = exp(alpha*u[i]);

      this->X[i] = this->Y[i] = (z[i] - mu)/ealphaui;

      this->tTransf[i] += sigmae*sigmae / (ealphaui*ealphaui);
    }
  }

  double get_lnDetD() const {
    double sum_u = 0;
    for(auto uu : u) sum_u += uu;
    return alpha*sum_u;
  }

  void do_pruning(int mode) {
    ppalgorithm.do_pruning(mode);
  }
};
};

#endif // POUMM_lnDetV_Q_1D_H_
