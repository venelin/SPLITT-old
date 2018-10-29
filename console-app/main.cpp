#include "SPLITT.h"
#include <iostream>
#include <fstream>


using namespace SPLITT;

template<class NameType>
struct NumericTraitData {
  // use const references to avoid copying of long vectors
  std::vector<NameType> const& names_;
  vec const& x_;
  NumericTraitData(
    std::vector<NameType> const& names,
    vec const& x): names_(names), x_(x) {}
};


template<class Tree>
class AbcPMM: public TraversalSpecification<Tree> {
  
public:
  typedef AbcPMM<Tree> MyType;
  typedef TraversalSpecification<Tree> BaseType;
  typedef Tree TreeType;
  typedef PostOrderTraversal<MyType> AlgorithmType;
  typedef vec ParameterType;
  typedef NumericTraitData<typename TreeType::NodeType> DataType;
  typedef vec StateType;
  
  double sigmae2, sigma2;
  vec x;
  vec a, b, c;
  
  AbcPMM(TreeType const& tree, DataType const& input_data):
    BaseType(tree) {
    
    if(input_data.x_.size() != this->ref_tree_.num_tips()) {
      std::ostringstream oss;
      oss<<"The vector x must be the same length as the number of tips ("<<
        this->ref_tree_.num_tips()<<"), but were"<<input_data.x_.size()<<".";
      throw std::invalid_argument(oss.str());
    } else {
      
      uvec ordNodes = this->ref_tree_.OrderNodes(input_data.names_);
      this->x = At(input_data.x_, ordNodes);
      this->a = vec(this->ref_tree_.num_nodes());
      this->b = vec(this->ref_tree_.num_nodes());
      this->c = vec(this->ref_tree_.num_nodes());
    }
  };
  
  StateType StateAtRoot() const {
    vec res(3);
    res[0] = a[this->ref_tree_.num_nodes() - 1];
    res[1] = b[this->ref_tree_.num_nodes() - 1];
    res[2] = c[this->ref_tree_.num_nodes() - 1];
    return res;
  };
  
  void SetParameter(ParameterType const& par) {
    if(par.size() != 2) {
      throw std::invalid_argument(
          "The par vector should be of length 2 with \
      elements corresponding to sigma2 and sigmae2.");
    }
    if(par[0] <= 0 || par[1] <= 0) {
      throw std::logic_error("The parameters sigma2 and sigmae2 should be positive.");
    }
    this->sigma2 = par[0];
    this->sigmae2 = par[1];
  }
  
  inline void InitNode(uint i) {
    
    if(i < this->ref_tree_.num_tips()) {
      a[i] = -0.5 / sigmae2;  
      b[i] = x[i] / sigmae2;
      c[i] = -0.5 * (x[i]*b[i] + log(2*G_PI*sigmae2));
    } else {
      a[i] = b[i] = c[i] = 0;
    }
  }
  
  inline void VisitNode(uint i) {
    double t = this->ref_tree_.LengthOfBranch(i);
    
    double d = 1 - 2*a[i]*sigma2*t;
    // the order is important here because for c[i] we use the previous values 
    // of a[i] and b[i].
    c[i] = c[i] - 0.5*log(d) + 0.5*b[i]*b[i]*sigma2*t/d;
    a[i] /= d;
    b[i] /= d;
  }
  
  inline void PruneNode(uint i, uint j) {
    a[j] = a[j] + a[i];
    b[j] = b[j] + b[i];
    c[j] = c[j] + c[i];
  }
  
};

typedef TraversalTask<AbcPMM<OrderedTree<uint, double>> > ParallelPruningAbcPMM;


int main() {
  
  std::fstream fin;
  
  std::cout<<"Hello from the SPLITT PMM example!"<<std::endl;
  std::cout<<"Reading the input tree, data and PMM parameters from the file test-input.txt..."<<std::endl;
  fin.open("test-input.txt", std::fstream::in);
  
  uint M, N;
  fin>>M>>N;
  uvec daughters(M-1);
  uvec parents(M-1);
  vec t(M-1);
  vec x(N);
  
  for(int i=0; i < M-1; i++) {
    fin>>daughters[i]>>parents[i]>>t[i];
    if(daughters[i] <= N) {
      fin>>x[i];
    }
  }
  
  double sigma, sigmae;
  
  fin>>sigma>>sigmae;
  
  vec param(2); 
  
  param[0] = sigma*sigma;
  param[1] = sigmae*sigmae;
  
  uvec node_names = Seq(uint(1), N);
  typename ParallelPruningAbcPMM::DataType data(node_names, x);
  
  ParallelPruningAbcPMM pruningTask(parents, daughters, t, data);
   
   
  vec abc = pruningTask.TraverseTree(param, 0);
  
  std::cout<<"Calculating the a, b, c at the root for sigma="<<sigma<<" and sigmae="<<sigmae<<std::endl;
  std::cout<<"a="<<abc[0]<<", b="<<abc[1]<<", c="<<abc[2]<<std::endl;
  std::cout<<"OpenMP-version: "<<pruningTask.algorithm().VersionOPENMP()<<std::endl;
  std::cout<<"Good bye!"<<std::endl;
  return 0;
}

