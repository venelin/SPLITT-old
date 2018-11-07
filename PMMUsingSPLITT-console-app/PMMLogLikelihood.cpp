#include "SPLITT.h"
#include <iostream>


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
  typedef vec ParameterType;
  typedef NumericTraitData<typename TreeType::NodeType> DataType;
  typedef vec StateType;
  typedef PostOrderTraversal<MyType> AlgorithmType;
  
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
  
  StateType StateAtRoot() const {
    vec res(3);
    res[0] = a[this->ref_tree_.num_nodes() - 1];
    res[1] = b[this->ref_tree_.num_nodes() - 1];
    res[2] = c[this->ref_tree_.num_nodes() - 1];
    return res;
  };
  
};

typedef TraversalTask< AbcPMM<OrderedTree<std::string, double>> > ParallelPruningAbcPMM;

int main(int argc, char* argv[]) {
  
  // Will be false if the program is run without any arguments. Pass any argument
  // after the program name and verbose will be set to true.
  bool verbose = (argc > 1); 
  
  // will be using std::string, std::vector, std::cin and std::cout quite a lot.
  using namespace std;
  
  cout<<"Hello from the SPLITT PMM example!"<<endl;
 
  cout<<"Reading the input tree..."<<endl;
  
  // Read the number of nodes and tips
  uint M, N;
  cin>>M>>N;
  
  // read the tree
  vector<string> daughters(M-1);
  vector<string> parents(M-1);
  vec t(M-1);
  
  
  for(int i=0; i < M-1; i++) {
    cin>>daughters[i]>>parents[i]>>t[i];
    if(verbose) {
      cout<<daughters[i]<<" "<<parents[i]<<" "<<t[i]<<endl;
    }
  }
  
  cout<<"Reading the trait data..."<<endl;
  // read the trait data
  vector<string> tip_names(N);
  vec x(N);
  
  
  for(int i = 0; i < N; i++) {
    cin>>tip_names[i]>>x[i];
    if(verbose) {
      cout<<tip_names[i]<<" "<<x[i]<<endl;
    }
  }
  
  // create the data-object
  typename ParallelPruningAbcPMM::DataType data(tip_names, x);
  
  // Create the TraversalTask object: this will create the OrderedTree and the 
  // AbcPMM objects in the memory.
  ParallelPruningAbcPMM pruningTask(parents, daughters, t, data);
  
  // The tree objects reorders the nodes, so that the tips are indexed from
  // 0 to N-1, the internal nodes are from N to M-2, and the root is M-1.
  if(verbose) {
    cout<<"Node-name order in the OrderedTree:"<<endl;
    for(int i = 0; i < pruningTask.tree().num_nodes(); i++) {
      cout<<i<<". "<<pruningTask.tree().FindNodeWithId(i)<<endl;
    } 
  }
  
  // Check if OPENMP is enabled:
  std::cout<<"OpenMP-version: "<<pruningTask.algorithm().VersionOPENMP()<<std::endl;
  
  // model parameters
  double gM, sigma, sigmae;
  vec param(2); 
  
  
  cout<<"Main loop..."<<endl;
  
  // do this loop as long as parameters can be read from the standard input
  while( cin>>gM>>sigma>>sigmae ) {
    param[0] = sigma*sigma;
    param[1] = sigmae*sigmae;
  
    cout<<"  Calculating a, b, c at the root for sigma="<<sigma<<" and sigmae="<<sigmae<<std::endl;
    vec abc = pruningTask.TraverseTree(param, 0);
    cout<<"  a="<<abc[0]<<", b="<<abc[1]<<", c="<<abc[2]<<std::endl;
    
    cout<<"  LL(gM="<<gM<<", sigma="<<sigma<<", sigmae="<<sigmae<<"): "<<abc[0]*gM*gM + abc[1]*gM + abc[2]<<endl;
  } 
  
  cout<<"Main loop done."<<endl;
  // Exit politely
  std::cout<<"Good bye!"<<std::endl;
  return 0;
}

