#ifndef ThreePointPOUMM_H_
#define ThreePointPOUMM_H_

#include "ThreePointUnivariate.h"
#include "NumericTraitData.h"

namespace ppa {

template<class Tree>
class ThreePointPOUMM: public ThreePointUnivariate<Tree> {

public:
  typedef Tree TreeType;
  typedef ThreePointUnivariate<TreeType> BaseType;
  typedef vec NodeStateType;
  typedef vec ParameterType;
  typedef NumericTraitData<typename TreeType::NodeType> InputDataType;

  // univariate trait vector
  ppa::vec z;
  ppa::vec h, u;
  double g0, alpha, theta, sigma, sigmae, e2alphaT, sum_u;
  double T; // tree height


  ThreePointPOUMM(
    TreeType const& tree, InputDataType const& input_data):
    BaseType(tree) {

    if(input_data.z_.size() != this->ref_tree_.num_tips() ||
       input_data.se_.size() != this->ref_tree_.num_tips()) {
      throw std::invalid_argument("The vectors z and se must be the same length as the number of tips.");
    } else {

      uvec ordNodes = this->ref_tree_.OrderNodes(input_data.names_);
      this->z = At(input_data.z_, ordNodes);
      vec X(this->ref_tree_.num_tips());
      this->set_X_and_Y(X, X);

      this->h = this->ref_tree_.CalculateHeights(0);

      this->T = *std::max_element(h.begin(), h.begin()+this->ref_tree_.num_tips());
      this->u = ppa::vec(this->ref_tree_.num_tips());
      for(int i = 0; i < this->ref_tree_.num_tips(); i++) u[i] = T - h[i];
      sum_u = 0;
      for(auto uu : u) sum_u += uu;
    }
  }

  void SetParameter(ParameterType const& par) {
    if(par.size() != 5) {
      throw std::invalid_argument(
          "The par vector should be of length 5 with \
      elements corresponding to g0, alpha, theta, sigma and sigmae.");
    }
    if(par[1] < 0 || par[3] < 0 || par[4] < 0) {
      throw std::logic_error("The parameters alpha, sigma and sigmae should be non-negative.");
    }
    this->g0 = par[0];
    this->alpha = par[1];
    this->theta = par[2];
    this->sigma = par[3];
    this->sigmae = par[4];
    this->e2alphaT = exp(-2*alpha*T);
  }

  NodeStateType StateAtRoot() const {
    vec res = BaseType::StateAtRoot();
    res.push_back(alpha*sum_u);
    return res;
  }

  inline void InitNode(uint i) {
    ThreePointUnivariate<TreeType>::InitNode(i);
    if(i < this->ref_tree_.num_nodes() - 1) {
      // if an internal node or a tip transform the branch length leading to this tip
      uint iParent = this->ref_tree_.FindIdOfParent(i);
      // tTransf[i] =
      //   sigma*sigma/(2*alpha) * ( (1-exp(-2*alpha*h[i]))*exp(-2*alpha*(T-h[i])) -
      //     (1-exp(-2*alpha*h[iParent]))*exp(-2*alpha*(T-h[iParent])) );

      double ealphahi = exp(alpha*h[i]);

      this->tTransf[i] = sigma*sigma/(2*alpha) *
        (e2alphaT*(ealphahi*ealphahi - exp(2*alpha*h[iParent])));

      if(i < this->ref_tree_.num_tips()) {
        //double mu = exp(-alpha*h[i])*g0 + (1-exp(-alpha*h[i]))*theta;
        double mu = g0 / ealphahi + (1 - 1/ealphahi)*theta;
        double ealphaui = exp(alpha*u[i]);
        this->X[i] = this->Y[i] = (z[i] - mu)/ealphaui;
        this->tTransf[i] += sigmae*sigmae / (ealphaui*ealphaui);
      }
    }
  }
};

typedef ParallelPruningUse<
  ThreePointPOUMM<ParallelPruningTree<uint, double>> > ParallelPruningThreePointPOUMM;

}

#endif // ThreePointPOUMM_H_
