#ifndef POUMM_lnDetV_Q_1D_H_
#define POUMM_lnDetV_Q_1D_H_

#include "ThreePointV_lnDetV_Q_1d.h"

namespace ppa {
class POUMM_lnDetV_Q_1d: public ThreePointV_lnDetV_Q_1d {
  ParallelPruningAlgorithm<POUMM_lnDetV_Q_1d> ppalgorithm;
public:
  // univariate trait vector
  ppa::vec z;
  ppa::vec h, u;
  double g0, alpha, theta, sigma, sigmae, e2alphaT;
  double T; // tree height


  POUMM_lnDetV_Q_1d(Tree const& tree, vec const& z):
    ThreePointV_lnDetV_Q_1d(tree), ppalgorithm(&this->pptree, this), z(z) {

    if(z.size() != pptree.N ) {
      throw std::invalid_argument("The trait vector must have N elements (N is the number of tips).");
    } else {
      vec z2 = z;
      for(int i = 0; i < pptree.N; ++i) {
        this->z[i] = z2[pptree.orderNodes[i]];
      }
      vec X(pptree.N);
      set_X_and_Y(X, X);
    }

    this->h = pptree.get_nodeHeights();

    this->T = *std::max_element(h.begin(), h.begin()+pptree.N);
    this->u = ppa::vec(pptree.N);
    for(int i = 0; i < pptree.N; i++) u[i] = T - h[i];
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
    uint iParent = pptree.parentNode[i];
    // tTransf[i] =
    //   sigma*sigma/(2*alpha) * ( (1-exp(-2*alpha*h[i]))*exp(-2*alpha*(T-h[i])) -
    //     (1-exp(-2*alpha*h[iParent]))*exp(-2*alpha*(T-h[iParent])) );

    double ealphahi = exp(alpha*h[i]);

    tTransf[i] = sigma*sigma/(2*alpha) *
      (e2alphaT*(ealphahi*ealphahi - exp(2*alpha*h[iParent])));

    if(i < pptree.N) {
      //double mu = exp(-alpha*h[i])*g0 + (1-exp(-alpha*h[i]))*theta;
      double mu = g0 / ealphahi + (1 - 1/ealphahi)*theta;

      double ealphaui = exp(alpha*u[i]);

      X[i] = Y[i] = (z[i] - mu)/ealphaui;

      tTransf[i] += sigmae*sigmae / (ealphaui*ealphaui);
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
