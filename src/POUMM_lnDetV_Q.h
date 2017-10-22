#ifndef POUMM_lnDetV_Q_H_
#define POUMM_lnDetV_Q_H_

namespace ppa {
#include "ThreePointV_lnDetV_Q.h"

class POUMM_lnDetV_Q: public ThreePointV_lnDetV_Q {
public:
  // univariate trait vector
  arma::vec z;
  ppa::vec h, u;
  double g0, alpha, theta, sigma, sigmae, e2alphaT;
  double T; // tree height


  POUMM_lnDetV_Q(Rcpp::List const& tree, arma::vec const& z):
    ThreePointV_lnDetV_Q(tree), z(z) {

    if(z.n_elem != N ) {
      Rcpp::stop("The trait vector must have N elements (N is the number of tips).");
    } else {
      arma::vec z2 = z;
      for(int i = 0; i < N; ++i) {
        this->z[i] = z2[orderNodes[i]];
      }
      arma::mat X(N, 1);
      set_X_and_Y(X, X);
    }

    this->h = get_nodeHeights();

    this->T = *std::max_element(h.begin(), h.begin()+N);
    this->u = ppa::vec(N);
    for(int i = 0; i < N; i++) u[i] = T - h[i];
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
    uint iParent = parentNode[i];
    // tTransf[i] =
    //   sigma*sigma/(2*alpha) * ( (1-exp(-2*alpha*h[i]))*exp(-2*alpha*(T-h[i])) -
    //     (1-exp(-2*alpha*h[iParent]))*exp(-2*alpha*(T-h[iParent])) );

    double ealphahi = exp(alpha*h[i]);

    tTransf[i] = sigma*sigma/(2*alpha) *
      (e2alphaT*(ealphahi*ealphahi - exp(2*alpha*h[iParent])));

    if(i < N) {
      //double mu = exp(-alpha*h[i])*g0 + (1-exp(-alpha*h[i]))*theta;
      double mu = g0 / ealphahi + (1 - 1/ealphahi)*theta;

      double ealphaui = exp(alpha*u[i]);

      X(i, 0) = Y(i, 0) = (z[i] - mu)/ealphaui;

      tTransf[i] += sigmae*sigmae / (ealphaui*ealphaui);
    }
  }

  double get_lnDetD() const {
    double sum_u = 0;
    for(auto uu : u) sum_u += uu;
    return alpha*sum_u;
  }
};
};
# endif //POUMM_lnDetV_Q_H_

