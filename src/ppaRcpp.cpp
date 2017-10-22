#include <RcppArmadillo.h>
#include <R_ext/Rdynload.h>

#include "POUMM_abc.h"
#include "POUMM_lnDetV_Q_1d.h"

// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppArmadillo)]]

// BEGIN: Needed for r-devel (R 3.4)
void R_init_ParallelPruning(DllInfo *info) {
  /* Register routines, allocate resources. */
  R_registerRoutines(info, NULL, NULL, NULL, NULL);
  R_useDynamicSymbols(info, TRUE);
}

void R_unload_ParallelPruning(DllInfo *info) {
  /* Release resources. */
}
// END Needed for r-devel (R 3.4)


// Construct a ppa::Tree from a ape::phylo object (represented as Rcpp::List)
ppa::Tree phylo(Rcpp::List const& tree) {
  if(!tree.inherits("phylo")) {
    Rcpp::stop("Input must be a phylo object.");
  }
  // 0-based indices
  arma::umat branches = tree["edge"];
  arma::uvec br_0 = branches.col(0) - 1;
  arma::uvec br_1 = branches.col(1) - 1;

  return ppa::Tree(Rcpp::as<Rcpp::CharacterVector>(tree["tip.label"]).size(),
                   ppa::uvec(br_0.begin(), br_0.end()),
                   ppa::uvec(br_1.begin(), br_1.end()),
                   Rcpp::as<ppa::vec>(tree["edge.length"]));
}

// Construct a ppa::Tree from number of tips, edge-matrix and edge.length vector
ppa::Tree phylo(uint N, arma::umat const& branches, arma::vec const& edgeLength) {
  arma::uvec br_0 = branches.col(0) - 1;
  arma::uvec br_1 = branches.col(1) - 1;
  return ppa::Tree(N,
                   ppa::uvec(br_0.begin(), br_0.end()),
                   ppa::uvec(br_1.begin(), br_1.end()),
                   ppa::vec(edgeLength.begin(), edgeLength.end()));
}

ppa::ParallelPruningTree* create_ParallelPruningTree(Rcpp::List const& tree) {
  return new ppa::ParallelPruningTree(phylo(tree));
}

ppa::POUMM_abc* create_POUMM_abc(Rcpp::List const& tree, ppa::vec const& z, ppa::vec const& se) {
  return new ppa::POUMM_abc(phylo(tree), z, se);
}

ppa::POUMM_lnDetV_Q_1d* create_POUMM_lnDetV_Q_1d(Rcpp::List const& tree, ppa::vec const& z) {
  return new ppa::POUMM_lnDetV_Q_1d(phylo(tree), z);
}

RCPP_MODULE(ParallelPruningTree) {
  Rcpp::class_<ppa::ParallelPruningTree>( "ParallelPruningTree" )
  .factory<Rcpp::List const&>( &create_ParallelPruningTree )
  .property("nLevels", &ppa::ParallelPruningTree::get_nLevels )
  .property("M", &ppa::ParallelPruningTree::get_M )
  .property("N", &ppa::ParallelPruningTree::get_N )
  .property("t", &ppa::ParallelPruningTree::get_t )
  .property("parentNode", &ppa::ParallelPruningTree::get_parentNode )
  .property("tipsVector", &ppa::ParallelPruningTree::get_tipsVector )
  .property("tipsVectorIndex", &ppa::ParallelPruningTree::get_tipsVectorIndex )
  .property("branchVector", &ppa::ParallelPruningTree::get_branchVector )
  .property("branchVectorIndex", &ppa::ParallelPruningTree::get_branchVectorIndex )
  .property("orderNodes", &ppa::ParallelPruningTree::get_orderNodes )
  .property("orderBranches", &ppa::ParallelPruningTree::get_orderBranches )
  .property("nodeHeights", &ppa::ParallelPruningTree::get_nodeHeights )
  ;
}

RCPP_MODULE(POUMM_abc) {
  Rcpp::class_<ppa::POUMM_abc>( "POUMM_abc" )
    .factory<Rcpp::List const&, ppa::vec const&, ppa::vec const&>(&create_POUMM_abc)
    .method( "do_pruning", &ppa::POUMM_abc::do_pruning )
    .method( "set_parameters", &ppa::POUMM_abc::set_parameters )
    .method( "abc", &ppa::POUMM_abc::get_abc )
    .field( "a", &ppa::POUMM_abc::a )
    .field( "b", &ppa::POUMM_abc::b )
    .field( "c", &ppa::POUMM_abc::c )
    .field( "se", &ppa::POUMM_abc::se )
    .field( "z", &ppa::POUMM_abc::z )
    .field( "sum_se2_sigmae2", &ppa::POUMM_abc::sum_se2_sigmae2 )
    .field( "talpha", &ppa::POUMM_abc::talpha )
    .field( "etalpha", &ppa::POUMM_abc::etalpha )
    .field( "gutalphasigma2", &ppa::POUMM_abc::gutalphasigma2 )
    .field( "e2talpha", &ppa::POUMM_abc::e2talpha )
    .field( "fe2talpha", &ppa::POUMM_abc::fe2talpha )
    .field( "alpha", &ppa::POUMM_abc::alpha )
    .field( "sigma", &ppa::POUMM_abc::sigma )
    .field( "sigmae", &ppa::POUMM_abc::sigmae )
    .field( "theta", &ppa::POUMM_abc::theta )
  ;
}

RCPP_MODULE(POUMM_lnDetV_Q_1d) {
  Rcpp::class_<ppa::ThreePointV_lnDetV_Q_1d>( "ThreePointV_lnDetV_Q_1d" )
    .property( "N", &ppa::ThreePointV_lnDetV_Q_1d::get_N )
    .property( "lnDetV", &ppa::ThreePointV_lnDetV_Q_1d::get_lnDetV )
    .property( "Q", &ppa::ThreePointV_lnDetV_Q_1d::get_Q )
    .field( "Q_all", &ppa::ThreePointV_lnDetV_Q_1d::Q )
    .field( "lnDetV_all", &ppa::ThreePointV_lnDetV_Q_1d::lnDetV )
    .field( "p", &ppa::ThreePointV_lnDetV_Q_1d::p )
    .field( "X", &ppa::ThreePointV_lnDetV_Q_1d::X )
    .field( "Y", &ppa::ThreePointV_lnDetV_Q_1d::Y )
    .field( "tTransf", &ppa::ThreePointV_lnDetV_Q_1d::tTransf )
    .field( "hat_mu_Y", &ppa::ThreePointV_lnDetV_Q_1d::hat_mu_Y )
    .field( "tilde_mu_X_prime", &ppa::ThreePointV_lnDetV_Q_1d::tilde_mu_X_prime )
  ;
  Rcpp::class_<ppa::POUMM_lnDetV_Q_1d>( "POUMM_lnDetV_Q_1d" )
    .derives<ppa::ThreePointV_lnDetV_Q_1d>("ThreePointV_lnDetV_Q_1d")
    .factory<Rcpp::List const&, ppa::vec const&>(&create_POUMM_lnDetV_Q_1d)
    .method( "set_parameters", &ppa::POUMM_lnDetV_Q_1d::set_parameters )
    .method( "do_pruning", &ppa::POUMM_lnDetV_Q_1d::do_pruning )
    .property( "lnDetD", &ppa::POUMM_lnDetV_Q_1d::get_lnDetD )
    .field( "h", &ppa::POUMM_lnDetV_Q_1d::h )
    .field( "u", &ppa::POUMM_lnDetV_Q_1d::u )
    .field( "alpha", &ppa::POUMM_lnDetV_Q_1d::alpha )
    .field( "sigma", &ppa::POUMM_lnDetV_Q_1d::sigma )
    .field( "sigmae", &ppa::POUMM_lnDetV_Q_1d::sigmae )
    .field( "g0", &ppa::POUMM_lnDetV_Q_1d::g0 )
    .field( "theta", &ppa::POUMM_lnDetV_Q_1d::theta )
    .field( "T", &ppa::POUMM_lnDetV_Q_1d::T )
    .field( "e2alphaT", &ppa::POUMM_lnDetV_Q_1d::e2alphaT )
  ;
}

// RCPP_MODULE(PBM_lnDetV_Q_1d) {
//   Rcpp::class_<ppa::ParallelPruningAlgorithm>( "ParallelPruningAlgorithm" )
//   .method( "do_pruning", &ppa::ParallelPruningAlgorithm::do_pruning )
//   .property("nLevels", &ppa::ParallelPruningAlgorithm::get_nLevels )
//   .property("M", &ppa::ParallelPruningAlgorithm::get_M )
//   .property("N", &ppa::ParallelPruningAlgorithm::get_N )
//   .property("t", &ppa::ParallelPruningAlgorithm::get_t )
//   .property("parentNode", &ppa::ParallelPruningAlgorithm::get_parentNode )
//   .property("tipsVector", &ppa::ParallelPruningAlgorithm::get_tipsVector )
//   .property("tipsVectorIndex", &ppa::ParallelPruningAlgorithm::get_tipsVectorIndex )
//   .property("branchVector", &ppa::ParallelPruningAlgorithm::get_branchVector )
//   .property("branchVectorIndex", &ppa::ParallelPruningAlgorithm::get_branchVectorIndex )
//   .property("orderNodes", &ppa::ParallelPruningAlgorithm::get_orderNodes )
//   .property("orderBranches", &ppa::ParallelPruningAlgorithm::get_orderBranches )
//   .property("nodeHeights", &ppa::ParallelPruningAlgorithm::get_nodeHeights )
//   ;
//   Rcpp::class_<ppa::ThreePointV_lnDetV_Q_1d>( "ThreePointV_lnDetV_Q_1d" )
//     .derives<ppa::ParallelPruningAlgorithm>("ParallelPruningAlgorithm")
//     .constructor<ppa::Tree const&>()
//     .method( "lnDetV", &ppa::ThreePointV_lnDetV_Q_1d::get_lnDetV )
//     .method( "Q", &ppa::ThreePointV_lnDetV_Q_1d::get_Q )
//     .field( "Q_all", &ppa::ThreePointV_lnDetV_Q_1d::Q )
//     .field( "lnDetV_all", &ppa::ThreePointV_lnDetV_Q_1d::lnDetV )
//     .field( "p", &ppa::ThreePointV_lnDetV_Q_1d::p )
//     .field( "X", &ppa::ThreePointV_lnDetV_Q_1d::X )
//     .field( "Y", &ppa::ThreePointV_lnDetV_Q_1d::Y )
//     .field( "tTransf", &ppa::ThreePointV_lnDetV_Q_1d::tTransf )
//     .field( "hat_mu_Y", &ppa::ThreePointV_lnDetV_Q_1d::hat_mu_Y )
//     .field( "tilde_mu_X_prime", &ppa::ThreePointV_lnDetV_Q_1d::tilde_mu_X_prime )
//   ;
//   Rcpp::class_<ppa::PBM_lnDetV_Q_1d>( "PBM_lnDetV_Q_1d" )
//     .derives<ppa::ThreePointV_lnDetV_Q_1d>("ThreePointV_lnDetV_Q_1d")
//     .constructor<ppa::Tree const&, ppa::vec const&>( )
//     .method( "set_parameters", &ppa::PBM_lnDetV_Q_1d::set_parameters )
//   ;
//
// }
//
// RCPP_MODULE(POUMM_lnDetV_Q) {
//   Rcpp::class_<ppa::ParallelPruningAlgorithm>( "ParallelPruningAlgorithm" )
//   .method( "do_pruning", &ppa::ParallelPruningAlgorithm::do_pruning )
//   .property("nLevels", &ppa::ParallelPruningAlgorithm::get_nLevels )
//   .property("M", &ppa::ParallelPruningAlgorithm::get_M )
//   .property("N", &ppa::ParallelPruningAlgorithm::get_N )
//   .property("t", &ppa::ParallelPruningAlgorithm::get_t )
//   .property("parentNode", &ppa::ParallelPruningAlgorithm::get_parentNode )
//   .property("tipsVector", &ppa::ParallelPruningAlgorithm::get_tipsVector )
//   .property("tipsVectorIndex", &ppa::ParallelPruningAlgorithm::get_tipsVectorIndex )
//   .property("branchVector", &ppa::ParallelPruningAlgorithm::get_branchVector )
//   .property("branchVectorIndex", &ppa::ParallelPruningAlgorithm::get_branchVectorIndex )
//   .property("orderNodes", &ppa::ParallelPruningAlgorithm::get_orderNodes )
//   .property("orderBranches", &ppa::ParallelPruningAlgorithm::get_orderBranches )
//   .property("nodeHeights", &ppa::ParallelPruningAlgorithm::get_nodeHeights )
//   ;
//   Rcpp::class_<ppa::ThreePointV_lnDetV_Q>( "ThreePointV_lnDetV_Q" )
//     .derives<ppa::ParallelPruningAlgorithm>("ParallelPruningAlgorithm")
//     .constructor<ppa::Tree const&>()
//     .method( "lnDetV", &ppa::ThreePointV_lnDetV_Q::get_lnDetV )
//     .method( "Q", &ppa::ThreePointV_lnDetV_Q::get_Q )
//     .field( "Q_all", &ppa::ThreePointV_lnDetV_Q::Q )
//     .field( "lnDetV_all", &ppa::ThreePointV_lnDetV_Q::lnDetV )
//     .field( "p", &ppa::ThreePointV_lnDetV_Q::p )
//     .field( "X", &ppa::ThreePointV_lnDetV_Q::X )
//     .field( "Y", &ppa::ThreePointV_lnDetV_Q::Y )
//     .field( "tTransf", &ppa::ThreePointV_lnDetV_Q::tTransf )
//     .field( "hat_mu_Y", &ppa::ThreePointV_lnDetV_Q::hat_mu_Y )
//     .field( "tilde_mu_X_prime", &ppa::ThreePointV_lnDetV_Q::tilde_mu_X_prime )
//   ;
//   Rcpp::class_<ppa::POUMM_lnDetV_Q>( "POUMM_lnDetV_Q" )
//     .derives<ppa::ThreePointV_lnDetV_Q>("ThreePointV_lnDetV_Q")
//     .constructor<ppa::Tree const&, ppa::vec const&>( )
//     .method( "set_parameters", &ppa::POUMM_lnDetV_Q::set_parameters )
//     .property( "lnDetD", &ppa::POUMM_lnDetV_Q::get_lnDetD )
//   ;
//
// }
//
