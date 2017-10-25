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
ppa::ParallelPruningTree<uint, double>* create_ParallelPruningTree(Rcpp::List const& tree) {
  arma::umat branches = tree["edge"];
  ppa::uvec br_0 = arma::conv_to<ppa::uvec>::from(branches.col(0));
  ppa::uvec br_1 = arma::conv_to<ppa::uvec>::from(branches.col(1));
  ppa::vec t = Rcpp::as<ppa::vec>(tree["edge.length"]);
  return new ppa::ParallelPruningTree<uint, double>(br_0, br_1, t);
}

ppa::POUMM_abc<uint>* create_POUMM_abc(Rcpp::List const& tree, ppa::vec const& z, ppa::vec const& se) {
  arma::umat branches = tree["edge"];
  ppa::uvec br_0 = arma::conv_to<ppa::uvec>::from(branches.col(0));
  ppa::uvec br_1 = arma::conv_to<ppa::uvec>::from(branches.col(1));
  ppa::vec t = Rcpp::as<ppa::vec>(tree["edge.length"]);
  uint N = Rcpp::as<Rcpp::CharacterVector>(tree["tip.label"]).size();
  return new ppa::POUMM_abc<uint>(br_0, br_1, t, ppa::seq(1, N), z, se);
}

ppa::POUMM_lnDetV_Q_1d<uint>* create_POUMM_lnDetV_Q_1d(Rcpp::List const& tree, ppa::vec const& z) {
  arma::umat branches = tree["edge"];
  ppa::uvec br_0 = arma::conv_to<ppa::uvec>::from(branches.col(0));
  ppa::uvec br_1 = arma::conv_to<ppa::uvec>::from(branches.col(1));
  ppa::vec t = Rcpp::as<ppa::vec>(tree["edge.length"]);
  uint N = Rcpp::as<Rcpp::CharacterVector>(tree["tip.label"]).size();
  return new ppa::POUMM_lnDetV_Q_1d<uint>(br_0, br_1, t, ppa::seq(1, N), z);
}

RCPP_MODULE(ParallelPruningTree) {
  Rcpp::class_<ppa::Tree<uint, double>> ( "Tree" )
  .property("M", &ppa::Tree<uint, double>::get_M )
  .property("N", &ppa::Tree<uint, double>::get_N )
  .property("t", &ppa::Tree<uint, double>::get_t )
  .method( "node", &ppa::Tree<uint, double>::get_node )
  .method( "nodes", &ppa::Tree<uint, double>::get_nodes )
  .method( "id", &ppa::Tree<uint, double>::get_id )
  ;
  Rcpp::class_<ppa::ParallelPruningTree<uint, double>>( "ParallelPruningTree" )
    .derives<ppa::Tree<uint, double>> ( "Tree" )
    .method( "orderNodes", &ppa::ParallelPruningTree<uint, double>::order_nodes )
    .factory<Rcpp::List const&>( &create_ParallelPruningTree )
    .property("nLevels", &ppa::ParallelPruningTree<uint, double>::get_nLevels )
    .property("parentNode", &ppa::ParallelPruningTree<uint, double>::get_parentNode )
    .property("tipsVector", &ppa::ParallelPruningTree<uint, double>::get_tipsVector )
    .property("tipsVectorIndex", &ppa::ParallelPruningTree<uint, double>::get_tipsVectorIndex )
    .property("branchVector", &ppa::ParallelPruningTree<uint, double>::get_branchVector )
    .property("branchVectorIndex", &ppa::ParallelPruningTree<uint, double>::get_branchVectorIndex )
    .property("orderNodeIds", &ppa::ParallelPruningTree<uint, double>::get_orderNodeIds )
    .property("orderBranches", &ppa::ParallelPruningTree<uint, double>::get_orderBranches )
    .property("nodeHeights", &ppa::ParallelPruningTree<uint, double>::get_nodeHeights )
  ;
}

RCPP_MODULE(POUMM_abc) {
  Rcpp::class_<ppa::POUMM_abc<uint>>( "POUMM_abc" )
  .factory<Rcpp::List const&, ppa::vec const&, ppa::vec const&>(&create_POUMM_abc)
  .method( "do_pruning", &ppa::POUMM_abc<uint>::do_pruning )
  .method( "set_parameters", &ppa::POUMM_abc<uint>::set_parameters )
  .method( "abc", &ppa::POUMM_abc<uint>::get_abc )
  .property( "nThreads", &ppa::POUMM_abc<uint>::get_nThreads )
  .property( "bestMinChunkSizeForHybrid", &ppa::POUMM_abc<uint>::get_bestMinChunkSizeForHybrid )
  .property( "pmaBestMode", &ppa::POUMM_abc<uint>::get_pmaBestMode )
  .property( "pmaTuningDurations", &ppa::POUMM_abc<uint>::get_pmaTuningDurations )
  .property( "pmaTuningMinChunkSizes", &ppa::POUMM_abc<uint>::get_pmaTuningMinChunkSizes )
  .field( "a", &ppa::POUMM_abc<uint>::a )
  .field( "b", &ppa::POUMM_abc<uint>::b )
  .field( "c", &ppa::POUMM_abc<uint>::c )
  .field( "se", &ppa::POUMM_abc<uint>::se )
  .field( "z", &ppa::POUMM_abc<uint>::z )
  .field( "sum_se2_sigmae2", &ppa::POUMM_abc<uint>::sum_se2_sigmae2 )
  .field( "talpha", &ppa::POUMM_abc<uint>::talpha )
  .field( "etalpha", &ppa::POUMM_abc<uint>::etalpha )
  .field( "gutalphasigma2", &ppa::POUMM_abc<uint>::gutalphasigma2 )
  .field( "e2talpha", &ppa::POUMM_abc<uint>::e2talpha )
  .field( "fe2talpha", &ppa::POUMM_abc<uint>::fe2talpha )
  .field( "alpha", &ppa::POUMM_abc<uint>::alpha )
  .field( "sigma", &ppa::POUMM_abc<uint>::sigma )
  .field( "sigmae", &ppa::POUMM_abc<uint>::sigmae )
  .field( "theta", &ppa::POUMM_abc<uint>::theta )
  ;
}

RCPP_MODULE(POUMM_lnDetV_Q_1d) {
  Rcpp::class_<ppa::ThreePointV_lnDetV_Q_1d<uint>>( "ThreePointV_lnDetV_Q_1d" )
  .property( "N", &ppa::ThreePointV_lnDetV_Q_1d<uint>::get_N )
  .property( "lnDetV", &ppa::ThreePointV_lnDetV_Q_1d<uint>::get_lnDetV )
  .property( "Q", &ppa::ThreePointV_lnDetV_Q_1d<uint>::get_Q )
  .field( "Q_all", &ppa::ThreePointV_lnDetV_Q_1d<uint>::Q )
  .field( "lnDetV_all", &ppa::ThreePointV_lnDetV_Q_1d<uint>::lnDetV )
  .field( "p", &ppa::ThreePointV_lnDetV_Q_1d<uint>::p )
  .field( "X", &ppa::ThreePointV_lnDetV_Q_1d<uint>::X )
  .field( "Y", &ppa::ThreePointV_lnDetV_Q_1d<uint>::Y )
  .field( "tTransf", &ppa::ThreePointV_lnDetV_Q_1d<uint>::tTransf )
  .field( "hat_mu_Y", &ppa::ThreePointV_lnDetV_Q_1d<uint>::hat_mu_Y )
  .field( "tilde_mu_X_prime", &ppa::ThreePointV_lnDetV_Q_1d<uint>::tilde_mu_X_prime )
  ;
  Rcpp::class_<ppa::POUMM_lnDetV_Q_1d<uint>>( "POUMM_lnDetV_Q_1d" )
    .derives<ppa::ThreePointV_lnDetV_Q_1d<uint>>("ThreePointV_lnDetV_Q_1d")
    .factory<Rcpp::List const&, ppa::vec const&>(&create_POUMM_lnDetV_Q_1d)
    .method( "set_parameters", &ppa::POUMM_lnDetV_Q_1d<uint>::set_parameters )
    .method( "do_pruning", &ppa::POUMM_lnDetV_Q_1d<uint>::do_pruning )
    .property( "lnDetD", &ppa::POUMM_lnDetV_Q_1d<uint>::get_lnDetD )
    .field( "h", &ppa::POUMM_lnDetV_Q_1d<uint>::h )
    .field( "u", &ppa::POUMM_lnDetV_Q_1d<uint>::u )
    .field( "alpha", &ppa::POUMM_lnDetV_Q_1d<uint>::alpha )
    .field( "sigma", &ppa::POUMM_lnDetV_Q_1d<uint>::sigma )
    .field( "sigmae", &ppa::POUMM_lnDetV_Q_1d<uint>::sigmae )
    .field( "g0", &ppa::POUMM_lnDetV_Q_1d<uint>::g0 )
    .field( "theta", &ppa::POUMM_lnDetV_Q_1d<uint>::theta )
    .field( "T", &ppa::POUMM_lnDetV_Q_1d<uint>::T )
    .field( "e2alphaT", &ppa::POUMM_lnDetV_Q_1d<uint>::e2alphaT )
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
