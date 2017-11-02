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

ppa::ParallelPruningTree<std::string, double>* create_ParallelPruningTree(
    std::vector<std::string> const& br_0,
    std::vector<std::string> const& br_1,
    std::vector<double> const& t) {

  return new ppa::ParallelPruningTree<std::string, double>(br_0, br_1, t);
}

ppa::ParallelPruningTreeFindChildren<uint, double>* create_ParallelPruningTreeFindChildren(Rcpp::List const& tree) {
  arma::umat branches = tree["edge"];
  ppa::uvec br_0 = arma::conv_to<ppa::uvec>::from(branches.col(0));
  ppa::uvec br_1 = arma::conv_to<ppa::uvec>::from(branches.col(1));
  ppa::vec t = Rcpp::as<ppa::vec>(tree["edge.length"]);
  return new ppa::ParallelPruningTreeFindChildren<uint, double>(br_0, br_1, t);
}

ppa::POUMM_abc<uint>* create_POUMM_abc(Rcpp::List const& tree, ppa::vec const& z, ppa::vec const& se) {
  arma::umat branches = tree["edge"];
  ppa::uvec br_0 = arma::conv_to<ppa::uvec>::from(branches.col(0));
  ppa::uvec br_1 = arma::conv_to<ppa::uvec>::from(branches.col(1));
  ppa::vec t = Rcpp::as<ppa::vec>(tree["edge.length"]);
  uint num_tips = Rcpp::as<Rcpp::CharacterVector>(tree["tip.label"]).size();
  return new ppa::POUMM_abc<uint>(br_0, br_1, t, ppa::Seq(1, num_tips), z, se);
}

ppa::POUMM_lnDetV_Q_1d<uint>* create_POUMM_lnDetV_Q_1d(Rcpp::List const& tree, ppa::vec const& z) {
  arma::umat branches = tree["edge"];
  ppa::uvec br_0 = arma::conv_to<ppa::uvec>::from(branches.col(0));
  ppa::uvec br_1 = arma::conv_to<ppa::uvec>::from(branches.col(1));
  ppa::vec t = Rcpp::as<ppa::vec>(tree["edge.length"]);
  uint num_tips = Rcpp::as<Rcpp::CharacterVector>(tree["tip.label"]).size();
  return new ppa::POUMM_lnDetV_Q_1d<uint>(br_0, br_1, t, ppa::Seq(1, num_tips), z);
}

RCPP_MODULE(ParallelPruningTreeStringNodes) {
  Rcpp::class_<ppa::Tree<std::string, double>> ( "TreeStringNodes" )
  .property("num_nodes", &ppa::Tree<std::string, double>::num_nodes )
  .property("num_tips", &ppa::Tree<std::string, double>::num_tips )
  .method("LengthOfBranch", &ppa::Tree<std::string, double>::LengthOfBranch )
  .method("FindNodeWithId", &ppa::Tree<std::string, double>::FindNodeWithId )
  .method("FindIdOfNode", &ppa::Tree<std::string, double>::FindIdOfNode )
  .method("FindIdOfParent", &ppa::Tree<std::string, double>::FindIdOfParent )
  ;
  Rcpp::class_<ppa::ParallelPruningTree<std::string, double>>( "ParallelPruningTreeStringNodes" )
    .derives<ppa::Tree<uint, double>> ( "TreeStringNodes" )
    .factory<std::vector<std::string> const&, std::vector<std::string> const&,std::vector<double> const&>( &create_ParallelPruningTree )
    .method("OrderNodes", &ppa::ParallelPruningTree<std::string, double>::OrderNodes )
  //.method("RangeIdPrune", &ppa::ParallelPruningTree<std::string, double>::RangeIdPrune )
  //.method("RangeIdUpdateParent", &ppa::ParallelPruningTree<std::string, double>::RangeIdUpdateParent )
    .method("CalculateHeights", &ppa::ParallelPruningTree<std::string, double>::CalculateHeights )
    .property("num_levels", &ppa::ParallelPruningTree<std::string, double>::num_levels )
    .property("num_parallel_ranges_prune", &ppa::ParallelPruningTree<std::string, double>::num_parallel_ranges_prune )
    .property("ranges_id_visit", &ppa::ParallelPruningTree<std::string, double>::ranges_id_visit )
    .property("ranges_id_prune", &ppa::ParallelPruningTree<std::string, double>::ranges_id_prune )
  ;
}

RCPP_MODULE(ParallelPruningTree) {
  Rcpp::class_<ppa::Tree<uint, double>> ( "Tree" )
  .property("num_nodes", &ppa::Tree<uint, double>::num_nodes )
  .property("num_tips", &ppa::Tree<uint, double>::num_tips )
  .method("LengthOfBranch", &ppa::Tree<uint, double>::LengthOfBranch )
  .method("FindNodeWithId", &ppa::Tree<uint, double>::FindNodeWithId )
  .method("FindIdOfNode", &ppa::Tree<uint, double>::FindIdOfNode )
  .method("FindIdOfParent", &ppa::Tree<uint, double>::FindIdOfParent )
  ;
  Rcpp::class_<ppa::ParallelPruningTree<uint, double>>( "ParallelPruningTree" )
    .derives<ppa::Tree<uint, double>> ( "Tree" )
    .factory<Rcpp::List const&>( &create_ParallelPruningTree )
    .method("OrderNodes", &ppa::ParallelPruningTree<uint, double>::OrderNodes )
    //.method("RangeIdPrune", &ppa::ParallelPruningTree<uint, double>::RangeIdPrune )
    //.method("RangeIdUpdateParent", &ppa::ParallelPruningTree<uint, double>::RangeIdUpdateParent )
    .method("CalculateHeights", &ppa::ParallelPruningTree<uint, double>::CalculateHeights )
    .property("num_levels", &ppa::ParallelPruningTree<uint, double>::num_levels )
    .property("ranges_id_visit", &ppa::ParallelPruningTree<uint, double>::ranges_id_visit )
    .property("ranges_id_prune", &ppa::ParallelPruningTree<uint, double>::ranges_id_prune )
  ;
}

RCPP_MODULE(ParallelPruningTreeFindChildren) {
  Rcpp::class_<ppa::Tree<uint, double>> ( "Tree" )
  .property("num_nodes", &ppa::Tree<uint, double>::num_nodes )
  .property("num_tips", &ppa::Tree<uint, double>::num_tips )
  .method("LengthOfBranch", &ppa::Tree<uint, double>::LengthOfBranch )
  .method("FindNodeWithId", &ppa::Tree<uint, double>::FindNodeWithId )
  .method("FindIdOfNode", &ppa::Tree<uint, double>::FindIdOfNode )
  .method("FindIdOfParent", &ppa::Tree<uint, double>::FindIdOfParent )
  ;
  Rcpp::class_<ppa::ParallelPruningTree<uint, double>>( "ParallelPruningTree" )
    .derives<ppa::Tree<uint, double>> ( "Tree" )
    .factory<Rcpp::List const&>( &create_ParallelPruningTree )
    .method("OrderNodes", &ppa::ParallelPruningTree<uint, double>::OrderNodes )
  //.method("RangeIdPrune", &ppa::ParallelPruningTree<uint, double>::RangeIdPrune )
  //.method("RangeIdUpdateParent", &ppa::ParallelPruningTree<uint, double>::RangeIdUpdateParent )
    .method("CalculateHeights", &ppa::ParallelPruningTree<uint, double>::CalculateHeights )
    .property("num_levels", &ppa::ParallelPruningTree<uint, double>::num_levels )
    .property("ranges_id_visit", &ppa::ParallelPruningTree<uint, double>::ranges_id_visit )
    .property("ranges_id_prune", &ppa::ParallelPruningTree<uint, double>::ranges_id_prune )
  ;
  Rcpp::class_<ppa::ParallelPruningTreeFindChildren<uint, double>>( "ParallelPruningTreeFindChildren" )
    .derives<ppa::ParallelPruningTree<uint, double>>( "ParallelPruningTree" )
    .factory<Rcpp::List const&>( &create_ParallelPruningTreeFindChildren )
    .method( "FindChildren", &ppa::ParallelPruningTreeFindChildren<uint, double>::FindChildren )
  ;
  }

RCPP_MODULE(POUMM_abc) {
  Rcpp::class_<ppa::POUMM_abc<uint>>( "POUMM_abc" )
  .factory<Rcpp::List const&, ppa::vec const&, ppa::vec const&>(&create_POUMM_abc)
  .method( "DoPruning", &ppa::POUMM_abc<uint>::DoPruning )
  .method( "set_parameters", &ppa::POUMM_abc<uint>::set_parameters )
  .method( "abc", &ppa::POUMM_abc<uint>::get_abc )
  .property( "ModeAuto", &ppa::POUMM_abc<uint>::ModeAuto )
  .property( "IndexMinSizeChunkVisit", &ppa::POUMM_abc<uint>::IndexMinSizeChunkVisit )
  .property( "IndexMinSizeChunkPrune", &ppa::POUMM_abc<uint>::IndexMinSizeChunkPrune )
  .property( "IsTuning", &ppa::POUMM_abc<uint>::IsTuning )
  .property( "num_threads", &ppa::POUMM_abc<uint>::num_threads )
  .property( "min_size_chunk_visit", &ppa::POUMM_abc<uint>::min_size_chunk_visit )
  .property( "min_size_chunk_prune", &ppa::POUMM_abc<uint>::min_size_chunk_prune )
  .property( "durations_tuning", &ppa::POUMM_abc<uint>::durations_tuning )
  .property( "fastest_step_tuning", &ppa::POUMM_abc<uint>::fastest_step_tuning )
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
  .property( "num_tips", &ppa::ThreePointV_lnDetV_Q_1d<uint>::num_tips )
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
    .method( "DoPruning", &ppa::POUMM_lnDetV_Q_1d<uint>::DoPruning )
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
