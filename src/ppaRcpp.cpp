#include <RcppArmadillo.h>
#include <R_ext/Rdynload.h>

#include "AbcPOUMM.h"
#include "ThreePointPOUMM.h"

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

ppa::Tree<uint, double>* CreateTree(Rcpp::List const& tree) {
  arma::umat branches = tree["edge"];
  ppa::uvec br_0 = arma::conv_to<ppa::uvec>::from(branches.col(0));
  ppa::uvec br_1 = arma::conv_to<ppa::uvec>::from(branches.col(1));
  ppa::vec t = Rcpp::as<ppa::vec>(tree["edge.length"]);
  return new ppa::Tree<uint, double>(br_0, br_1, t);
}

RCPP_MODULE(Tree1) {
  Rcpp::class_<ppa::Tree<uint, double>> ( "Tree1" )
  .factory<Rcpp::List const&>( &CreateTree )
  .property("num_nodes", &ppa::Tree<uint, double>::num_nodes )
  .property("num_tips", &ppa::Tree<uint, double>::num_tips )
  .method("LengthOfBranch", &ppa::Tree<uint, double>::LengthOfBranch )
  .method("FindNodeWithId", &ppa::Tree<uint, double>::FindNodeWithId )
  .method("FindIdOfNode", &ppa::Tree<uint, double>::FindIdOfNode )
  .method("FindIdOfParent", &ppa::Tree<uint, double>::FindIdOfParent )
  .method( "FindChildren", &ppa::Tree<uint, double>::FindChildren )
  ;
}

ppa::PruningTree<uint, double>* CreatePruningTree(Rcpp::List const& tree) {
  arma::umat branches = tree["edge"];
  ppa::uvec br_0 = arma::conv_to<ppa::uvec>::from(branches.col(0));
  ppa::uvec br_1 = arma::conv_to<ppa::uvec>::from(branches.col(1));
  ppa::vec t = Rcpp::as<ppa::vec>(tree["edge.length"]);
  return new ppa::PruningTree<uint, double>(br_0, br_1, t);
}


RCPP_MODULE(PruningTree1) {
  Rcpp::class_<ppa::Tree<uint, double>> ( "Tree1" )
  .factory<Rcpp::List const&>( &CreateTree )
  .property("num_nodes", &ppa::Tree<uint, double>::num_nodes )
  .property("num_tips", &ppa::Tree<uint, double>::num_tips )
  .method("LengthOfBranch", &ppa::Tree<uint, double>::LengthOfBranch )
  .method("FindNodeWithId", &ppa::Tree<uint, double>::FindNodeWithId )
  .method("FindIdOfNode", &ppa::Tree<uint, double>::FindIdOfNode )
  .method("FindIdOfParent", &ppa::Tree<uint, double>::FindIdOfParent )
  .method( "FindChildren", &ppa::Tree<uint, double>::FindChildren )
  .method("OrderNodes", &ppa::Tree<uint, double>::OrderNodes )
  ;
  Rcpp::class_<ppa::PruningTree<uint, double>>( "PruningTree1" )
    .derives<ppa::Tree<uint, double>> ( "Tree1" )
    .factory<Rcpp::List const&>( &CreatePruningTree )
    .method("CalculateHeights", &ppa::PruningTree<uint, double>::CalculateHeights )
    .property("num_levels", &ppa::PruningTree<uint, double>::num_levels )
    .property("ranges_id_visit", &ppa::PruningTree<uint, double>::ranges_id_visit )
    .property("ranges_id_prune", &ppa::PruningTree<uint, double>::ranges_id_prune )
  ;
}

ppa::PruningTree<std::string, double>* CreatePruningTreeStringNodes(
    std::vector<std::string> const& br_0,
    std::vector<std::string> const& br_1,
    std::vector<double> const& t) {

  return new ppa::PruningTree<std::string, double>(br_0, br_1, t);
}

RCPP_MODULE(PruningTreeStringNodes) {
  Rcpp::class_<ppa::Tree<std::string, double>> ( "TreeStringNodes" )
  .property("num_nodes", &ppa::Tree<std::string, double>::num_nodes )
  .property("num_tips", &ppa::Tree<std::string, double>::num_tips )
  .method("LengthOfBranch", &ppa::Tree<std::string, double>::LengthOfBranch )
  .method("FindNodeWithId", &ppa::Tree<std::string, double>::FindNodeWithId )
  .method("FindIdOfNode", &ppa::Tree<std::string, double>::FindIdOfNode )
  .method("FindIdOfParent", &ppa::Tree<std::string, double>::FindIdOfParent )
  .method( "FindChildren", &ppa::Tree<std::string, double>::FindChildren )
  .method("OrderNodes", &ppa::PruningTree<std::string, double>::OrderNodes )
  ;
  Rcpp::class_<ppa::PruningTree<std::string, double>>( "PruningTreeStringNodes" )
    .derives<ppa::Tree<std::string, double>> ( "TreeStringNodes" )
    .factory<std::vector<std::string> const&, std::vector<std::string> const&,std::vector<double> const&>( &CreatePruningTreeStringNodes )
    .method("CalculateHeights", &ppa::PruningTree<std::string, double>::CalculateHeights )
    .property("num_levels", &ppa::PruningTree<std::string, double>::num_levels )
    .property("num_parallel_ranges_prune", &ppa::PruningTree<std::string, double>::num_parallel_ranges_prune )
    .property("ranges_id_visit", &ppa::PruningTree<std::string, double>::ranges_id_visit )
    .property("ranges_id_prune", &ppa::PruningTree<std::string, double>::ranges_id_prune )
  ;
}

ppa::ParallelPruningThreePointPOUMM* CreateParallelPruningThreePointPOUMM(
    Rcpp::List const& tree, ppa::vec const& z, ppa::vec const& se) {
  arma::umat branches = tree["edge"];
  ppa::uvec br_0 = arma::conv_to<ppa::uvec>::from(branches.col(0));
  ppa::uvec br_1 = arma::conv_to<ppa::uvec>::from(branches.col(1));
  ppa::vec t = Rcpp::as<ppa::vec>(tree["edge.length"]);
  uint num_tips = Rcpp::as<Rcpp::CharacterVector>(tree["tip.label"]).size();
  ppa::uvec node_names = ppa::Seq(1, num_tips);
  typename ppa::ParallelPruningThreePointPOUMM::InputDataType data(node_names, z, se);
  return new ppa::ParallelPruningThreePointPOUMM(br_0, br_1, t, data);
}

RCPP_EXPOSED_CLASS_NODECL(ppa::ParallelPruningThreePointPOUMM::TreeType)
RCPP_EXPOSED_CLASS_NODECL(ppa::ParallelPruningThreePointPOUMM::PruningSpecType)
RCPP_EXPOSED_CLASS_NODECL(ppa::ParallelPruningThreePointPOUMM::PruningAlgorithmType)


RCPP_MODULE(ParallelPruningThreePointPOUMM) {
  Rcpp::class_<ppa::ParallelPruningThreePointPOUMM::TreeType::Tree> ( "Tree" )
  .property("num_nodes", &ppa::ParallelPruningThreePointPOUMM::TreeType::Tree::num_nodes )
  .property("num_tips", &ppa::ParallelPruningThreePointPOUMM::TreeType::Tree::num_tips )
  .method("LengthOfBranch", &ppa::ParallelPruningThreePointPOUMM::TreeType::Tree::LengthOfBranch )
  .method("FindNodeWithId", &ppa::ParallelPruningThreePointPOUMM::TreeType::Tree::FindNodeWithId )
  .method("FindIdOfNode", &ppa::ParallelPruningThreePointPOUMM::TreeType::Tree::FindIdOfNode )
  .method("FindIdOfParent", &ppa::ParallelPruningThreePointPOUMM::TreeType::Tree::FindIdOfParent )
  .method("OrderNodes", &ppa::ParallelPruningThreePointPOUMM::TreeType::Tree::OrderNodes )
  ;
  Rcpp::class_<ppa::ParallelPruningThreePointPOUMM::TreeType>( "PruningTree" )
    .derives<ppa::ParallelPruningThreePointPOUMM::TreeType::Tree> ( "Tree" )
    .method("CalculateHeights", &ppa::ParallelPruningThreePointPOUMM::TreeType::CalculateHeights )
    .property("num_levels", &ppa::ParallelPruningThreePointPOUMM::TreeType::num_levels )
    .property("ranges_id_visit", &ppa::ParallelPruningThreePointPOUMM::TreeType::ranges_id_visit )
    .property("ranges_id_prune", &ppa::ParallelPruningThreePointPOUMM::TreeType::ranges_id_prune )
  ;
  Rcpp::class_<ppa::ParallelPruningThreePointPOUMM::PruningSpecType::BaseType>( "ThreePointUnivariate" )
  .field( "Q", &ppa::ParallelPruningThreePointPOUMM::PruningSpecType::BaseType::Q )
  .field( "lnDetV_all", &ppa::ParallelPruningThreePointPOUMM::PruningSpecType::BaseType::lnDetV )
  .field( "p", &ppa::ParallelPruningThreePointPOUMM::PruningSpecType::BaseType::p )
  .field( "X", &ppa::ParallelPruningThreePointPOUMM::PruningSpecType::BaseType::X )
  .field( "Y", &ppa::ParallelPruningThreePointPOUMM::PruningSpecType::BaseType::Y )
  .field( "tTransf", &ppa::ParallelPruningThreePointPOUMM::PruningSpecType::BaseType::tTransf )
  .field( "hat_mu_Y", &ppa::ParallelPruningThreePointPOUMM::PruningSpecType::BaseType::hat_mu_Y )
  .field( "tilde_mu_X_prime", &ppa::ParallelPruningThreePointPOUMM::PruningSpecType::BaseType::tilde_mu_X_prime )
  ;
  Rcpp::class_<ppa::ParallelPruningThreePointPOUMM::PruningSpecType>( "ThreePointPOUMM" )
    .derives<ppa::ParallelPruningThreePointPOUMM::PruningSpecType::BaseType>("ThreePointUnivariate")
    .field( "h", &ppa::ParallelPruningThreePointPOUMM::PruningSpecType::h )
    .field( "u", &ppa::ParallelPruningThreePointPOUMM::PruningSpecType::u )
    .field( "alpha", &ppa::ParallelPruningThreePointPOUMM::PruningSpecType::alpha )
    .field( "sigma", &ppa::ParallelPruningThreePointPOUMM::PruningSpecType::sigma )
    .field( "sigmae", &ppa::ParallelPruningThreePointPOUMM::PruningSpecType::sigmae )
    .field( "g0", &ppa::ParallelPruningThreePointPOUMM::PruningSpecType::g0 )
    .field( "theta", &ppa::ParallelPruningThreePointPOUMM::PruningSpecType::theta )
    .field( "T", &ppa::ParallelPruningThreePointPOUMM::PruningSpecType::T )
    .field( "e2alphaT", &ppa::ParallelPruningThreePointPOUMM::PruningSpecType::e2alphaT )
  ;
  Rcpp::class_<ppa::ParallelPruningThreePointPOUMM>( "ParallelPruningThreePointPOUMM" )
    .factory<Rcpp::List const&, ppa::vec const&, ppa::vec const&>(&CreateParallelPruningThreePointPOUMM)
    .method( "DoPruning", &ppa::ParallelPruningThreePointPOUMM::DoPruning )
    .property( "tree", &ppa::ParallelPruningThreePointPOUMM::tree )
    .property( "spec", &ppa::ParallelPruningThreePointPOUMM::spec )
    .property( "algorithm", &ppa::ParallelPruningThreePointPOUMM::algorithm )
  ;
}

ppa::ParallelPruningAbcPOUMM* CreateParallelPruningAbcPOUMM(
    Rcpp::List const& tree, ppa::vec const& z, ppa::vec const& se) {
  arma::umat branches = tree["edge"];
  ppa::uvec br_0 = arma::conv_to<ppa::uvec>::from(branches.col(0));
  ppa::uvec br_1 = arma::conv_to<ppa::uvec>::from(branches.col(1));
  ppa::vec t = Rcpp::as<ppa::vec>(tree["edge.length"]);
  uint num_tips = Rcpp::as<Rcpp::CharacterVector>(tree["tip.label"]).size();
  ppa::uvec node_names = ppa::Seq(1, num_tips);
  typename ppa::ParallelPruningAbcPOUMM::InputDataType data(node_names, z, se);
  return new ppa::ParallelPruningAbcPOUMM(br_0, br_1, t, data);
}



// template <> SEXP Rcpp::wrap(ppa::ParallelPruningAbcPOUMM::TreeType const& tree) {
//   using namespace Rcpp;
//   typedef ppa::ParallelPruningAbcPOUMM::TreeType A;
//   Rcpp::XPtr<A> xp( &tree, true ) ; // copy and mark as finalizable
//   Function maker=Environment::Rcpp_namespace()[ "cpp_object_maker"];
//   return maker ( typeid(A).name() , xp );
// }
//
// template <> SEXP Rcpp::wrap(ppa::ParallelPruningAbcPOUMM::PruningSpecType const& spec) {
//   using namespace Rcpp;
//   typedef ppa::ParallelPruningAbcPOUMM::PruningSpecType A;
//   Rcpp::XPtr<A> xp( &spec, true ) ; // copy and mark as finalizable
//   Function maker=Environment::Rcpp_namespace()[ "cpp_object_maker"];
//   return maker ( typeid(A).name() , xp );
// }
//
// template <> SEXP Rcpp::wrap(ppa::ParallelPruningAbcPOUMM::PruningAlgorithmType & algorithm) {
//   using namespace Rcpp;
//   typedef ppa::ParallelPruningAbcPOUMM::PruningAlgorithmType A;
//   Rcpp::XPtr<A> xp( new A(algorithm), true ) ; // copy and mark as finalizable
//   Function maker=Environment::Rcpp_namespace()[ "cpp_object_maker"];
//   return maker ( typeid(A).name() , xp );
// }

RCPP_EXPOSED_CLASS_NODECL(ppa::ParallelPruningAbcPOUMM::PruningSpecType)
RCPP_EXPOSED_CLASS_NODECL(ppa::ParallelPruningAbcPOUMM::PruningAlgorithmType)

RCPP_MODULE(ParallelPruningAbcPOUMM) {
  Rcpp::class_<ppa::ParallelPruningAbcPOUMM::TreeType::Tree> ( "Tree" )
  .property("num_nodes", &ppa::ParallelPruningAbcPOUMM::TreeType::Tree::num_nodes )
  .property("num_tips", &ppa::ParallelPruningAbcPOUMM::TreeType::Tree::num_tips )
  .method("LengthOfBranch", &ppa::ParallelPruningAbcPOUMM::TreeType::Tree::LengthOfBranch )
  .method("FindNodeWithId", &ppa::ParallelPruningAbcPOUMM::TreeType::Tree::FindNodeWithId )
  .method("FindIdOfNode", &ppa::ParallelPruningAbcPOUMM::TreeType::Tree::FindIdOfNode )
  .method("FindIdOfParent", &ppa::ParallelPruningAbcPOUMM::TreeType::Tree::FindIdOfParent )
  .method("OrderNodes", &ppa::ParallelPruningAbcPOUMM::TreeType::Tree::OrderNodes )
  ;
  Rcpp::class_<ppa::ParallelPruningAbcPOUMM::TreeType>( "PruningTree" )
    .derives<ppa::ParallelPruningAbcPOUMM::TreeType::Tree> ( "Tree" )

  //.method("RangeIdPrune", &ppa::PruningTree<uint, double>::RangeIdPrune )
  //.method("RangeIdUpdateParent", &ppa::PruningTree<uint, double>::RangeIdUpdateParent )
    .method("CalculateHeights", &ppa::ParallelPruningAbcPOUMM::TreeType::CalculateHeights )
    .property("num_levels", &ppa::ParallelPruningAbcPOUMM::TreeType::num_levels )
    .property("ranges_id_visit", &ppa::ParallelPruningAbcPOUMM::TreeType::ranges_id_visit )
    .property("ranges_id_prune", &ppa::ParallelPruningAbcPOUMM::TreeType::ranges_id_prune )
  ;
  Rcpp::class_<ppa::ParallelPruningAbcPOUMM::PruningAlgorithmType> ( "ParallelPruning" )
  .property( "VersionOPENMP", &ppa::ParallelPruningAbcPOUMM::PruningAlgorithmType::VersionOPENMP )
  .method( "ModeAutoStep", &ppa::ParallelPruningAbcPOUMM::PruningAlgorithmType::ModeAutoStep )
  .property( "ModeAutoCurrent", &ppa::ParallelPruningAbcPOUMM::PruningAlgorithmType::ModeAutoCurrent )
  .property( "IsTuning", &ppa::ParallelPruningAbcPOUMM::PruningAlgorithmType::IsTuning )
  .property( "num_threads", &ppa::ParallelPruningAbcPOUMM::PruningAlgorithmType::num_threads )
  .property( "min_size_chunk_visit", &ppa::ParallelPruningAbcPOUMM::PruningAlgorithmType::min_size_chunk_visit )
  .property( "min_size_chunk_prune", &ppa::ParallelPruningAbcPOUMM::PruningAlgorithmType::min_size_chunk_prune )
  .property( "durations_tuning", &ppa::ParallelPruningAbcPOUMM::PruningAlgorithmType::durations_tuning )
  .property( "fastest_step_tuning", &ppa::ParallelPruningAbcPOUMM::PruningAlgorithmType::fastest_step_tuning )
  ;
  Rcpp::class_<ppa::ParallelPruningAbcPOUMM::PruningSpecType> ( "PruningSpec" )
  .field( "a", &ppa::ParallelPruningAbcPOUMM::PruningSpecType::a )
  .field( "b", &ppa::ParallelPruningAbcPOUMM::PruningSpecType::b )
  .field( "c", &ppa::ParallelPruningAbcPOUMM::PruningSpecType::c )
  .field( "se", &ppa::ParallelPruningAbcPOUMM::PruningSpecType::se )
  .field( "z", &ppa::ParallelPruningAbcPOUMM::PruningSpecType::z )
  .field( "alpha", &ppa::ParallelPruningAbcPOUMM::PruningSpecType::alpha )
  .field( "sigma2", &ppa::ParallelPruningAbcPOUMM::PruningSpecType::sigma2 )
  .field( "sigmae2", &ppa::ParallelPruningAbcPOUMM::PruningSpecType::sigmae2 )
  .field( "theta", &ppa::ParallelPruningAbcPOUMM::PruningSpecType::theta )
  ;
  Rcpp::class_<ppa::ParallelPruningAbcPOUMM>( "ParallelPruningAbcPOUMM" )
  .factory<Rcpp::List const&, ppa::vec const&, ppa::vec const&>(&CreateParallelPruningAbcPOUMM)
  .method( "DoPruning", &ppa::ParallelPruningAbcPOUMM::DoPruning )
  .property( "tree", &ppa::ParallelPruningAbcPOUMM::tree )
  .property( "spec", &ppa::ParallelPruningAbcPOUMM::spec )
  .property( "algorithm", &ppa::ParallelPruningAbcPOUMM::algorithm )
  ;
}
