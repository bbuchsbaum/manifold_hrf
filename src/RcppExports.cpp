#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;

// declarations from distance_knn.cpp
arma::mat pairwise_distances_cpp(const arma::mat& M);
Rcpp::List knn_search_cpp(const arma::mat& data, const arma::mat& query, int k);

extern "C" {

SEXP _manifoldhrf_pairwise_distances_cpp(SEXP MSEXP) {
    BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< arma::mat >::type M(MSEXP);
    rcpp_result_gen = Rcpp::wrap(pairwise_distances_cpp(M));
    return rcpp_result_gen;
    END_RCPP
}

SEXP _manifoldhrf_knn_search_cpp(SEXP dataSEXP, SEXP querySEXP, SEXP kSEXP) {
    BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< arma::mat >::type data(dataSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type query(querySEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(knn_search_cpp(data, query, k));
    return rcpp_result_gen;
    END_RCPP
}

}

static const R_CallMethodDef CallEntries[] = {
    {"_manifoldhrf_pairwise_distances_cpp", (DL_FUNC) &_manifoldhrf_pairwise_distances_cpp, 1},
    {"_manifoldhrf_knn_search_cpp", (DL_FUNC) &_manifoldhrf_knn_search_cpp, 3},
    {NULL, NULL, 0}
};

extern "C" void R_init_manifoldhrf(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
