// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393
#include <RcppArmadillo.h>
#include <Rcpp.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;


// matprod
mat matprod(mat amat, mat bmat);
RcppExport SEXP _quasipoissonprobs_matprod(SEXP amatSEXP, SEXP bmatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< mat >::type amat(amatSEXP);
    Rcpp::traits::input_parameter< mat >::type bmat(bmatSEXP);
    rcpp_result_gen = Rcpp::wrap(matprod(amat, bmat));
    return rcpp_result_gen;
END_RCPP
}
// mattridotprod
mat mattridotprod(mat amat, mat bmat, mat cmat);
RcppExport SEXP _quasipoissonprobs_mattridotprod(SEXP amatSEXP, SEXP bmatSEXP, SEXP cmatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< mat >::type amat(amatSEXP);
    Rcpp::traits::input_parameter< mat >::type bmat(bmatSEXP);
    Rcpp::traits::input_parameter< mat >::type cmat(cmatSEXP);
    rcpp_result_gen = Rcpp::wrap(mattridotprod(amat, bmat, cmat));
    return rcpp_result_gen;
END_RCPP
}
// mattridotprod2
mat mattridotprod2(vec amat, mat bmat, mat cmat);
RcppExport SEXP _quasipoissonprobs_mattridotprod2(SEXP amatSEXP, SEXP bmatSEXP, SEXP cmatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< vec >::type amat(amatSEXP);
    Rcpp::traits::input_parameter< mat >::type bmat(bmatSEXP);
    Rcpp::traits::input_parameter< mat >::type cmat(cmatSEXP);
    rcpp_result_gen = Rcpp::wrap(mattridotprod2(amat, bmat, cmat));
    return rcpp_result_gen;
END_RCPP
}
// lambdakfactorial_alt3
mat lambdakfactorial_alt3(vec L, vec eL, int maxy);
RcppExport SEXP _quasipoissonprobs_lambdakfactorial_alt3(SEXP LSEXP, SEXP eLSEXP, SEXP maxySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< vec >::type L(LSEXP);
    Rcpp::traits::input_parameter< vec >::type eL(eLSEXP);
    Rcpp::traits::input_parameter< int >::type maxy(maxySEXP);
    rcpp_result_gen = Rcpp::wrap(lambdakfactorial_alt3(L, eL, maxy));
    return rcpp_result_gen;
END_RCPP
}
// obs_min_expected
mat obs_min_expected(vec L, int maxy);
RcppExport SEXP _quasipoissonprobs_obs_min_expected(SEXP LSEXP, SEXP maxySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< vec >::type L(LSEXP);
    Rcpp::traits::input_parameter< int >::type maxy(maxySEXP);
    rcpp_result_gen = Rcpp::wrap(obs_min_expected(L, maxy));
    return rcpp_result_gen;
END_RCPP
}


static const R_CallMethodDef CallEntries[] = {
    {"_quasipoissonprobs_matprod", (DL_FUNC) &_quasipoissonprobs_matprod, 2},
    {"_quasipoissonprobs_mattridotprod", (DL_FUNC) &_quasipoissonprobs_mattridotprod, 3},
    {"_quasipoissonprobs_mattridotprod2", (DL_FUNC) &_quasipoissonprobs_mattridotprod2, 3},
    {"_quasipoissonprobs_lambdakfactorial_alt3", (DL_FUNC) &_quasipoissonprobs_lambdakfactorial_alt3, 3},
    {"_quasipoissonprobs_obs_min_expected", (DL_FUNC) &_quasipoissonprobs_obs_min_expected, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_quasipoissonprobs(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
