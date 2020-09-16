// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// san_timediscrete_c
List san_timediscrete_c(const List C_init, const NumericVector p_S, const NumericVector p_0, const NumericVector p_R, const NumericVector p_A, const NumericVector p_N, const NumericVector p_D, const IntegerVector steps);
RcppExport SEXP _SANsimulatoR_san_timediscrete_c(SEXP C_initSEXP, SEXP p_SSEXP, SEXP p_0SEXP, SEXP p_RSEXP, SEXP p_ASEXP, SEXP p_NSEXP, SEXP p_DSEXP, SEXP stepsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const List >::type C_init(C_initSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type p_S(p_SSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type p_0(p_0SEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type p_R(p_RSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type p_A(p_ASEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type p_N(p_NSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type p_D(p_DSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type steps(stepsSEXP);
    rcpp_result_gen = Rcpp::wrap(san_timediscrete_c(C_init, p_S, p_0, p_R, p_A, p_N, p_D, steps));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_SANsimulatoR_san_timediscrete_c", (DL_FUNC) &_SANsimulatoR_san_timediscrete_c, 8},
    {NULL, NULL, 0}
};

RcppExport void R_init_SANsimulatoR(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
