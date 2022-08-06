// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// NeedForSpeed1
List NeedForSpeed1(SEXP D, SEXP S, SEXP pD, SEXP pS, SEXP nrow, SEXP N, SEXP N_len, SEXP ssq, SEXP B, SEXP overlaps, SEXP overlaps_P);
RcppExport SEXP _EVALFQ_NeedForSpeed1(SEXP DSEXP, SEXP SSEXP, SEXP pDSEXP, SEXP pSSEXP, SEXP nrowSEXP, SEXP NSEXP, SEXP N_lenSEXP, SEXP ssqSEXP, SEXP BSEXP, SEXP overlapsSEXP, SEXP overlaps_PSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type D(DSEXP);
    Rcpp::traits::input_parameter< SEXP >::type S(SSEXP);
    Rcpp::traits::input_parameter< SEXP >::type pD(pDSEXP);
    Rcpp::traits::input_parameter< SEXP >::type pS(pSSEXP);
    Rcpp::traits::input_parameter< SEXP >::type nrow(nrowSEXP);
    Rcpp::traits::input_parameter< SEXP >::type N(NSEXP);
    Rcpp::traits::input_parameter< SEXP >::type N_len(N_lenSEXP);
    Rcpp::traits::input_parameter< SEXP >::type ssq(ssqSEXP);
    Rcpp::traits::input_parameter< SEXP >::type B(BSEXP);
    Rcpp::traits::input_parameter< SEXP >::type overlaps(overlapsSEXP);
    Rcpp::traits::input_parameter< SEXP >::type overlaps_P(overlaps_PSEXP);
    rcpp_result_gen = Rcpp::wrap(NeedForSpeed1(D, S, pD, pS, nrow, N, N_len, ssq, B, overlaps, overlaps_P));
    return rcpp_result_gen;
END_RCPP
}
// NeedForSpeed2
List NeedForSpeed2(SEXP D, SEXP pD, SEXP nrow, SEXP N, SEXP N_len, SEXP B, SEXP overlaps, SEXP overlaps_P);
RcppExport SEXP _EVALFQ_NeedForSpeed2(SEXP DSEXP, SEXP pDSEXP, SEXP nrowSEXP, SEXP NSEXP, SEXP N_lenSEXP, SEXP BSEXP, SEXP overlapsSEXP, SEXP overlaps_PSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type D(DSEXP);
    Rcpp::traits::input_parameter< SEXP >::type pD(pDSEXP);
    Rcpp::traits::input_parameter< SEXP >::type nrow(nrowSEXP);
    Rcpp::traits::input_parameter< SEXP >::type N(NSEXP);
    Rcpp::traits::input_parameter< SEXP >::type N_len(N_lenSEXP);
    Rcpp::traits::input_parameter< SEXP >::type B(BSEXP);
    Rcpp::traits::input_parameter< SEXP >::type overlaps(overlapsSEXP);
    Rcpp::traits::input_parameter< SEXP >::type overlaps_P(overlaps_PSEXP);
    rcpp_result_gen = Rcpp::wrap(NeedForSpeed2(D, pD, nrow, N, N_len, B, overlaps, overlaps_P));
    return rcpp_result_gen;
END_RCPP
}
// pvalue
NumericVector pvalue(SEXP a, SEXP b);
RcppExport SEXP _EVALFQ_pvalue(SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type a(aSEXP);
    Rcpp::traits::input_parameter< SEXP >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(pvalue(a, b));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_EVALFQ_NeedForSpeed1", (DL_FUNC) &_EVALFQ_NeedForSpeed1, 11},
    {"_EVALFQ_NeedForSpeed2", (DL_FUNC) &_EVALFQ_NeedForSpeed2, 8},
    {"_EVALFQ_pvalue", (DL_FUNC) &_EVALFQ_pvalue, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_EVALFQ(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}