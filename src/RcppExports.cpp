// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// LinCDECdf
NumericVector LinCDECdf(const NumericMatrix& cellProb, const NumericVector& y, const NumericVector& splitPointY);
RcppExport SEXP _LinCDE_LinCDECdf(SEXP cellProbSEXP, SEXP ySEXP, SEXP splitPointYSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type cellProb(cellProbSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type splitPointY(splitPointYSEXP);
    rcpp_result_gen = Rcpp::wrap(LinCDECdf(cellProb, y, splitPointY));
    return rcpp_result_gen;
END_RCPP
}
// LinCDEQuantiles
NumericMatrix LinCDEQuantiles(const NumericMatrix& cellProb, const NumericVector& probs, const NumericVector& splitPointY);
RcppExport SEXP _LinCDE_LinCDEQuantiles(SEXP cellProbSEXP, SEXP probsSEXP, SEXP splitPointYSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type cellProb(cellProbSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type probs(probsSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type splitPointY(splitPointYSEXP);
    rcpp_result_gen = Rcpp::wrap(LinCDEQuantiles(cellProb, probs, splitPointY));
    return rcpp_result_gen;
END_RCPP
}
// LinCDESplit
NumericVector LinCDESplit(const NumericMatrix& X, const IntegerVector& yIndex, const NumericMatrix& cellProb, const NumericMatrix& z, const NumericMatrix& covMatrixInv, const NumericMatrix& splitPoint, int numberBin);
RcppExport SEXP _LinCDE_LinCDESplit(SEXP XSEXP, SEXP yIndexSEXP, SEXP cellProbSEXP, SEXP zSEXP, SEXP covMatrixInvSEXP, SEXP splitPointSEXP, SEXP numberBinSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type yIndex(yIndexSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type cellProb(cellProbSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type z(zSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type covMatrixInv(covMatrixInvSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type splitPoint(splitPointSEXP);
    Rcpp::traits::input_parameter< int >::type numberBin(numberBinSEXP);
    rcpp_result_gen = Rcpp::wrap(LinCDESplit(X, yIndex, cellProb, z, covMatrixInv, splitPoint, numberBin));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_hello_world
List rcpp_hello_world();
RcppExport SEXP _LinCDE_rcpp_hello_world() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(rcpp_hello_world());
    return rcpp_result_gen;
END_RCPP
}