// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// getSij
NumericMatrix getSij(const NumericMatrix& expmu, const NumericVector& expdelta, const IntegerVector& cdindex);
RcppExport SEXP _BLPestimatoR_getSij(SEXP expmuSEXP, SEXP expdeltaSEXP, SEXP cdindexSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type expmu(expmuSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type expdelta(expdeltaSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type cdindex(cdindexSEXP);
    rcpp_result_gen = Rcpp::wrap(getSij(expmu, expdelta, cdindex));
    return rcpp_result_gen;
END_RCPP
}
// getSjtMod
NumericVector getSjtMod(const NumericMatrix& expmu, const NumericVector& expdelta, const int& nprodt, const int& startpos, const NumericVector& weights);
RcppExport SEXP _BLPestimatoR_getSjtMod(SEXP expmuSEXP, SEXP expdeltaSEXP, SEXP nprodtSEXP, SEXP startposSEXP, SEXP weightsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type expmu(expmuSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type expdelta(expdeltaSEXP);
    Rcpp::traits::input_parameter< const int& >::type nprodt(nprodtSEXP);
    Rcpp::traits::input_parameter< const int& >::type startpos(startposSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type weights(weightsSEXP);
    rcpp_result_gen = Rcpp::wrap(getSjtMod(expmu, expdelta, nprodt, startpos, weights));
    return rcpp_result_gen;
END_RCPP
}
// getExpMu
NumericMatrix getExpMu(const NumericMatrix& theta2Matrix, const NumericMatrix& qv, const NumericMatrix& Xrandom, const IntegerVector& cdid, const NumericMatrix& demographics);
RcppExport SEXP _BLPestimatoR_getExpMu(SEXP theta2MatrixSEXP, SEXP qvSEXP, SEXP XrandomSEXP, SEXP cdidSEXP, SEXP demographicsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type theta2Matrix(theta2MatrixSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type qv(qvSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type Xrandom(XrandomSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type cdid(cdidSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type demographics(demographicsSEXP);
    rcpp_result_gen = Rcpp::wrap(getExpMu(theta2Matrix, qv, Xrandom, cdid, demographics));
    return rcpp_result_gen;
END_RCPP
}
// getDelta
List getDelta(const NumericMatrix& theta2, const NumericVector& deltaOld, const IntegerVector& cdid, const IntegerVector& cdindex, const NumericMatrix& Xrandom, const NumericVector& obsshare, const double& innerCrit, const int& innerMaxit, const int& printLevel, const NumericMatrix& indices, const NumericMatrix& nodesRcMktShape, const NumericMatrix& nodesDemMktShape, const NumericVector& weights);
RcppExport SEXP _BLPestimatoR_getDelta(SEXP theta2SEXP, SEXP deltaOldSEXP, SEXP cdidSEXP, SEXP cdindexSEXP, SEXP XrandomSEXP, SEXP obsshareSEXP, SEXP innerCritSEXP, SEXP innerMaxitSEXP, SEXP printLevelSEXP, SEXP indicesSEXP, SEXP nodesRcMktShapeSEXP, SEXP nodesDemMktShapeSEXP, SEXP weightsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type theta2(theta2SEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type deltaOld(deltaOldSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type cdid(cdidSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type cdindex(cdindexSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type Xrandom(XrandomSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type obsshare(obsshareSEXP);
    Rcpp::traits::input_parameter< const double& >::type innerCrit(innerCritSEXP);
    Rcpp::traits::input_parameter< const int& >::type innerMaxit(innerMaxitSEXP);
    Rcpp::traits::input_parameter< const int& >::type printLevel(printLevelSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type indices(indicesSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type nodesRcMktShape(nodesRcMktShapeSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type nodesDemMktShape(nodesDemMktShapeSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type weights(weightsSEXP);
    rcpp_result_gen = Rcpp::wrap(getDelta(theta2, deltaOld, cdid, cdindex, Xrandom, obsshare, innerCrit, innerMaxit, printLevel, indices, nodesRcMktShape, nodesDemMktShape, weights));
    return rcpp_result_gen;
END_RCPP
}
// dstddelta_c
arma::mat dstddelta_c(arma::mat& sijt, arma::mat& weights);
RcppExport SEXP _BLPestimatoR_dstddelta_c(SEXP sijtSEXP, SEXP weightsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type sijt(sijtSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type weights(weightsSEXP);
    rcpp_result_gen = Rcpp::wrap(dstddelta_c(sijt, weights));
    return rcpp_result_gen;
END_RCPP
}
// dstdtheta_c
arma::mat dstdtheta_c(arma::mat& sijt_arma, const NumericMatrix& indices, arma::mat& xt_arma, arma::mat& qvt_arma, arma::mat& dt_arma, arma::mat& weights_arma);
RcppExport SEXP _BLPestimatoR_dstdtheta_c(SEXP sijt_armaSEXP, SEXP indicesSEXP, SEXP xt_armaSEXP, SEXP qvt_armaSEXP, SEXP dt_armaSEXP, SEXP weights_armaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type sijt_arma(sijt_armaSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type indices(indicesSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type xt_arma(xt_armaSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type qvt_arma(qvt_armaSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type dt_arma(dt_armaSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type weights_arma(weights_armaSEXP);
    rcpp_result_gen = Rcpp::wrap(dstdtheta_c(sijt_arma, indices, xt_arma, qvt_arma, dt_arma, weights_arma));
    return rcpp_result_gen;
END_RCPP
}
// jacob_c
NumericMatrix jacob_c(NumericMatrix& sij, const NumericMatrix& indices, const List& blp_data, const List& blp_parameters, const List& blp_integration, const int& printLevel);
RcppExport SEXP _BLPestimatoR_jacob_c(SEXP sijSEXP, SEXP indicesSEXP, SEXP blp_dataSEXP, SEXP blp_parametersSEXP, SEXP blp_integrationSEXP, SEXP printLevelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix& >::type sij(sijSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type indices(indicesSEXP);
    Rcpp::traits::input_parameter< const List& >::type blp_data(blp_dataSEXP);
    Rcpp::traits::input_parameter< const List& >::type blp_parameters(blp_parametersSEXP);
    Rcpp::traits::input_parameter< const List& >::type blp_integration(blp_integrationSEXP);
    Rcpp::traits::input_parameter< const int& >::type printLevel(printLevelSEXP);
    rcpp_result_gen = Rcpp::wrap(jacob_c(sij, indices, blp_data, blp_parameters, blp_integration, printLevel));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_BLPestimatoR_getSij", (DL_FUNC) &_BLPestimatoR_getSij, 3},
    {"_BLPestimatoR_getSjtMod", (DL_FUNC) &_BLPestimatoR_getSjtMod, 5},
    {"_BLPestimatoR_getExpMu", (DL_FUNC) &_BLPestimatoR_getExpMu, 5},
    {"_BLPestimatoR_getDelta", (DL_FUNC) &_BLPestimatoR_getDelta, 13},
    {"_BLPestimatoR_dstddelta_c", (DL_FUNC) &_BLPestimatoR_dstddelta_c, 2},
    {"_BLPestimatoR_dstdtheta_c", (DL_FUNC) &_BLPestimatoR_dstdtheta_c, 6},
    {"_BLPestimatoR_jacob_c", (DL_FUNC) &_BLPestimatoR_jacob_c, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_BLPestimatoR(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
