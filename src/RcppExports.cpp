// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// log_n_stop_C
void log_n_stop_C(std::string message, std::string m_details);
RcppExport SEXP _posir_log_n_stop_C(SEXP messageSEXP, SEXP m_detailsSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type message(messageSEXP);
    Rcpp::traits::input_parameter< std::string >::type m_details(m_detailsSEXP);
    log_n_stop_C(message, m_details);
    return R_NilValue;
END_RCPP
}
// n_traj_simu_1D_C
arma::mat n_traj_simu_1D_C(arma::mat const& Z, arma::ivec const& Intdgrid, bool is_standard);
RcppExport SEXP _posir_n_traj_simu_1D_C(SEXP ZSEXP, SEXP IntdgridSEXP, SEXP is_standardSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat const& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< arma::ivec const& >::type Intdgrid(IntdgridSEXP);
    Rcpp::traits::input_parameter< bool >::type is_standard(is_standardSEXP);
    rcpp_result_gen = Rcpp::wrap(n_traj_simu_1D_C(Z, Intdgrid, is_standard));
    return rcpp_result_gen;
END_RCPP
}
// n_traj_simu_2D_C
arma::vec n_traj_simu_2D_C(arma::mat const& Z, arma::ivec const& Intdgrid, bool is_standard);
RcppExport SEXP _posir_n_traj_simu_2D_C(SEXP ZSEXP, SEXP IntdgridSEXP, SEXP is_standardSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat const& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< arma::ivec const& >::type Intdgrid(IntdgridSEXP);
    Rcpp::traits::input_parameter< bool >::type is_standard(is_standardSEXP);
    rcpp_result_gen = Rcpp::wrap(n_traj_simu_2D_C(Z, Intdgrid, is_standard));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_posir_log_n_stop_C", (DL_FUNC) &_posir_log_n_stop_C, 2},
    {"_posir_n_traj_simu_1D_C", (DL_FUNC) &_posir_n_traj_simu_1D_C, 3},
    {"_posir_n_traj_simu_2D_C", (DL_FUNC) &_posir_n_traj_simu_2D_C, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_posir(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
