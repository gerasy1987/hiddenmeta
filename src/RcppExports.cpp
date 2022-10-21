// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// move_elements
std::vector<int> move_elements(std::vector<int> x, int old_index, int new_index);
RcppExport SEXP _hiddenmeta_move_elements(SEXP xSEXP, SEXP old_indexSEXP, SEXP new_indexSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<int> >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type old_index(old_indexSEXP);
    Rcpp::traits::input_parameter< int >::type new_index(new_indexSEXP);
    rcpp_result_gen = Rcpp::wrap(move_elements(x, old_index, new_index));
    return rcpp_result_gen;
END_RCPP
}
// choose_cpp
uint64_t choose_cpp(uint64_t n, uint64_t k);
RcppExport SEXP _hiddenmeta_choose_cpp(SEXP nSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< uint64_t >::type n(nSEXP);
    Rcpp::traits::input_parameter< uint64_t >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(choose_cpp(n, k));
    return rcpp_result_gen;
END_RCPP
}
// combn_cpp
arma::mat combn_cpp(std::vector<int> x, int K);
RcppExport SEXP _hiddenmeta_combn_cpp(SEXP xSEXP, SEXP KSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<int> >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    rcpp_result_gen = Rcpp::wrap(combn_cpp(x, K));
    return rcpp_result_gen;
END_RCPP
}
// rdirichlet_cpp
std::vector<double> rdirichlet_cpp(std::vector<double> alpha);
RcppExport SEXP _hiddenmeta_rdirichlet_cpp(SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<double> >::type alpha(alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(rdirichlet_cpp(alpha));
    return rcpp_result_gen;
END_RCPP
}
// table_cpp
std::vector<int> table_cpp(std::vector<int> x);
RcppExport SEXP _hiddenmeta_table_cpp(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<int> >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(table_cpp(x));
    return rcpp_result_gen;
END_RCPP
}
// rep_times
std::vector<int> rep_times(std::vector<int> x, int n);
RcppExport SEXP _hiddenmeta_rep_times(SEXP xSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<int> >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(rep_times(x, n));
    return rcpp_result_gen;
END_RCPP
}
// rep_each
std::vector<int> rep_each(std::vector<int> x, int n);
RcppExport SEXP _hiddenmeta_rep_each(SEXP xSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<int> >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(rep_each(x, n));
    return rcpp_result_gen;
END_RCPP
}
// gen_range
std::vector<int> gen_range(int from, int to);
RcppExport SEXP _hiddenmeta_gen_range(SEXP fromSEXP, SEXP toSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type from(fromSEXP);
    Rcpp::traits::input_parameter< int >::type to(toSEXP);
    rcpp_result_gen = Rcpp::wrap(gen_range(from, to));
    return rcpp_result_gen;
END_RCPP
}
// int_vec_insert
std::vector<int> int_vec_insert(std::vector<int> vec, std::vector<int> vals, std::vector<int> pos);
RcppExport SEXP _hiddenmeta_int_vec_insert(SEXP vecSEXP, SEXP valsSEXP, SEXP posSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<int> >::type vec(vecSEXP);
    Rcpp::traits::input_parameter< std::vector<int> >::type vals(valsSEXP);
    Rcpp::traits::input_parameter< std::vector<int> >::type pos(posSEXP);
    rcpp_result_gen = Rcpp::wrap(int_vec_insert(vec, vals, pos));
    return rcpp_result_gen;
END_RCPP
}
// mat_to_mat_insert
arma::mat mat_to_mat_insert(arma::mat old_m, arma::mat new_m, std::vector<int> new_rows, std::vector<int> new_cols, std::vector<int> old_rows, std::vector<int> old_cols);
RcppExport SEXP _hiddenmeta_mat_to_mat_insert(SEXP old_mSEXP, SEXP new_mSEXP, SEXP new_rowsSEXP, SEXP new_colsSEXP, SEXP old_rowsSEXP, SEXP old_colsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type old_m(old_mSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type new_m(new_mSEXP);
    Rcpp::traits::input_parameter< std::vector<int> >::type new_rows(new_rowsSEXP);
    Rcpp::traits::input_parameter< std::vector<int> >::type new_cols(new_colsSEXP);
    Rcpp::traits::input_parameter< std::vector<int> >::type old_rows(old_rowsSEXP);
    Rcpp::traits::input_parameter< std::vector<int> >::type old_cols(old_colsSEXP);
    rcpp_result_gen = Rcpp::wrap(mat_to_mat_insert(old_m, new_m, new_rows, new_cols, old_rows, old_cols));
    return rcpp_result_gen;
END_RCPP
}
// lt_permute
std::vector<std::vector<int>> lt_permute(List link_list, std::vector<int> wave, std::vector<int> name);
RcppExport SEXP _hiddenmeta_lt_permute(SEXP link_listSEXP, SEXP waveSEXP, SEXP nameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type link_list(link_listSEXP);
    Rcpp::traits::input_parameter< std::vector<int> >::type wave(waveSEXP);
    Rcpp::traits::input_parameter< std::vector<int> >::type name(nameSEXP);
    rcpp_result_gen = Rcpp::wrap(lt_permute(link_list, wave, name));
    return rcpp_result_gen;
END_RCPP
}
// lt_gibbs
std::vector<int> lt_gibbs(DataFrame data, arma::mat y_samp, std::vector<int> strata, int n_strata, int n_waves, int total, int chain_samples, List priors, List param_init);
RcppExport SEXP _hiddenmeta_lt_gibbs(SEXP dataSEXP, SEXP y_sampSEXP, SEXP strataSEXP, SEXP n_strataSEXP, SEXP n_wavesSEXP, SEXP totalSEXP, SEXP chain_samplesSEXP, SEXP priorsSEXP, SEXP param_initSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame >::type data(dataSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type y_samp(y_sampSEXP);
    Rcpp::traits::input_parameter< std::vector<int> >::type strata(strataSEXP);
    Rcpp::traits::input_parameter< int >::type n_strata(n_strataSEXP);
    Rcpp::traits::input_parameter< int >::type n_waves(n_wavesSEXP);
    Rcpp::traits::input_parameter< int >::type total(totalSEXP);
    Rcpp::traits::input_parameter< int >::type chain_samples(chain_samplesSEXP);
    Rcpp::traits::input_parameter< List >::type priors(priorsSEXP);
    Rcpp::traits::input_parameter< List >::type param_init(param_initSEXP);
    rcpp_result_gen = Rcpp::wrap(lt_gibbs(data, y_samp, strata, n_strata, n_waves, total, chain_samples, priors, param_init));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_hiddenmeta_move_elements", (DL_FUNC) &_hiddenmeta_move_elements, 3},
    {"_hiddenmeta_choose_cpp", (DL_FUNC) &_hiddenmeta_choose_cpp, 2},
    {"_hiddenmeta_combn_cpp", (DL_FUNC) &_hiddenmeta_combn_cpp, 2},
    {"_hiddenmeta_rdirichlet_cpp", (DL_FUNC) &_hiddenmeta_rdirichlet_cpp, 1},
    {"_hiddenmeta_table_cpp", (DL_FUNC) &_hiddenmeta_table_cpp, 1},
    {"_hiddenmeta_rep_times", (DL_FUNC) &_hiddenmeta_rep_times, 2},
    {"_hiddenmeta_rep_each", (DL_FUNC) &_hiddenmeta_rep_each, 2},
    {"_hiddenmeta_gen_range", (DL_FUNC) &_hiddenmeta_gen_range, 2},
    {"_hiddenmeta_int_vec_insert", (DL_FUNC) &_hiddenmeta_int_vec_insert, 3},
    {"_hiddenmeta_mat_to_mat_insert", (DL_FUNC) &_hiddenmeta_mat_to_mat_insert, 6},
    {"_hiddenmeta_lt_permute", (DL_FUNC) &_hiddenmeta_lt_permute, 3},
    {"_hiddenmeta_lt_gibbs", (DL_FUNC) &_hiddenmeta_lt_gibbs, 9},
    {NULL, NULL, 0}
};

RcppExport void R_init_hiddenmeta(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}