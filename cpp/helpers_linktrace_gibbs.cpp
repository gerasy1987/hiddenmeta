#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List lt_permute(DataFrame data){
  NumericVector wave = data["rds_wave"];
  LogicalVector n = wave == 1;
  int n_inital = sum(n);
  IntegerVector t = table(wave);
  int n_waves = t.length();

  List wave_samples(n_waves);
  NumericVector name = data["name"];
  bool replace = false;
  NumericVector s0 = sample(name, n_inital, replace);
  wave_samples[0] = s0;

  Function c_unlist("unlist");

  for(int i = 1; i < n_waves; i++){
    List l = data["links_list"];
    NumericVector prev_val =  wave_samples[i - 1];
    LogicalVector get_elem =  in(name, prev_val);
    List l_i = l[get_elem];
    NumericVector set1_temp = c_unlist(l_i);
    NumericVector set1 = unique(set1_temp);
    List set2_temp = wave_samples[Range(0,i - 1)];
    NumericVector set2 = c_unlist(set2_temp);
    NumericVector res_i = setdiff(set1,set2);
    wave_samples[i] = res_i;
  }
  return wave_samples;
}


// [[Rcpp::export]]
List lt_gibbs(DataFrame data, IntegerMatrix y_samp, IntegerVector strata, int n_strata,
              int n_waves, int total, int chain_samples, int chain_burnin, List priors,
              List param_init) {

  Function c_unlist("unlist");
  List data_p_waves = lt_permute(data);
  DataFrame data_p = clone(data);

  IntegerVector n_p(n_waves);

  for(int i = 0; i < n_waves; i++){
    n_p[i] = data_p_waves[i].size();
  }

  List data_p_reorder(n_waves);
  data_p_reorder[0] = Range(1,n_p[0]);

  for(int i = 1; i < n_waves; i++){
    data_p_reorder[i] = Range(sum(head(n_p,i)) + 1, sum(head(n_p, i + 1)));
  }

  NumericMatrix l(chain_samples,n_strata);
  List b(chain_samples);
  IntegerVector n(chain_samples);

  l[0,_] = param_init["l_0"];
  b[0] = param_init["b_0"];
  n[0] = param_init["n_0"];

  int prior_n = priors["p_n"];
  double prior_l = priors["p_l"];
  int prior_b = priors["p_b"];

  for(int t = 1; t < chain_samples; t++){

    IntegerVector data_p_strata = clone(data_p["strata"]);
    IntegerVector rows = c_unlist(head(data_p_waves, n_waves - 1));
    IntegerVector rows_pull = clone(rows) - 1;

    IntegerVector strata_t(rows_pull.size());

    for(int i = 0; i < rows_pull.size(); i++){
      int pos = rows_pull[i];
      strata_t[i] = data_p_strata[pos];
    }

    IntegerVector strata_count = table(strata_t);

    NumericVector no_link_init(1, n_strata);

    for(int i = 0; i < n_strata; i++){
      for(int j = 0; j < n_strata; j++){
        NumericMatrix bi = b[t - 1];
        no_link_init[i] = no_link_init[i] * (1 - bi(j,i))^strata_count[j];
      }
    }

    NumericVector no_link_l = l(t - 1,_) * no_link_init;
    double no_link = sum(no_link_l);

    int nn_0 =

  }
}





