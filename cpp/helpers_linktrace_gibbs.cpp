#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector sv_math(NumericVector x, double num, String operation){
  NumericVector ret(x.size());
  for(int i = 0; i < ret.size(); i++){

    if(operation == "multiply"){
      ret[i] = x[i] * num;
    }

    if(operation == "add"){
      ret[i] = x[i] + num;
    }

    if(operation == "subtract"){
      ret[i] = x[i] - num;
    }

    if(operation == "divide"){
      ret[i] = x[i] / num;
    }

  }
  return ret;
}

// [[Rcpp::export]]
IntegerVector sv_math_int(IntegerVector x, double num, String operation){
  IntegerVector ret(x.size());
  for(int i = 0; i < ret.size(); i++){

    if(operation == "multiply"){
      ret[i] = x[i] * num;
    }

    if(operation == "add"){
      ret[i] = x[i] + num;
    }

    if(operation == "subtract"){
      ret[i] = x[i] - num;
    }

    if(operation == "divide"){
      ret[i] = x[i] / num;
    }

  }
  return ret;
}

// [[Rcpp::export]]
IntegerVector int_vec_insert(IntegerVector vec, IntegerVector vals, IntegerVector pos){
  for(int i = 0; i < vals.size(); i++){
    int pos_i = pos[i];
    vec[pos_i] = vals[i];
  }
  return vec;
}

// [[Rcpp::export]]
IntegerMatrix matrix_extract(){}

// [[Rcpp::export]]
IntegerMatrix matrix_insert(){}


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
  Function c_rep("rep");
  Function c_expand_grid("expand.grid");

  List data_p_waves = lt_permute(data);
  DataFrame data_p = data;

  IntegerVector n_p(n_waves);

  for(int i = 0; i < n_waves; i++){
    IntegerVector npi = data_p_waves[i];
    n_p[i] = npi.size();
  }

  List data_p_reorder(n_waves);
  data_p_reorder[0] = Range(1,n_p[0]);

  for(int i = 1; i < n_waves; i++){
    data_p_reorder[i] = Range(sum(head(n_p,i)) + 1, sum(head(n_p, i + 1)));
  }

  NumericMatrix l(chain_samples,n_strata);
  List b(chain_samples);
  IntegerVector n(chain_samples);

  l(0,_) = as<NumericVector>(param_init["l_0"]);
  b[0] = as<NumericMatrix>(param_init["b_0"]);
  n[0] = as<int>(param_init["n_0"]);

  int prior_n = as<int>(priors["p_n"]);
  NumericVector prior_l = as<NumericVector>(priors["p_l"]);
  int prior_b = as<int>(priors["p_b"]);

  for(int t = 1; t < chain_samples; t++){

    // generate new N

    IntegerVector data_p_strata = data_p["strata"];
    IntegerVector rows = c_unlist(head(data_p_waves, n_waves - 1));
    IntegerVector rows_pull = clone(rows) - 1;

    IntegerVector strata_t(rows_pull.size());

    for(int i = 0; i < rows_pull.size(); i++){
      int pos = rows_pull[i];
      strata_t[i] = data_p_strata[pos];
    }

    IntegerVector strata_count = table(strata_t);

    NumericVector no_link_init = c_rep(1, n_strata);

    for(int i = 0; i < n_strata; i++){
      for(int j = 0; j < n_strata; j++){
        NumericMatrix bi = b[t - 1];
        no_link_init[i] = no_link_init[i] * pow((1 - bi(j,i)),strata_count[j]);
      }
    }

    NumericVector no_link_l = l(t - 1,_) * no_link_init;
    double no_link = sum(no_link_l);

    IntegerVector n_0 = data_p_waves[0];
    int nn_0 = n_0.size();
    IntegerVector nn_vec = c_unlist(data_p_waves);
    int nn = nn_vec.size();

    IntegerVector n_post_range = Range(nn, total * 5);

    NumericVector n_sample_prob_vec(n_post_range.size());

    for(int i = 0; i < n_post_range.size(); i++){
     IntegerVector r_i = Range(n_post_range[i] + 1 - nn, n_post_range[i] - nn_0);
     n_sample_prob_vec[i] =  sum(log(r_i)) +
       (n_post_range[i] - nn) * log(no_link) - prior_n * log(n_post_range[i]);
    }

    double n_sample_prob_max = max(n_sample_prob_vec);
    NumericVector n_sample_prob = exp(sv_math(n_sample_prob_vec, n_sample_prob_max, "subtract"));
    n[t] = sample(n_post_range, 1, false, n_sample_prob)[0];

    // assign strata to non sampled units

    IntegerVector data_p_reorder_unlist = c_unlist(data_p_reorder);
    IntegerVector n_range = Range(1, n[t]);
    IntegerVector not_sampled = setdiff(n_range, data_p_reorder_unlist);
    not_sampled = sv_math_int(not_sampled, 1, "subtract");

    IntegerVector stratum(n[t]);

    IntegerVector ins_pos_us = sv_math_int(data_p_reorder_unlist, 1, "subtract");
    IntegerVector data_p_waves_unlist = c_unlist(data_p_waves);
    rows_pull = clone(data_p_waves_unlist) - 1;
    IntegerVector ins_val_us(rows_pull.size());

    for(int i = 0; i < rows_pull.size(); i++){
      int pos = rows_pull[i];
      ins_val_us[i] = data_p_strata[pos];
    }

    stratum = int_vec_insert(stratum, ins_val_us, ins_pos_us);

    IntegerVector strat = Range(1, n_strata);
    NumericVector pstrat = no_link_l / no_link;
    IntegerVector stratsamp = sample(strat, not_sampled.size(), true, pstrat);

    stratum = int_vec_insert(stratum, stratsamp, not_sampled);
    DataFrame link_comb = c_expand_grid(Range(1,n_waves), Range(1,n_waves));
    IntegerVector g1 = link_comb[0];
    IntegerVector g2 = link_comb[1];

    for(int i = 0; i < g1.size - 1; i++){

    }

  }

  List ret(1);
  ret[0] = n;
  return ret;
}





