// [[Rcpp::export]]
IntegerVector int_vec_insert(IntegerVector vec,
                             IntegerVector vals,
                             IntegerVector pos){

  for(int i = 0; i < vals.size(); i++){
    int pos_i = pos[i];
    vec[pos_i] = vals[i];
  }
  return vec;
}


// [[Rcpp::export]]
IntegerMatrix int_mat_insert(IntegerMatrix old_m,
                             IntegerMatrix new_m,
                             IntegerVector new_rows,
                             IntegerVector new_cols,
                             IntegerVector old_rows,
                             IntegerVector new_rows){

  for(int i = 0; i < new_rows.lenght(); i++){
    for(int j = 0; j < new_cols.length(); j++){
      new_m(new_rows[i] - 1, new_cols[j] - 1) = old_m(old_rows[i] - 1, old_cols[j] - 1);
    }
  }
  return new_m;
}


// [[Rcpp::export]]
List lt_permute(DataFrame data){

  Function c_unlist("unlist");

  NumericVector wave = data["rds_wave"];
  int n_inital = sum(wave == 1);
  int n_waves = table(wave).length();

  List wave_samples(n_waves);
  NumericVector s0 = sample(as<NumericVector>(data["name"]), n_inital, false, R_NilValue);
  wave_samples[0] = s0;

  for(int i = 1; i < n_waves; i++){
    List l = data["links_list"];
    LogicalVector get_elem =  in(as<NumericVector>(data["name"]),
                                 as<NumericVector>(wave_samples[i - 1]));
    List l_i = l[get_elem];
    NumericVector set1 = unique(as<NumericVector>(c_unlist(l_i)));
    NumericVector set2 = c_unlist(wave_samples[Range(0,i - 1)]);
    wave_samples[i] = setdiff(set1,set2);
  }
  return wave_samples;
}



// [[Rcpp::export]]
List lt_gibbs(DataFrame data, IntegerMatrix y_samp, IntegerVector strata, int n_strata,
              int n_waves, int total, int chain_samples, int chain_burnin, List priors,
              List param_init) {

  Function c_unlist("unlist");
  Function c_expand_grid("expand.grid");

  // permute data
  List data_p_waves = lt_permute(data);
  DataFrame data_p = data;

  // reordering samples to estimate N
  IntegerVector n_p(n_waves);

  for(int i = 0; i < n_waves; i++){
    n_p[i] = as<IntegerVector>(data_p_waves[i]).size();
  }

  List data_p_reorder(n_waves);
  data_p_reorder[0] = Range(1,n_p[0]);

  for(int i = 1; i < n_waves; i++){
    data_p_reorder[i] = Range(sum(head(n_p,i)) + 1, sum(head(n_p, i + 1)));
  }

  //assign seeds
  NumericMatrix l(chain_samples,n_strata);
  List b(chain_samples);
  IntegerVector n(chain_samples);

  l(0,_) = as<NumericVector>(param_init["l_0"]);
  b[0] = as<NumericMatrix>(param_init["b_0"]);
  n[0] = as<int>(param_init["n_0"]);

  int prior_n = as<int>(priors["p_n"]);
  NumericVector prior_l = as<NumericVector>(priors["p_l"]);
  int prior_b = as<int>(priors["p_b"]);

  // begin MCMC
  for(int t = 1; t < chain_samples; t++){

    //##################
    //# generate new N #
    //##################

    //get number of units in each strata
    IntegerVector data_p_strata = data_p["strata"];
    IntegerVector rows = c_unlist(head(data_p_waves, n_waves - 1));
    IntegerVector rows_pull = clone(rows) - 1;

    IntegerVector strata_t(rows_pull.size());

    for(int i = 0; i < rows_pull.size(); i++){
      strata_t[i] = data_p_strata[rows_pull[i]];
    }

    IntegerVector strata_count = table(strata_t);

    //get p(no link between strata)
    NumericVector no_link_init = rep(1, n_strata);

    for(int i = 0; i < n_strata; i++){
      for(int j = 0; j < n_strata; j++){
        NumericMatrix bi = b[t - 1];
        no_link_init[i] = no_link_init[i] * pow((1 - bi(j,i)),strata_count[j]);
      }
    }

    NumericVector no_link_l = l(t - 1,_) * no_link_init;
    double no_link = sum(no_link_l);

    int nn_0 = as<IntegerVector>(data_p_waves[0]).size();
    int nn = as<IntegerVector>(c_unlist(data_p_waves)).size();

    IntegerVector n_post_range = Range(nn, total * 5);

    NumericVector n_sample_prob_vec(n_post_range.size());

    for(int i = 0; i < n_post_range.size(); i++){
     IntegerVector r_i = Range(n_post_range[i] + 1 - nn, n_post_range[i] - nn_0);
     n_sample_prob_vec[i] =  sum(log(r_i)) +
       (n_post_range[i] - nn) * log(no_link) - prior_n * log(n_post_range[i]);
    }

    NumericVector n_sample_prob = exp(n_sample_prob_vec - max(n_sample_prob_vec));
    n[t] = sample(n_post_range, 1, false, n_sample_prob)[0];


    //#######################
    //# generate new lambda #
    //#######################

    // assign strata to non sampled units

    // get indices of non sampled units
    IntegerVector n_range = Range(1, n[t]);
    IntegerVector not_sampled = setdiff(n_range, as<IntegerVector>(c_unlist(data_p_reorder))) - 1;

    IntegerVector stratum(n[t]);

    //fill stratum vector with strata of sampled units
    IntegerVector ins_pos_us = c_unlist(data_p_reorder) - 1;
    IntegerVector rows_pull = as<IntegerVector>(c_unlist(data_p_waves)) - 1;
    IntegerVector ins_val_us(rows_pull.size());

    for(int i = 0; i < rows_pull.size(); i++){
      ins_val_us[i] = data_p_strata[rows_pull[i]];
    }

    stratum = int_vec_insert(stratum, ins_val_us, ins_pos_us);

    //fill stratum vector with strata of non sampled units
    IntegerVector strat = Range(1, n_strata);
    NumericVector pstrat = no_link_l / no_link;
    IntegerVector stratsamp = sample(strat, not_sampled.size(), true, pstrat);
    stratum = int_vec_insert(stratum, stratsamp, not_sampled);


    // populate link matrix for reordered sample
    DataFrame link_comb = c_expand_grid(Range(0,n_waves - 1), Range(0,n_waves - 1));
    IntegerVector g1 = clone(link_comb[0]);
    IntegerVector g2 = clone(link_comb[1]);

    // for matrix replacement use matrix indeces of values (i.e. col 0 & row 2 == 2)
    // (col * nrow) + row (-1)?
    for(int i = 0; i < g1.size(); i++){

    }

  }

  List ret(1);
  ret[0] = n;
  return ret;
}





