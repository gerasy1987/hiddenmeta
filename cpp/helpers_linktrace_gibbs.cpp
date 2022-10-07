#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// n choose k helper
// [[Rcpp::export]]
uint64_t choose_cpp(uint64_t n, uint64_t k) {
  if(k == 0) return 1;
  return (n * choose_cpp(n - 1, k - 1)) / k;
}

// helper to generate combinations of x k at a time (cpp implementation of combn)
// [[Rcpp::export]]
arma::mat combn_cpp(std::vector<int> x, int K) {
  int N = x.size();

  std::string bitmask(K, 1);
  bitmask.resize(N, 0);

  uint64_t n_combos = choose_cpp(N,K);
  arma::mat results = arma::mat(n_combos, K);
  uint64_t row_position = 0;
  do {
    uint64_t col_position = 0;
    for (int i = 0; i < N; ++i)  {
      if (bitmask[i]) {
        results(row_position, col_position) = x[i];
        col_position++;
      }
    }
    row_position++;
  } while (std::prev_permutation(bitmask.begin(), bitmask.end()));
  return results;
}


// cpp implementation of rep
// [[Rcpp::export]]
std::vector<int> rep_times(std::vector<int> x, int n){

  std::vector<int> ret;

  for(int i = 0; i < n; i++){
    ret.insert(ret.end(),x.begin(),x.end());
  }

  return ret;
}

// cpp implementation of rep with each argument
// [[Rcpp::export]]
std::vector<int> rep_each(std::vector<int> x, int n){

  std::vector<int> ret;

  for(int i = 0; i < x.size(); i++){
    for(int j = 0; j < n; j++){
      ret.push_back(x[i]);
    }
  }

  return ret;
}

// helper to generate range of integers
// [[Rcpp::export]]
std::vector<int> gen_range(int from,
                           int to){
  std::vector<int> ret;
  ret.push_back(from);

  for(int i = 0; i < to - from; i++){
    ret.push_back(ret[i] + 1);
  }

  return ret;
}

// helper to do random access insertion for vectors
// [[Rcpp::export]]
std::vector<int> int_vec_insert(std::vector<int> vec,
                                std::vector<int> vals,
                                std::vector<int> pos){

  for(int i = 0; i < vals.size(); i++){
    vec[pos[i]] = vals[i];
  }

  return vec;
}

// helper to do random access insertions of vector elements into matrix
// [[Rcpp::export]]
IntegerMatrix int_mat_insert(IntegerMatrix m,
                             IntegerVector col,
                             IntegerVector row,
                             IntegerVector vals){

  int nrow = m.nrow();
  IntegerVector pos = ((col - 1) * nrow) + row - 1;

  for(int i = 0; i < pos.length(); i++){
    m[pos[i]] = vals[i];
  }

  return m;
}


// helper to do random access insertion of matrix into another matrix
// [[Rcpp::export]]
arma::mat mat_to_mat_insert(arma::mat old_m,
                            arma::mat new_m,
                            std::vector<int> new_rows,
                            std::vector<int> new_cols,
                            std::vector<int> old_rows,
                            std::vector<int> old_cols){

  for(int i = 0; i < new_rows.size(); i++){
    for(int j = 0; j < new_cols.size(); j++){
      new_m(new_rows[i] - 1, new_cols[j] - 1) = old_m(old_rows[i] - 1, old_cols[j] - 1);
    }
  }
  return new_m;
}

// helper to flatten matrix
// [[Rcpp::export]]
NumericVector mat_by_mat(NumericMatrix m,
                         IntegerVector row,
                         IntegerVector col){
  NumericVector res(row.length());

  for(int i = 0; i < res.length(); i++){
    res[i] = m(row[i] - 1, col[i] - 1);
  }

  return res;
}


// helper to turn R matrix into double vector of vectors
// [[Rcpp::export]]
std::vector<std::vector<double>> m_to_v_double(NumericMatrix m){

  std::vector<std::vector<double>> v;

  for(int i = 0; i < m.nrow(); i++){

    std::vector<double> row_i;

    for(int j = 0; j < m.ncol(); j++){
      row_i.push_back(m(i,j));
    }

    v.push_back(row_i);
  }

  return v;

}


// helper to permute sampling data
// [[Rcpp::export]]
std::vector<std::vector<int>> lt_permute(List link_list,
                                         std::vector<int> wave,
                                         std::vector<int> name){

  int n_inital = std::count(wave.begin(), wave.end(),1);
  sort(wave.begin(), wave.end());
  wave.erase(unique(wave.begin(), wave.end()), wave.end());
  int n_waves = wave.size();

  std::vector<std::vector<int>> wave_samples;
  std::vector<int> s_0(name);
  std::random_shuffle(s_0.begin(), s_0.end());
  s_0.resize(n_inital);
  wave_samples.push_back(s_0);

  for(int i = 1; i < n_waves; i++){
    std::vector<int> set1;
    std::vector<int> wave_samples_i = wave_samples[i-1];

    for(int j = 0; j < wave_samples_i.size(); j++){
      std::vector<int> l_j = link_list[wave_samples_i[j] - 1];
      set1.insert(set1.end(), l_j.begin(), l_j.end());
    }

    set1.erase(unique(set1.begin(), set1.end()),set1.end());
    sort(set1.begin(), set1.end());

    std::vector<int> set2;
    for(int k = 0; k < i; k++){
      std::vector<int> wave_sample_k = wave_samples[k];
      set2.insert(set2.end(), wave_sample_k.begin(), wave_sample_k.end());
    }
    sort(set2.begin(), set2.end());

    std::vector<int> ret;
    std::set_difference(set1.begin(),set1.end(),set2.begin(), set2.end(),
                        std::inserter(ret, ret.end()));

    ret.erase(unique(ret.begin(), ret.end()), ret.end());
    wave_samples.push_back(ret);
  }

  return wave_samples;
}



// gibbs sampler
// [[Rcpp::export]]
List lt_gibbs(DataFrame data,
              arma::mat y_samp,
              std::vector<int> strata,
              int n_strata,
              int n_waves,
              int total,
              int chain_samples,
              int chain_burnin,
              List priors,
              List param_init) {

  //Function c_unlist("unlist");
  //Function c_expand_grid("expand.grid");
  //Function c_combn("combn");
  //Function c_setdiff("setdiff");
  //Function c_which("which");
  //Function c_rdirichlet("rdirichlet");
  Function c_sample("sample");

  // permute data
  std::vector<std::vector<int>> data_p_waves = lt_permute(data["links_list"],data["rds_wave"],data["name"]);
  std::vector<int> data_p_waves_id;

  for(int i = 0; i < data_p_waves.size(); i++){
    for(int j = 0; j < data_p_waves[i].size(); j++){
      data_p_waves_id.push_back(data_p_waves[i][j] - 1);
    }
  }

  // reordering samples to estimate N
  std::vector<int> n_p;

  for(int i = 0; i < n_waves; i++){
    n_p.push_back(data_p_waves[i].size());
  }

  std::vector<std::vector<int>>  data_p_reorder;
  data_p_reorder.push_back(gen_range(1,n_p[0]));

  for(int i = 1; i < n_waves; i++){

    int from_i = std::accumulate(n_p.begin(), n_p.begin() + i, 0) + 1;
    int to_i = std::accumulate(n_p.begin(), n_p.begin() + i + 1, 0);
    data_p_reorder.push_back(gen_range(from_i,to_i));
  }

  std::vector<int> data_p_reorder_id;

  for(int i = 0; i < data_p_reorder.size(); i++){
    for(int j = 0; j < data_p_reorder[i].size(); j++){
      data_p_reorder_id.push_back(data_p_reorder[i][j] - 1);
    }
  }


  //assign seeds
  std::vector<std::vector<double>> l;
  std::vector<std::vector<std::vector<double>>> b;
  std::vector<int> n;

  l.push_back(param_init["l_0"]);
  b.push_back(m_to_v_double(as<NumericMatrix>(param_init["b_0"])));
  n.push_back(param_init["n_0"]);

  int prior_n = priors["p_n"];
  std::vector<double> prior_l = priors["p_l"];
  int prior_b = priors["p_b"];

  // begin MCMC
  for(int t = 1; t < chain_samples; t++){

    //##################
    //# generate new N #
    //##################

    //get number of units in each strata
    std::vector<int> data_p_strata = data["strata"];

    std::vector<int> strata_t;

    for(int i = 0; i < data_p_waves.size() - 1; i++){
      std::vector<int> strata_ti;

      for(int j = 0; j < data_p_waves[i].size(); j ++){
        strata_ti.push_back(data_p_strata[data_p_waves[i][j] - 1]);
      }

      strata_t.insert(strata_t.end(), strata_ti.begin(), strata_ti.end());
    }

    std::vector<int> strat(strata_t);
    sort(strat.begin(),strat.end());
    strat.erase(unique(strat.begin(),strat.end() ),strat.end());

    std::vector<int> strata_count;

    for(int i = 0; i < strat.size(); i++){
      strata_count.push_back(std::count(strata_t.begin(),strata_t.end(), strat[i]));
    }

    //get p(no link between strata)
    std::vector<double> no_link_init(n_strata, 1);

    for(int i = 0; i < n_strata; i++){
      for(int j = 0; j < n_strata; j++){
        no_link_init[i] = no_link_init[i] * std::pow((1 - b[t-1][j][i]),strata_count[j]);
      }
    }

    std::vector<double> no_link_l;
    std::transform (l[t-1].begin(),l[t-1].end(),
                    no_link_init.begin(),
                    std::back_inserter(no_link_l),
                    std::multiplies<double>());

    double no_link = std::accumulate(no_link_l.begin(), no_link_l.end(), 0.0);

    int nn_0 = data_p_waves[0].size();
    int nn = std::accumulate(n_p.begin(),n_p.end(),0);

    std::vector<int> n_post_range = gen_range(nn, total * 5);

    std::vector<double> n_sample_prob_vec;

    for(int i = 0; i < n_post_range.size(); i++){
     std::vector<int> r_i = gen_range(n_post_range[i] + 1 - nn, n_post_range[i] - nn_0);

     std::vector<double> log_r_i;
     for(int j = 0; j < r_i.size(); j++){
       log_r_i.push_back(log(r_i[j]));
     }

     double s_log_r_i = std::accumulate(log_r_i.begin(), log_r_i.end(), 0.0);

     n_sample_prob_vec.push_back(
       s_log_r_i + (n_post_range[i] - nn) * log(no_link) - prior_n * log(n_post_range[i])
     );
    }

    double max_n_sample_prob_vec = *std::max_element(n_sample_prob_vec.begin(), n_sample_prob_vec.end());

    std::vector<double> n_sample_prob;

    for(int i = 0; i < n_sample_prob_vec.size(); i++){
      n_sample_prob.push_back(exp(n_sample_prob_vec[i] - max_n_sample_prob_vec));
    }

    std::vector<int> nt = as<std::vector<int>>(c_sample(n_post_range, 1, false, n_sample_prob));
    n[t] = nt[0];

    //#######################
    //# generate new lambda #
    //#######################

    // assign strata to non sampled units

    // get indices of non sampled units
    std::vector<int> n_range = gen_range(0,n[t] - 1);
    std::vector<int> not_sampled;
    std::vector<int> data_p_waves_id_sort = data_p_waves_id;
    sort(data_p_waves_id_sort.begin(),data_p_waves_id_sort.end());
    std::set_difference(n_range.begin(),n_range.end(),data_p_waves_id_sort.begin(), data_p_waves_id_sort.end(),
                        std::inserter(not_sampled, not_sampled.end()));

    std::vector<int> stratum;
    stratum.resize(n[t]);

    //fill stratum vector with strata of sampled units
    std::vector<int> ins_val_us;

    for(int i = 0; i < data_p_waves_id.size(); i++){
      ins_val_us.push_back(data_p_strata[data_p_waves_id[i]]);
    }

    stratum = int_vec_insert(stratum,
                             ins_val_us,
                             data_p_reorder_id);

    //fill stratum vector with strata of non sampled units
    std::vector<int> strat_s = gen_range(1,n_strata);
    std::vector<double> pstrat;

    for(int i = 0; i < no_link_l.size(); i++){
      pstrat.push_back(no_link_l[i]/no_link);
    }

    std::vector<int> stratsamp = as<std::vector<int>>(c_sample(strat_s, not_sampled.size(), true, pstrat));

    stratum = int_vec_insert(stratum,
                             stratsamp,
                             not_sampled);


    // populate link matrix for reordered sample
    std::vector<int> g1 = rep_times(gen_range(0,n_waves - 1), n_waves);
    std::vector<int> g2 = rep_each(gen_range(0, n_waves - 1), n_waves);

    arma::mat y = arma::mat(n[t],n[t]);

    for(int i = 0; i < g1.size() - 1; i++){
      y = mat_to_mat_insert(y_samp,
                            y,
                            data_p_reorder[g1[i]],
                            data_p_reorder[g2[i]],
                            data_p_waves[g1[i]],
                            data_p_waves[g2[i]]);
    }


    std::vector<int>

    IntegerVector lp_1 = c_unlist(data_p_reorder[Range(0,data_p_reorder.size() - 2)]);
    IntegerVector lp_2 = Range(1,n[t]);
    IntegerMatrix link_pairs =  c_combn(c_setdiff(lp_2,lp_1),2);
    link_pairs  = transpose(link_pairs);

    int n_pairs = link_pairs.nrow();
    NumericVector link_prob = runif(n_pairs);

    IntegerVector assigned = c_which(link_prob < mat_by_mat(b[t],
                                                            stratum[link_pairs(_,0)],
                                                            stratum[link_pairs(_,1)]));

    IntegerVector link_pairs_0 = link_pairs(_,0);
    IntegerVector link_pairs_1 = link_pairs(_,1);
    IntegerVector vals (link_pairs_0.length(),1);

    y = int_mat_insert(y,
                       link_pairs_1[assigned - 1],
                       link_pairs_0[assigned - 1],
                       vals);

    y = int_mat_insert(y,
                       link_pairs_0[assigned - 1],
                       link_pairs_0[assigned - 1],
                       vals);


    // new lambda
    NumericVector strata_count_num = as<NumericVector>(table(stratum));
    IntegerVector strate_count_int = table(stratum);
    NumericMatrix dirichlet = c_rdirichlet(1, strata_count_num + prior_l);
    l(t) = dirichlet(0,_);

    //#####################
    //# generate new beta #
    //#####################



  }

  List ret(3);
  ret[0] = data_p_reorder;
  ret[1] = data_p_waves;
  ret[2] = ins_pos_us;
  return ret;
}






