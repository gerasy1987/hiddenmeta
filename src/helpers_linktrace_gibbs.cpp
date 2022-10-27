#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

//' cpp helper to move vector elements to new indices
//'
//' @param x integer vector of elements to be shuffled
//' @param old_index integer original index of element to be moved
//' @param new_index integer index element should be moved to
//' @return shuffled vector
//' @keywords internal
// [[Rcpp::export]]
std::vector<int> move_elements(std::vector<int> x, int old_index, int new_index){

  if(old_index > new_index){
    std::rotate(x.rend() - old_index - 1, x.rend() - old_index, x.rend() - new_index);
  } else {
    std::rotate(x.begin() + old_index, x.begin() + old_index + 1, x.begin() + new_index + 1);
  }

  return x;
}

//' n choose k helper for combn_cpp
//'
//' @param n integer number of elements to choose from
//' @param k integer number of elements to choose
//' @return integer number of ways to choose k items out of n
//' @keywords internal
// [[Rcpp::export]]
uint64_t choose_cpp(uint64_t n, uint64_t k) {
  if(k == 0) return 1;
  return (n * choose_cpp(n - 1, k - 1)) / k;
}

//' helper to generate combinations of x k at a time (cpp implementation of combn)
//'
//' @param x integer vector of elements to combine
//' @param K integer order of combinations
//' @return a matrix of k-wise combinations of elements in x
//' @keywords internal
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

//' helper to generate single draw from dirichlet distribution
//'
//' @param alpha double vector of alpha parameters
//' @return a double vector containing a single draw from a dirichlet distribution
//' @keywords internal
// [[Rcpp::export]]
std::vector<double> rdirichlet_cpp(std::vector<double> alpha){

  std::vector<double> vec(alpha.size());

  double sum = 0.0;

  for(int i = 0; i < alpha.size(); i++){

    if(alpha[i] > 0.0){
      vec[i] = as<double>(Rcpp::rgamma(1,alpha[i],(1/alpha[i])));
    } else {
      vec[i] = 0.0;
    }
    sum += vec[i];
  }

  for(int i = 0; i < vec.size(); i++){
    vec[i] = vec[i]/sum;
  }

  return vec;
}

//' cpp implementation of R table
//'
//' @param x an integer vector of elements to be counted
//' @return an integer vector of counts of unique elements in x (sorted in ascending order of elements in x)
//' @keywords internal
// [[Rcpp::export]]
std::vector<int> table_cpp(std::vector<int> x){

  std::vector<int> vec(x);
  sort(vec.begin(),vec.end());
  vec.erase(unique(vec.begin(),vec.end()),vec.end());

  std::vector<int> vec_count;

  for(int i = 0; i < vec.size(); i++){
    vec_count.push_back(std::count(x.begin(),x.end(), vec[i]));
  }

  return vec_count;
}

//' cpp implementation of standard R rep
//'
//' @param x integer vector to be repeated
//' @param n integer number of repetitions
//' @return integer vector x repeated n times
//' @keywords internal
// [[Rcpp::export]]
std::vector<int> rep_times(std::vector<int> x, int n){

  std::vector<int> ret;

  for(int i = 0; i < n; i++){
    ret.insert(ret.end(),x.begin(),x.end());
  }

  return ret;
}

//' cpp implementation of R rep with each argument
//'
//' @param x integer vector of elements to be repeated
//' @param n integer number of repetitions
//' @return integer vector with each element in x repeated n times in order
//' @keywords internal
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

//' helper to generate range of consecutive integers
//'
//' @param from integer first integer
//' @param to last integer
//' @return integer vector of consecutive integers between from and to inclusive
//' @keywords internal
// [[Rcpp::export]]
std::vector<int> gen_range(int from,
                           int to){

  int l = (to - from) + 1;
  std::vector<int> vec(l);
  std::iota(vec.begin(),vec.end(),from);
  return vec;

}

//' helper to do random access insertion for vectors
//'
//' @param vec integer vector to insert values into
//' @param vals integer vector of values to insert into vec
//' @param pos integer vector of index positions in vec to insert vals into
//' @return integer vector with vals inserted at pos
//' @keywords internal
// [[Rcpp::export]]
std::vector<int> int_vec_insert(std::vector<int> vec,
                                std::vector<int> vals,
                                std::vector<int> pos){

  for(int i = 0; i < vals.size(); i++){
    vec[pos[i]] = vals[i];
  }

  return vec;
}


//' helper to do random access insertion of matrix into another matrix
//'
//' @param old_m matrix to take values from
//' @param new_m matrix to insert values into
//' @param new_rows integer vector of row indeces of new_m to insert values into
//' @param new_cols integer vector of column indeces of new_m to insert values into
//' @param old_rows integer vector of row indeces of old_m to take values from
//' @param old_cols integer vector of column indeces of old_m to take values from
//' @return matrix new_m with values from old_m in specified positions
//' @keywords internal
// [[Rcpp::export]]
arma::mat mat_to_mat_insert(arma::mat old_m,
                            arma::mat new_m,
                            std::vector<int> new_rows,
                            std::vector<int> new_cols,
                            std::vector<int> old_rows,
                            std::vector<int> old_cols){

  for(int i = 0; i < new_rows.size(); i++){
    for(int j = 0; j < new_cols.size(); j++){
      new_m(new_rows[i] - 1 , new_cols[j] - 1) = old_m(old_rows[i] - 1 , old_cols[j] - 1);
    }
  }

  return new_m;
}


//' helper to permute sampling data
//'
//' @param link_list List holding indices of linked units for each unit
//' @param wave integer vector holding the sampling wave of each unit
//' @param name integer vector holding the name of each unit
//' @return vector of integer vectors holding permuted sampling waves
//' @keywords internal
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


//' Link-tracing Gibbs sampler
//' @param data DataFrame pass through of sammple data created in get_study_est_linktrace
//' @param y_samp matrix pass through of adjecency matrix generated in get_study_est_linktrace
//' @param strata integer stratum id of each unit
//' @param n_strata integer number of unique strata
//' @param n_waves integer number of sampling waves
//' @param total integer total size of the population
//' @param chain_samples integer number of samples per MCMC chain
//' @param priors List of priors for n, lambda and beta
//' @param param_init List of initial values for n,lambda and beta
//' @return integer vector of samples from the n distribution
//' @keywords internal
// [[Rcpp::export]]
std::vector<int> lt_gibbs(DataFrame data,
                          arma::mat y_samp,
                          std::vector<int> strata,
                          int n_strata,
                          int n_waves,
                          int total,
                          int chain_samples,
                          List priors,
                          List param_init) {

  Function cpp_sample("sample");

  // permute data
  std::vector<std::vector<int>> data_p_waves;
  int last_wave = 0;

  while(last_wave < 1){

    data_p_waves = lt_permute(data["links_list"],data["rds_wave"],data["name"]);
    last_wave = data_p_waves[n_waves - 1].size();

  }

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
  arma::mat l(chain_samples,n_strata);
  arma::cube b(n_strata,n_strata,chain_samples);
  std::vector<int> n;

  l.row(0) = arma::conv_to<arma::rowvec>::from(as<std::vector<double>>(param_init["l_0"]));
  b.slice(0) = as<arma::mat>(param_init["b_0"]);
  n.push_back(as<int>(param_init["n_0"]));

  int prior_n = priors["p_n"];
  std::vector<double> prior_l = priors["p_l"];
  int prior_b = priors["p_b"];


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
  std::vector<int> strata_count = table_cpp(strata_t);

  // begin MCMC
  for(int t = 1; t < chain_samples; t++){

    //##################
    //# generate new N #
    //##################

    //get p(no link between strata)
    std::vector<double> no_link_init(n_strata, 1);

    for(int i = 0; i < n_strata; i++){
      for(int j = 0; j < n_strata; j++){
        no_link_init[i] = no_link_init[i] * std::pow((1 - b.slice(t-1)(j,i)),strata_count[j]);
      }
    }

    std::vector<double> no_link_l;
    std::vector<double> lt = arma::conv_to<std::vector<double>>::from(l.row(t-1));
    std::transform (lt.begin(),lt.end(),
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

    std::vector<int> nt = as<std::vector<int>>(cpp_sample(n_post_range, 1, false, n_sample_prob));
    n.push_back(nt[0]);

    //#######################
    //# generate new lambda #
    //#######################

    // assign strata to non sampled units

    // get indices of non sampled units
    std::vector<int> n_range = gen_range(0,n[t] - 1);
    std::vector<int> not_sampled;
    std::set_difference(n_range.begin(),n_range.end(),data_p_reorder_id.begin(), data_p_reorder_id.end(),
                        std::inserter(not_sampled, not_sampled.end()));

    std::vector<int> stratum(n[t]);

    // fill stratum vector with strata of sampled units
    std::vector<int> ins_val_us;

    for(int i = 0; i < data_p_waves_id.size(); i++){
      ins_val_us.push_back(data_p_strata[data_p_waves_id[i]]);
    }

    stratum = int_vec_insert(stratum,
                             ins_val_us,
                             data_p_reorder_id);

    // fill stratum vector with strata of non sampled units
    std::vector<int> strat_s = gen_range(1,n_strata);
    std::vector<double> pstrat;

    for(int i = 0; i < no_link_l.size(); i++){
      pstrat.push_back(no_link_l[i]/no_link);
    }

    std::vector<int> stratsamp = as<std::vector<int>>(cpp_sample(strat_s, not_sampled.size(), true, pstrat));

    stratum = int_vec_insert(stratum,
                             stratsamp,
                             not_sampled);


    // populate link matrix for reordered sample
    // generate empty matrix and fill known link pairs
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

    //generate unkown pairs
    std::vector<int> lp_1;
    for(int i = 0; i < data_p_reorder.size() - 1; i++){
      lp_1.insert(lp_1.end(), data_p_reorder[i].begin(), data_p_reorder[i].end());
    }
    sort(lp_1.begin(),lp_1.end());

    std::vector<int> lp_2 = gen_range(1,n[t]);

    std::vector<int> lp;
    std::set_difference(lp_2.begin(),lp_2.end(),lp_1.begin(),lp_1.end(),
                        std::inserter(lp,lp.end()));


    // if an unknown pair exists add links based on link probability
    if(lp.size() > 1){

      arma::mat link_pairs =  combn_cpp(lp,2);
      int n_pairs = link_pairs.n_rows;

      std::vector<double> link_prob(n_pairs);
      std::generate(link_prob.begin(), link_prob.end(), [](){
        return (double)std::rand() / (double)RAND_MAX;
      });

      for(int i = 0; i < n_pairs; i++){

        int id_1 = link_pairs(i,0) - 1;
        int id_2 = link_pairs(i,1) - 1;

        double link_prob_i = b.slice(t - 1)(stratum[id_1] - 1, stratum[id_2] - 1);

        if(link_prob_i > link_prob[i]){
          y(id_1,id_2) = 1;
          y(id_2,id_1) = 1;
        }

      }

    }

    // new lambda
    std::vector<int> strata_count_int = table_cpp(stratum);

    // if a certain stratum was not sampled we need to add 0 to the count
    if(strata_count_int.size() < n_strata){

      std::vector<int> s = gen_range(1,n_strata);
      std::vector<int> s_new(stratum);
      sort(s_new.begin(),s_new.end());
      std::vector<int> s_miss;
      std::set_difference(s.begin(),s.end(),s_new.begin(),s_new.end(),
                          std::inserter(s_miss, s_miss.end()));

      for(int i = 0; i < s_miss.size(); i++){
        strata_count_int.push_back(0);
      }

      for(int i = 0; i < s_miss.size(); i++){
        strata_count_int = move_elements(strata_count_int,
                                         strata_count_int.size() - (s_miss.size() - i),
                                         s_miss[i] - 1);
      }
    }

    std::vector<double> alphas;
    std::transform(strata_count_int.begin(), strata_count_int.end(),
                   prior_l.begin(),std::back_inserter(alphas),
                   std::plus<double>());

    l.row(t) = arma::conv_to<arma::rowvec>::from(rdirichlet_cpp(alphas));

    //#####################
    //# generate new beta #
    //#####################

    // count links between strata
    arma::mat strata_link_count(n_strata,n_strata);
    arma::mat node_pairs = combn_cpp(gen_range(1,n[t]),2);
    int n_pairs_b = node_pairs.n_rows;

    for(int i = 0; i < n_pairs_b; i++){

      int np_1 = node_pairs(i,0) - 1;
      int np_2 = node_pairs(i,1) - 1;

      int c_1 = y(np_1,np_2);
      int c_2 = y(np_2,np_1);

      int stratum_1 = stratum[np_1] - 1;
      int stratum_2 = stratum[np_2] - 1;

      strata_link_count(stratum_1,stratum_2) = strata_link_count(stratum_1,stratum_2) + c_1;
      strata_link_count(stratum_2,stratum_1) = strata_link_count(stratum_2,stratum_1) + c_2;

    }

    arma::mat b_i(n_strata,n_strata);

    for(int i = 0; i < n_strata; i++){

      double shape_1 = strata_link_count(i,i) + prior_b;
      double shape_2 = choose_cpp(strata_count_int[i],2) - strata_link_count(i,i) + prior_b;
      b_i(i,i) = as<double>(Rcpp::rbeta(1,shape_1,shape_2));

    }

    for(int i = 0; i < n_strata - 1; i++){
      for(int j = i + 1; j < n_strata; j++){

        double shape_1 = strata_link_count(i,j) + strata_link_count(j,i) + prior_b;
        double shape_2 = strata_count_int[i] * strata_count_int[j] - strata_link_count(i,j) - strata_link_count(j,i) + prior_b;
        double beta_i = as<double>(Rcpp::rbeta(1,shape_1,shape_2));

        b_i(i,j) = beta_i;
        b_i(j,i) = beta_i;

      }
    }

    b.slice(t) = b_i;
  }

  return n;
}



//[[Rcpp::export]]

std::vector<int> lt_gibbs_test(DataFrame data,
                               arma::mat y_samp,
                               std::vector<int> strata,
                               int n_strata,
                               int n_waves,
                               int total,
                               int chain_samples,
                               List priors,
                               List param_init) {

  Function cpp_sample("sample");

  //********************************
  //* Assign priors and inital vals*
  //********************************
  arma::mat l(chain_samples,n_strata);
  arma::cube b(n_strata,n_strata,chain_samples);
  std::vector<int> n;

  l.row(0) = arma::conv_to<arma::rowvec>::from(as<std::vector<double>>(param_init["l_0"]));
  b.slice(0) = as<arma::mat>(param_init["b_0"]);
  n.push_back(as<int>(param_init["n_0"]));

  int prior_n = priors["p_n"];
  std::vector<double> prior_l = priors["p_l"];
  int prior_b = priors["p_b"];

  //***************
  //* permute data*
  //***************
  std::vector<std::vector<int>> data_p_waves;
  int last_wave = 0;

  while(last_wave < 1){
    data_p_waves = lt_permute(data["links_list"],data["rds_wave"],data["name"]);
    last_wave = data_p_waves[n_waves - 1].size();
  }

  std::vector<int> data_p_waves_id;

  for(int i = 0; i < data_p_waves.size(); i++){
    for(int j = 0; j < data_p_waves[i].size(); j++){
      data_p_waves_id.push_back(data_p_waves[i][j] - 1);
    }
  }

  //*****************
  //* re-index units*
  //*****************
  std::vector<int> n_p(n_waves);

  for(int i = 0; i < n_waves; i++){
    n_p[i] = data_p_waves[i].size();
  }

  std::vector<std::vector<int>>  data_p_reorder(n_waves);
  data_p_reorder[0] = gen_range(1,n_p[0]);

  for(int i = 1; i < n_waves; i++){
    int from_i = std::accumulate(n_p.begin(), n_p.begin() + i, 0) + 1;
    int to_i = std::accumulate(n_p.begin(), n_p.begin() + i + 1, 0);
    data_p_reorder[i] = gen_range(from_i,to_i);
  }

  std::vector<int> data_p_reorder_id;

  for(int i = 0; i < data_p_reorder.size(); i++){
    for(int j = 0; j < data_p_reorder[i].size(); j++){
      data_p_reorder_id.push_back(data_p_reorder[i][j] - 1);
    }
  }

  //*****************************************************
  //* for wave 1:n-1 get number of units in each stratum*
  //*****************************************************
  std::vector<int> data_p_strata = data["strata"];
  int first_n = std::accumulate(n_p.begin(), n_p.end() - 1, 0);
  std::vector<int> strata_t(first_n);

  for(int i = 0; i < first_n; i++){
    strata_t[i] = data_p_strata[data_p_waves_id[i]];
  }

  std::vector<int> strata_count = table_cpp(strata_t);

  //******************************
  //* get strata of sampled units*
  //******************************
  std::vector<int> stratum_sampled(data_p_waves_id.size());

  for(int i = 0; i < data_p_waves_id.size(); i++){
    stratum_sampled[i] = data_p_strata[data_p_waves_id[i]];
  }

  //***********************************************
  //* fill reo-rdered link matrix with known pairs*
  //***********************************************
  int n_units = std::accumulate(n_p.begin(), n_p.end(), 0);
  arma::mat y_known(n_units,n_units);

  std::vector<int> g1 = rep_times(gen_range(0,n_waves - 1), n_waves);
  std::vector<int> g2 = rep_each(gen_range(0, n_waves - 1), n_waves);

  for(int i = 0; i < g1.size() - 1; i++){
    y_known = mat_to_mat_insert(y_samp,
                                y_known,
                                data_p_reorder[g1[i]],
                                data_p_reorder[g2[i]],
                                data_p_waves[g1[i]],
                                data_p_waves[g2[i]]);
  }


  // begin MCMC
  for(int t = 1; t < chain_samples; t++){

    //*****************
    //* generate new N*
    //*****************

    //get p(no link between strata)
    std::vector<double> no_link_init(n_strata, 1);

    for(int i = 0; i < n_strata; i++){
      for(int j = 0; j < n_strata; j++){
        no_link_init[i] = no_link_init[i] * std::pow((1 - b.slice(t-1)(j,i)),strata_count[j]);
      }
    }

    std::vector<double> no_link_l;
    std::vector<double> lt = arma::conv_to<std::vector<double>>::from(l.row(t-1));
    std::transform (lt.begin(),lt.end(),
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

    std::vector<int> nt = as<std::vector<int>>(cpp_sample(n_post_range, 1, false, n_sample_prob));
    n.push_back(nt[0]);

    //**********************
    //* generate new lambda*
    //**********************

    // assign strata to non sampled units
    // get indices of non sampled units
    std::vector<int> not_sampled = gen_range(data_p_reorder_id.back() + 1,n[t] - 1);

    //fill stratum vector with strata of sampled units
    std::vector<int> stratum(stratum_sampled);

    // fill stratum vector with strata of non sampled units
    std::vector<int> strat_s = gen_range(1,n_strata);
    std::vector<double> pstrat;

    for(int i = 0; i < no_link_l.size(); i++){
      pstrat.push_back(no_link_l[i]/no_link);
    }

    std::vector<int> stratsamp = as<std::vector<int>>(cpp_sample(strat_s, not_sampled.size(), true, pstrat));

    for(int i = 0; i < stratsamp.size(); i++){
      stratum.push_back(stratsamp[i]);
    }

    // populate link matrix for reordered sample
    arma::mat y(y_known);

    if(n[t] > y.n_rows){
      y.resize(n[t],n[t]);
    }

    //generate unkown pairs
    std::vector<int> lp_1(data_p_reorder_id.begin(), data_p_reorder_id.begin() + first_n);
    std::vector<int> lp_2 = gen_range(0,n[t] - 1);

    std::vector<int> lp;
    std::set_difference(lp_2.begin(),lp_2.end(),lp_1.begin(),lp_1.end(),
                        std::inserter(lp,lp.end()));

    // if an unknown pair exists add links based on link probability
    if(lp.size() > 1){

      arma::mat link_pairs =  combn_cpp(lp,2);
      int n_pairs = link_pairs.n_rows;

      std::vector<double> link_prob(n_pairs);
      std::generate(link_prob.begin(), link_prob.end(), [](){
        return (double)std::rand() / (double)RAND_MAX;
      });

      for(int i = 0; i < n_pairs; i++){

        int id_1 = link_pairs(i,0);
        int id_2 = link_pairs(i,1);

        double link_prob_i = b.slice(t - 1)(stratum[id_1] - 1, stratum[id_2] - 1);

        if(link_prob_i > link_prob[i]){
          y(id_1,id_2) = 1;
          y(id_2,id_1) = 1;
        }

      }

    }

    // new lambda
    std::vector<int> strata_count_int = table_cpp(stratum);

    // if a certain stratum was not sampled we need to add 0 to the count
    if(strata_count_int.size() < n_strata){

      std::vector<int> s = gen_range(1,n_strata);
      std::vector<int> s_new(stratum);
      sort(s_new.begin(),s_new.end());
      std::vector<int> s_miss;
      std::set_difference(s.begin(),s.end(),s_new.begin(),s_new.end(),
                          std::inserter(s_miss, s_miss.end()));

      for(int i = 0; i < s_miss.size(); i++){
        strata_count_int.push_back(0);
      }

      for(int i = 0; i < s_miss.size(); i++){
        strata_count_int = move_elements(strata_count_int,
                                         strata_count_int.size() - (s_miss.size() - i),
                                         s_miss[i] - 1);
      }
    }

    std::vector<double> alphas;
    std::transform(strata_count_int.begin(), strata_count_int.end(),
                   prior_l.begin(),std::back_inserter(alphas),
                   std::plus<double>());

    l.row(t) = arma::conv_to<arma::rowvec>::from(rdirichlet_cpp(alphas));

    //#####################
    //# generate new beta #
    //#####################

    // count links between strata
    arma::mat strata_link_count(n_strata,n_strata);
    arma::mat node_pairs = combn_cpp(gen_range(1,n[t]),2);
    int n_pairs_b = node_pairs.n_rows;

    for(int i = 0; i < n_pairs_b; i++){

      int np_1 = node_pairs(i,0) - 1;
      int np_2 = node_pairs(i,1) - 1;

      int c_1 = y(np_1,np_2);
      int c_2 = y(np_2,np_1);

      int stratum_1 = stratum[np_1] - 1;
      int stratum_2 = stratum[np_2] - 1;

      strata_link_count(stratum_1,stratum_2) = strata_link_count(stratum_1,stratum_2) + c_1;
      strata_link_count(stratum_2,stratum_1) = strata_link_count(stratum_2,stratum_1) + c_2;

    }

    arma::mat b_i(n_strata,n_strata);

    for(int i = 0; i < n_strata; i++){

      double shape_1 = strata_link_count(i,i) + prior_b;
      double shape_2 = choose_cpp(strata_count_int[i],2) - strata_link_count(i,i) + prior_b;
      b_i(i,i) = as<double>(Rcpp::rbeta(1,shape_1,shape_2));

    }

    for(int i = 0; i < n_strata - 1; i++){
      for(int j = i + 1; j < n_strata; j++){

        double shape_1 = strata_link_count(i,j) + strata_link_count(j,i) + prior_b;
        double shape_2 = strata_count_int[i] * strata_count_int[j] - strata_link_count(i,j) - strata_link_count(j,i) + prior_b;
        double beta_i = as<double>(Rcpp::rbeta(1,shape_1,shape_2));

        b_i(i,j) = beta_i;
        b_i(j,i) = beta_i;

      }
    }

    b.slice(t) = b_i;
  }

  return n;
}





