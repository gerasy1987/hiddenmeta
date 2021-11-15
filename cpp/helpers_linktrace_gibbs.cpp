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



List lt_gibbs() {

}



