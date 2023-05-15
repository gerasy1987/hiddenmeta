data {
  int<lower=0> N;
  int<lower=0> K;
  int<lower=0> n_ests;
  real<lower=0> alpha_prior_mean[N];
  real<lower=0> alpha_prior_sd[N];
  real<lower=0> bias_prior_mean[K];
  real<lower=0> bias_prior_sd[K];
  int<lower=1> study[n_ests];
  int<lower=1> design[n_ests];
  vector<lower=0>[n_ests] ests;
  vector<lower=0>[n_ests] ses;
}
parameters {
  vector<lower=0.01>[K] bias;
  vector<lower=0>[N] alpha;
}

model {

  target += normal_lpdf(ests | bias[design].*alpha[study], ses);
  alpha ~ normal(alpha_prior_mean, alpha_prior_sd);
  bias ~ lognormal(log(bias_prior_mean), bias_prior_sd);

}
