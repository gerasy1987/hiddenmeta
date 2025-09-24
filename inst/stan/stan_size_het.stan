data {
  int<lower=0> N;
  int<lower=0> K;
  int<lower=0> n_ests;
  vector[N] covariate;
  real<lower=0> alpha_prior_mean[N];
  real<lower=0> alpha_prior_sd[N];
  real<lower=0> bias_a_prior_mean[K];
  real<lower=0> bias_a_prior_sd[K];
  real<lower=0> bias_b_prior_mean[K];
  real<lower=0> bias_b_prior_sd[K];
  int<lower=0> study[n_ests];
  int<lower=0> design[n_ests];
  vector<lower=0>[n_ests] ests;
  vector<lower=0>[n_ests] ses;
}
parameters {
  vector[K] bias_a;
  vector[K] bias_b;
  vector<lower=0>[N] alpha;
}
transformed parameters {
  vector<lower=0.01>[n_ests] bias;

  for (x in 1:n_ests)
      bias[x] =  exp(bias_a[design[x]]  +  bias_b[design[x]] * covariate[study[x]] - 1);

}
model {

  target += normal_lpdf(ests | bias .* alpha[study], ses);
  alpha ~ normal(alpha_prior_mean, alpha_prior_sd);
  bias_a ~ normal(bias_a_prior_mean, bias_a_prior_sd);
  bias_b ~ normal(bias_b_prior_mean, bias_b_prior_sd);

}
