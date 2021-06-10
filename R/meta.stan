data {
  int<lower=0> N;
  int<lower=0> K;
  int<lower=0,upper=N> n1;
  int<lower=0,upper=N> n2;
  int<lower=0,upper=N> n3;
  int<lower=0,upper=N> n4;
  int<lower=0,upper=N> n5;
  int<lower=0,upper=N> n6;
  int<lower=0,upper=N> n7;
  int<lower=0,upper=N> n8;
  int<lower=0,upper=N> observed1[n1];
  int<lower=0,upper=N> observed2;
  int<lower=0,upper=N> observed3[n3];
  int<lower=0,upper=N> observed4[n4];
  int<lower=0,upper=N> observed5;
  int<lower=0,upper=N> observed6;
  int<lower=0,upper=N> observed7;
  int<lower=0,upper=N> observed8;
  real<lower=0> est1[n1];
  real<lower=0> est2;
  real<lower=0> est3[n3];
  real<lower=0> est4[n4];
  real<lower=0> est5;
  real<lower=0> est6;
  real<lower=0> est7;
  real<lower=0> est8;
  real<lower=0> est_se1[n1];
  real<lower=0> est_se2;
  real<lower=0> est_se3[n3];
  real<lower=0> est_se4[n4];
  real<lower=0> est_se5;
  real<lower=0> est_se6;
  real<lower=0> est_se7;
  real<lower=0> est_se8;
}
parameters {
  real error[K];
  vector<lower=0>[N] alpha;
}
model {
  target += normal_lpdf(est1| error[1] + alpha[observed1], est_se1);
  target += normal_lpdf(est2| error[2] + alpha[observed2], est_se2);
  target += normal_lpdf(est3| error[3] + alpha[observed3], est_se3);
  target += normal_lpdf(est4| error[4] + alpha[observed4], est_se4);
  target += normal_lpdf(est5| error[5] + alpha[observed5], est_se5);
  target += normal_lpdf(est6| error[6] + alpha[observed6], est_se6);
  target += normal_lpdf(est7| error[7] + alpha[observed7], est_se7);
  target += normal_lpdf(est8| error[8] + alpha[observed8], est_se8);
  error ~ normal(0, 100);
  alpha ~ normal(0, 100);
}
