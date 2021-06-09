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
  real<lower=0> est2[n2];
  real<lower=0> est3[n3];
  real<lower=0> est4[n4];
  real<lower=0> est5[n5];
  real<lower=0> est6[n6];
  real<lower=0> est7[n7];
  real<lower=0> est8[n8];
  real<lower=0> est1_se[n1];
  real<lower=0> est2_se[n2];
  real<lower=0> est3_se[n3];
  real<lower=0> est4_se[n4];
  real<lower=0> est5_se[n5];
  real<lower=0> est6_se[n6];
  real<lower=0> est7_se[n7];
  real<lower=0> est8_se[n8];
}
parameters {
  real error[K];
  vector<lower=0>[N] alpha;
}
model {
  target += normal_lpdf(est1| error[1] + alpha[observed1], est1_se);
  target += normal_lpdf(est2| error[2] + alpha[observed2], est1_se);
  target += normal_lpdf(est3| error[3] + alpha[observed3], est1_se);
  target += normal_lpdf(est4| error[4] + alpha[observed4], est1_se);
  target += normal_lpdf(est5| error[5] + alpha[observed5], est1_se);
  target += normal_lpdf(est6| error[6] + alpha[observed6], est1_se);
  target += normal_lpdf(est7| error[7] + alpha[observed7], est1_se);
  target += normal_lpdf(est8| error[8] + alpha[observed8], est1_se);
  error ~ normal(0, 100);
  alpha ~ normal(0, 100);
}
