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
  int<lower=0,upper=N> observed2[n2];
  int<lower=0,upper=N> observed3;
  int<lower=0,upper=N> observed4[n4];
  int<lower=0,upper=N> observed5;
  int<lower=0,upper=N> observed6;
  int<lower=0,upper=N> observed7;
  int<lower=0,upper=N> observed8;
  real<lower=0> est1[n1];
  real<lower=0> est2[n2];
  real<lower=0> est3;
  real<lower=0> est4[n4];
  real<lower=0> est5;
  real<lower=0> est6;
  real<lower=0> est7;
  real<lower=0> est8;
  real<lower=0> est_se1[n1];
  real<lower=0> est_se2[n2];
  real<lower=0> est_se3;
  real<lower=0> est_se4[n4];
  real<lower=0> est_se5;
  real<lower=0> est_se6;
  real<lower=0> est_se7;
  real<lower=0> est_se8;
}
parameters {
  real<lower=0.01> rel_bias[K-1];
  vector<lower=1>[N] alpha;
}
model {
  target += normal_lpdf(est1 | alpha[observed1], est_se1);
  target += normal_lpdf(est2 | rel_bias[1] * alpha[observed2], est_se2);
  target += normal_lpdf(est3 | rel_bias[2] * alpha[observed3], est_se3);
  target += normal_lpdf(est4 | rel_bias[3] * alpha[observed4], est_se4);
  target += normal_lpdf(est5 | rel_bias[4] * alpha[observed5], est_se5);
  target += normal_lpdf(est6 | rel_bias[5] * alpha[observed6], est_se6);
  target += normal_lpdf(est7 | rel_bias[6] * alpha[observed7], est_se7);
  target += normal_lpdf(est8 | rel_bias[7] * alpha[observed8], est_se8);
  rel_bias ~ normal(1, 10);
  alpha ~ normal(300, 100);
}
generated quantities {
 matrix<lower=0>[K,K] rel_biases;
 vector<lower=0>[K] rel_bias_full;

 rel_bias_full[1] = 1;

 for (i in 2:(K - 1)) {
   rel_bias_full[i] = rel_bias[i - 1];
 }


 for (i in 1:K) {
   for (j in 1:K) {
     rel_biases[i,j] = rel_bias_full[i]/rel_bias_full[j];
   }
 }

}
