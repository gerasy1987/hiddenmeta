data {
    int<lower=0> N;   // number of studies
    int<lower=0> K;   // max number of estimator (estimator-sampling pairs)
    // number of studies with estimator (estimator-sampling pair)
    int<lower=0,upper=N> N1;
    int<lower=0,upper=N> N2;
    int<lower=0,upper=N> N3;
    // ids of studies with specific estimator (estimator-sampling pair)
    int<lower=0,upper=N> observed1[N1];
    int<lower=0,upper=N> observed2[N2];
    int<lower=0,upper=N> observed3[N3];
    // parameter estimates
    real<lower=0,upper=1> est1[N1];
    real<lower=0,upper=1> est2[N2];
    real<lower=0,upper=1> est3[N3];
    // estimated standard errors of parameter estimates
    real<lower=0> est1_sd[N1];
    real<lower=0> est2_sd[N2];
    real<lower=0> est3_sd[N3];
  }
  parameters {
    // (additive) error factor for each estimator/estimator-sampling pair
    real<lower=-1,upper=1> error[K];
    // prevalence estimate for each study
    vector<lower=0,upper=1>[N] alpha;
    // need to add Sigma to allow for interdependence of errors across estimators
    // or studies
  }

  model
    target += normal_lpmf(est1 | error[1] + alpha[observed1], est1_sd);
    target += normal_lpmf(est2 | error[2] + alpha[observed2], est2_sd);
    target += normal_lpmf(est3 | error[3] + alpha[observed3], est3_sd);
  }
