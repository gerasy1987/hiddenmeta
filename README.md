
# hiddenmeta ![R](https://github.com/gerasy1987/hiddenmeta/workflows/R/badge.svg) [![codecov](https://codecov.io/gh/gerasy1987/hiddenmeta/branch/main/graph/badge.svg?token=ZG9A64Q0A1)](https://codecov.io/gh/gerasy1987/hiddenmeta) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://github.com/gerasy1987/hiddenmeta/blob/main/LICENSE)

## Overview

### Install and load packages

``` r
install.packages("DeclareDesign")
install.packages("igraph")

devtools::install_github("gerasy1987/hiddenmeta", build_vignettes = TRUE)
```

``` r
library(hiddenmeta)
library(DeclareDesign)
library(igraph)
```

### Step 1. Provide study design features

``` r
## STUDY 1
study_1 <- 
  list(
    # total population size for one study
    N = 2000,
    # number of groups
    # (K-th group is hidden population we are interested in)
    K = 2,
    # probability of membership in each of the groups (prev_K[K] is the true prevalence)
    prev_K = c(known = .3, hidden = .1),
    # correlation matrix of group memberships
    rho_K = .05,
    
    # block edge probabilities depending on group memberships
    # 1 - list of in- and out-group probability of links for each group
    # 2 - probability of link between in- and out-group members
    p_edge_within = list(known = c(0.05, 0.05), hidden = c(0.05, 0.9)),
    p_edge_between = list(known = 0.05, hidden = 0.01),
    
    # probability of visibility (show-up) for each group
    p_visibility = list(known = .99, hidden = .7),
    
    # probability of service utilization in hidden population
    # for service multiplier
    add_groups = list(service_use = 0.3, 
                      loc_1 = 0.3, loc_2 = 0.1, loc_3 = 0.2, 
                      known_2 = 0.1, known_3 = 0.2),
    
    # RDS parameters
    n_seed = 20,
    n_coupons = 3,
    target_type = "sample",
    target_n_rds = 100,
    
    # prop sampling parameters
    target_n_pps = 400,
    
    # TLS sampling parameters
    target_n_tls = 1
  )
```

### Step 2. Declare study population

``` r
population_study <-
  do.call(what = declare_population,
          args = c(handler = get_study_population, study_1[1:8]))

set.seed(19872312)
( example_pop <- population_study() )
```

### Step 3. Declare all relevant study sampling procedures

The sampling procedures are additive in a sense that each procedure
appends several columns relevant to the sampling procedure and
particular draw based on population simulation, but does not change the
study population data frame (unless you specify
`drop_nonsampled = TRUE`).

``` r
rds_study <- 
  do.call(declare_sampling,
          c(handler = sample_rds, 
            sampling_variable = "rds",
            drop_nonsampled = FALSE, study_1[9:12]))

set.seed(19872312)
draw_data(population_study + 
            rds_study)
```

``` r
pps_study <- 
  do.call(declare_sampling,
          c(handler = sample_pps, 
            sampling_variable = "pps",
            drop_nonsampled = FALSE, study_1[13]))

set.seed(19872312)
draw_data(population_study + 
            rds_study + pps_study)
```

``` r
tls_study <- 
  do.call(declare_sampling,
          c(handler = sample_tls, 
            sampling_variable = "tls",
            drop_nonsampled = FALSE, study_1[14]))

set.seed(19872312)
draw_data(population_study + 
            rds_study + pps_study + tls_study)
```

### Step 4. Declare study level estimands

``` r
study_estimands <- 
  declare_estimand(handler = get_study_estimands)

set.seed(19872312)
draw_estimands(population_study + 
                 rds_study + pps_study + tls_study + 
                 study_estimands)
```

### Step 5. Declare estimators used in the study

``` r
estimator_sspse <- declare_estimator(handler = get_study_est_sspse, label = "sspse")
estimator_ht <- declare_estimator(handler = get_study_est_ht, label = "ht")
estimator_chords <- declare_estimator(type = "integrated",
                                      handler = get_study_est_chords, label = "chords")
estimator_nsum <- declare_estimator(handler = get_study_est_nsum, label = "nsum")

set.seed(19872312)
draw_estimates(population_study +
                 rds_study + pps_study + tls_study +
                 study_estimands +
                 estimator_sspse + estimator_ht + estimator_chords + estimator_nsum)
```

### Step 6. Diagnose study design

``` r
study_diagnosands <-
  declare_diagnosands(
    mean_estimand = mean(estimand),
    mean_estimate = mean(estimate),
    sd_estimate = sd(estimate),
    bias = mean(estimate - estimand),
    rmse = sqrt(mean((estimate - estimand) ^ 2))
  )

diagnose_design(
  population_study + 
    rds_study + pps_study + 
    study_estimands + 
    estimator_sspse + estimator_ht + estimator_chords + estimator_nsum, 
  diagnosands = study_diagnosands,
  sims = 10,
  bootstrap_sims = 10)
```
