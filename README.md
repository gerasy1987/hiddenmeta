
# hiddenmeta ![R](https://github.com/gerasy1987/hiddenmeta/workflows/R/badge.svg) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://github.com/gerasy1987/hiddenmeta/blob/main/LICENSE)

## Overview

``` r
if (!require(pacman)) install.packages("pacman")

pacman::p_load_current_gh("gerasy1987/hiddenmeta")

pacman::p_load(
  tidyverse, knitr, magrittr, # data management/plotting/printing
  DeclareDesign, # Design declaration
  # rstan, # Stan
  RcppAlgos, # speed up combinatorics
  networkreporting, NSUM, # NSUM methods
  sspse, RDS, chords # RDS+ methods
)

# number of studies
N_studies <- 2

## STUDY 1
study_1 <- 
  list(
    # total population size for one study
    N = 1000,
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
    add_groups = list(p_service = 0.3, 
                      loc_1 = 0.3, loc_2 = 0.1, loc_3 = 0.2),
    
    # RDS parameters
    n_seed = 5,
    n_coupons = 3,
    target_type = "sample",
    target_n_rds = 50,
    
    # prop sampling parameters
    target_n_prop = 400,
    
    # TLS sampling parameters
    target_n_tls = 1
  )


## STUDY 2
study_2 <- 
  list(
    N = 2000,
    K = 3,
    prev_K = c(known_1 = .2, known_2 = .2, hidden = .1),
    rho_K = c(.05, .01, .01),
    p_edge_within = list(known_1 = c(0.05, 0.05),
                         known_2 = c(0.05, 0.05),
                         hidden = c(0.05, 0.5)),
    p_edge_between = list(known_1 = 0.05, 
                          known_2 = 0.05, 
                          hidden = 0.01),
    p_visibility = list(known_1 = .99,
                        known_2 = .99,
                        hidden = .5),
    add_groups = list(service_use = 0.2, 
                      loc_1 = 0.6, loc_2 = 0.5, loc_3 = 0.4),
    
    # RDS parameters
    n_seed = 20,
    n_coupons = 3,
    target_type = "sample",
    target_n_rds = 120,
    
    # prop sampling parameters
    target_n_prop = 800,
    
    # TLS sampling parameters
    target_n_tls = 2
  )


# SINGLE STUDY --------------------------------------------------------------------------------


population_study <-
  do.call(what = declare_population,
          args = c(handler = get_study_population, study_1[1:8]))

# population_study()

# population_study <- 
#   do.call(what = declare_population, 
#           args = c(handler = get_study_population, study_2[1:8]))
# 
# # population_study()

rds_study <- 
  do.call(declare_sampling,
          c(handler = sample_rds, 
            sampling_variable = "rds",
            drop_nonsampled = FALSE, study_2[9:12]))

# rds_study(population_study())

prop_study <- 
  do.call(declare_sampling,
          c(handler = sample_prop, 
            sampling_variable = "prop",
            drop_nonsampled = FALSE, study_2[13]))

# rds_study(prop_study(population_study()))

tls_study <- 
  do.call(declare_sampling,
          c(handler = sample_tls, 
            sampling_variable = "tls",
            drop_nonsampled = FALSE, study_2[14]))

# tls_study(rds_study(prop_study(population_study())))

draw_data(population_study + rds_study + prop_study + tls_study)


estimands <- 
  declare_estimand(handler = get_study_estimands)

set.seed(19872312)
draw_estimands(population_study + rds_study + prop_study + tls_study + estimands)
draw_estimands(population_study + estimands)
```
