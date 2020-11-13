
# hiddenmeta ![R](https://github.com/gerasy1987/hiddenmeta/workflows/R/badge.svg) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://github.com/gerasy1987/hiddenmeta/blob/main/LICENSE)

## Overview

``` r
if (!require(pacman)) install.packages("pacman")

pacman::p_load(
  tidyverse, knitr, magrittr, # data management/plotting/printing
  DeclareDesign, # Design declaration
  # rstan, # Stan
  RcppAlgos, # speed up combinatorics
  bindata, # MV binomial dsitribution
  networkreporting, NSUM, # NSUM methods 
  sspse, RDS, chords, # RDS+ methods
  statnet, network, sna, # network packages 
  igraph, ggnetwork, 
  rgl, htmltools # network visualization
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
    p_visibility = list(hidden = 1, known = 1),
    
    # probability of service utilization in hidden population
    # for service multiplier
    p_service = 0.3,
    
    # rds parameters
    n_seed_hidden = 5,
    n_seed_other = 0,
    n_coupons = 3,
    target_size_rds = 50,
    
    # Prop sampling parameters
    target_size_prop = 400
  )


## STUDY 2
study_2 <- 
  list(
    N = 2000,
    K = 2,
    prev_K = c(known = .2, hidden = .1),
    rho_K = .05,
    p_edge_within = list(known = c(0.05, 0.05), hidden = c(0.05, 0.5)),
    p_edge_between = list(known = 0.05, hidden = 0.01),
    p_visibility = list(hidden = 1, known = 1),
    p_service = 0.2,
    n_seed_hidden = 20,
    n_seed_other = 0,
    n_coupons = 3,
    target_size_rds = 120,
    target_size_prop = 800
  )

# pop_args <- 
#   list(study_1 = study_1, study_2 = study_2)

pop_args <- 
  list(study_1 = study_1)


population <- 
  declare_population(
    handler = get_pop_network, 
    pop_args = pop_args
  )

population()

hpop_rds <- 
  declare_sampling(
    handler = get_hpop_sample, 
    sampling_variable = "rds",
    drop_nonsampled = FALSE,
    sampling_args = pop_args,
    sampling_strategy = "sample_prop") 

set.seed(19872312)
# hpop_rds(population())
draw_data(population + hpop_rds)
```
