
# hiddenmeta ![R](https://github.com/gerasy1987/hiddenmeta/workflows/R/badge.svg) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://github.com/gerasy1987/hiddenmeta/blob/main/LICENSE)

## Overview

``` r
if (!require(pacman)) install.packages("pacman")

pacman::p_load(
  tidyverse, knitr, magrittr, # data management/plotting/printing
  DeclareDesign, # Design declaration
  RcppAlgos, # speed up combinatorics
  networkreporting, NSUM, # NSUM methods 
  sspse, RDS, chords, # RDS+ methods
  igraph, 
  # ggnetwork, 
  # rstan, # Stan
  # statnet, network, sna, # network packages 
  # rgl, htmltools # network visualization
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
    # total population size for one study
    N = 2000,
    # number of groups
    # (K-th group is hidden population we are interested in)
    K = 2,
    # probability of membership in each of the groups (prev_K[K] is the true prevalence)
    prev_K = c(known = .2, hidden = .1),
    # correlation matrix of group memberships
    rho_K = .05,
    
    # block edge probabilities depending on group memberships
    # 1 - list of in- and out-group probability of links for each group
    # 2 - probability of link between in- and out-group members
    p_edge_within = list(known = c(0.05, 0.05), hidden = c(0.05, 0.5)),
    p_edge_between = list(known = 0.05, hidden = 0.01),
    
    # probability of visibility (show-up) for each group
    p_visibility = list(hidden = 1, known = 1),
    
    # probability of service utilization in hidden population
    # for service multiplier
    p_service = 0.2,
    
    # rds parameters
    n_seed_hidden = 20,
    n_seed_other = 0,
    n_coupons = 3,
    target_size_rds = 120,
    
    # Prop sampling parameters
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

hpop_rds <- 
  declare_sampling(
    handler = get_hpop_sample, 
    sampling_variable = "rds",
    drop_nonsampled = FALSE,
    sampling_args = pop_args,
    sampling_strategy = "sample_rds") 

pop_prop <- 
  declare_sampling(
    handler = get_hpop_sample, 
    sampling_variable = "prop",
    drop_nonsampled = FALSE,
    sampling_args = pop_args,
    sampling_strategy = "sample_prop") 

set.seed(19872312)
# hpop_rds(population())
# pop_prop(population())

draw_data(population + hpop_rds + pop_prop)

study_estimands <- 
  declare_estimand(handler = get_hpop_estimands)

# get_hpop_estimands(draw_data(population + hpop_rds + pop_prop))

draw_estimands(population + hpop_rds + pop_prop + study_estimands)
```
