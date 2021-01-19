
# hiddenmeta ![R](https://github.com/gerasy1987/hiddenmeta/workflows/R/badge.svg) [![codecov](https://codecov.io/gh/gerasy1987/hiddenmeta/branch/main/graph/badge.svg?token=ZG9A64Q0A1)](https://codecov.io/gh/gerasy1987/hiddenmeta) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://github.com/gerasy1987/hiddenmeta/blob/main/LICENSE)

## Install and load packages

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

## Step 1. Provide study design features

``` r
## STUDY 1
study_1 <- 
  list(
    network_handler = sim_block_network,
    network_handler_args = 
      list(N = 2000, K = 2, prev_K = c(known = .3, hidden = .1), rho_K = .05,
           p_edge_within = list(known = c(0.05, 0.05), hidden = c(0.05, 0.9)),
           p_edge_between = list(known = 0.05, hidden = 0.01),
           directed = FALSE),
    
    group_names = c("known", "hidden"),
    
    # probability of visibility (show-up) for each group
    p_visibility = list(known = .99, hidden = .7),
    
    # probability of service utilization in hidden population
    # for service multiplier
    add_groups = list(p_service = 0.3, 
                      loc_1 = 0.3, loc_2 = 0.1, loc_3 = 0.2, 
                      known_2 = 0.1, known_2 = 0.2),
    
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

## Step 2. Declare study population

``` r
study_population <-
  do.call(what = declare_population,
          args = c(handler = get_study_population, study_1[1:5]))

set.seed(19872312)
( example_pop <- study_population() )
```

### Show the network

``` r
g <-
  example_pop %$% {
    igraph::graph_from_adj_list(links,
                                mode = "all") %>%
      igraph::set_vertex_attr("name", value = name) %>%
      igraph::set_vertex_attr("type", value = type)
  }

igraph::V(g)$color <-
  plyr::mapvalues(igraph::V(g)$type,
                  from = unique(igraph::V(g)$type),
                  to = grDevices::palette.colors(n = length(unique(igraph::V(g)$type)), palette = "Set 3"))

plot(g,
     layout = igraph::layout_on_grid(g, dim = 2, width = 100),
     vertex.size = 3, vertex.dist = 2, vertex.label = NA, edge.width = 1,
     edge.arrow.size = .2, edge.curved = .2)

legend(x = -1, y = -1.2,
       legend = c("none", "known only", "hidden only", "both"),
       pt.bg = grDevices::palette.colors(n = length(unique(igraph::V(g)$type)), palette = "Set 3"),
       pch = 21, col = "#777777", pt.cex = 2, cex = 1.5, bty = "o", ncol = 2)
```

## Step 3. Declare all relevant study sampling procedures

The sampling procedures are additive in a sense that each procedure
appends several columns relevant to the sampling procedure and
particular draw based on population simulation, but does not change the
study population data frame (unless you specify
`drop_nonsampled = TRUE`).

``` r
study_sample_rds <- 
  do.call(declare_sampling,
          c(handler = sample_rds, 
            sampling_variable = "rds",
            drop_nonsampled = FALSE, study_1[6:9]))

set.seed(19872312)
draw_data(study_population + 
            study_sample_rds)
```

``` r
study_sample_pps <- 
  do.call(declare_sampling,
          c(handler = sample_pps, 
            sampling_variable = "pps",
            drop_nonsampled = FALSE, study_1[10]))

set.seed(19872312)
draw_data(study_population + 
            study_sample_rds + study_sample_pps)
```

``` r
study_sample_tls <- 
  do.call(declare_sampling,
          c(handler = sample_tls, 
            sampling_variable = "tls",
            drop_nonsampled = FALSE, study_1[11]))

set.seed(19872312)
draw_data(study_population + 
            study_sample_rds + study_sample_pps + study_sample_tls)
```

## Step 4. Declare study level estimands

``` r
study_estimands <- 
  declare_estimand(handler = get_study_estimands)

set.seed(19872312)
draw_estimands(study_population + 
                 study_sample_rds + study_sample_pps + study_sample_tls + 
                 study_estimands)
```

## Step 5. Declare estimators used in the study

``` r
estimator_sspse <- declare_estimator(handler = get_study_est_sspse, label = "sspse")
estimator_ht <- declare_estimator(handler = get_study_est_ht, label = "ht")
estimator_chords <- declare_estimator(type = "integrated",
                                      handler = get_study_est_chords, label = "chords")
estimator_nsum <- declare_estimator(handler = get_study_est_nsum, label = "nsum")

set.seed(19872312)
draw_estimates(study_population +
                 study_sample_rds + study_sample_pps + study_sample_tls +
                 study_estimands +
                 estimator_sspse + estimator_ht + estimator_chords + estimator_nsum)
```

## Step 6. Diagnose study design

``` r
study_diagnosands <-
  declare_diagnosands(
    mean_estimand = mean(estimand),
    mean_estimate = mean(estimate),
    sd_estimate = sd(estimate),
    mean_se = mean(se),
    bias = mean(estimate - estimand),
    rmse = sqrt(mean((estimate - estimand) ^ 2))
  )

diagnose_design(
  study_population + 
    study_sample_rds + study_sample_pps + 
    study_estimands + 
    estimator_sspse + estimator_ht + estimator_chords + estimator_nsum, 
  diagnosands = study_diagnosands,
  sims = 10,
  bootstrap_sims = 10)
```
