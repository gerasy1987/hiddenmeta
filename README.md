
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
    pop = 
      list(
        handler = get_study_population,
        
        # network structure setup
        network_handler = sim_block_network,
        network_handler_args = 
          list(N = 2000, K = 2, prev_K = c(known = .3, hidden = .1), rho_K = .05,
               p_edge_within = list(known = c(0.05, 0.05), hidden = c(0.05, 0.7)),
               p_edge_between = list(known = 0.05, hidden = 0.01),
               directed = FALSE),
        
        # groups
        group_names = c("known", "hidden"),
        
        # probability of visibility (show-up) for each group
        p_visible = list(known = 1, hidden = .7),
        
        # probability of service utilization in hidden population
        # for service multiplier
        add_groups = 
          list(service_use = "rbinom(n(), 1, 0.25 + hidden * 0.05)",
         "purrr::map_df(hidden, ~ sapply( `names<-`(c(.2, .1, .3, .2), paste0('loc_', 1:4)), function(add) rbinom(length(.x), 1, 0.1 + .x * add)))",
         known_2 = 0.1, known_3 = 0.2)
      ),
    sample = 
      list(
        rds = list(handler = sample_rds,
                   # RDS parameters
                   sampling_variable = "rds",
                   n_seed = 20,
                   n_coupons = 3,
                   target_type = "sample",
                   target_n_rds = 100),
        tls = list(handler = sample_tls,
                   sampling_variable = "tls",
                   # TLS sampling parameters
                   target_n_clusters = 3,
                   target_n_tls = 100,
                   cluster = paste0("loc_", 1:4)),
        pps = list(handler = sample_pps,
                   sampling_variable = "pps",
                   # prop sampling parameters
                   sampling_frame = NULL,
                   strata = NULL,
                   cluster = NULL,
                   target_n_pps = 400)
      ),
    inquiries = list(handler = get_study_estimands,
                     known_pattern = "^known", 
                     hidden_pattern = "^hidden$"),
    estimators = 
      list(
        rds = 
          list(sspse = list(handler = get_study_est_sspse,
                            prior_mean = 200,
                            mcmc_params = list(interval = 5, burnin = 2000, samplesize = 500),
                            rds_prefix = "rds", 
                            label = "rds_sspse"),
               chords = list(handler = get_study_est_chords, 
                             type = "mle",
                             seed_condition = "rds_from == -999",
                             n_boot = 100,
                             rds_prefix = "rds",
                             label = "chords"),
               multiplier = list(handler = get_study_est_multiplier, 
                                 service_var = "service_use",
                                 seed_condition = "rds_from == -999",
                                 n_boot = 100,
                                 rds_prefix = "rds",
                                 label = "rds_multi")),
        tls =
          list(ht = list(handler = get_study_est_ht,
                         prefix = "tls",
                         label = "tls_ht"),
               nsum = list(handler = get_study_est_nsum,
                           known = c("known", "known_2", "known_3"),
                           hidden = "hidden_visible_out",
                           survey_design = ~ tls_cluster,
                           n_boot = 100,
                           prefix = "tls",
                           label = "tls_nsum")),
        pps = 
          list(ht = list(handler = get_study_est_ht,
                         prefix = "pps",
                         label = "pps_ht"),
               nsum = list(handler = get_study_est_nsum,
                           known = c("known", "known_2", "known_3"),
                           hidden = "hidden_visible_out",
                           survey_design = ~ pps_cluster + strata(pps_strata),
                           n_boot = 100,
                           prefix = "pps",
                           label = "pps_nsum")),
        all = 
          list(recap = list(handler = get_study_est_recapture,
                            capture_vars = paste0("loc_", 1:4),
                            model = "Mt",
                            hidden_condition = "hidden == 1",
                            label = "all_recap"),
               recap2 = list(handler = get_study_est_recapture,
                             capture_vars = c("rds", "pps"),
                             model = "Mt",
                             hidden_condition = "hidden == 1",
                             label = "samp_recap"))
      )
  )
```

## Step 2. Declare study population

``` r
study_population <-
  do.call(what = declare_population,
          args = study_1$pop)

set.seed(19872312)
( example_pop <- study_population() )
#> # A tibble: 2,000 x 51
#>     name type  known hidden links    service_use loc_1 loc_2 loc_3 loc_4 known_2
#>    <int> <chr> <int>  <int> <chr>          <int> <int> <int> <int> <int>   <int>
#>  1     1 00        0      0 166;107…           0     0     0     0     0       0
#>  2     2 00        0      0 84;199;…           0     0     0     0     0       0
#>  3     3 00        0      0 143;167…           1     1     0     0     0       0
#>  4     4 00        0      0 272;300…           0     1     0     0     1       0
#>  5     5 00        0      0 29;536;…           0     0     1     0     0       0
#>  6     6 00        0      0 745;103…           1     0     1     0     0       0
#>  7     7 00        0      0 66;933;…           0     0     0     0     0       0
#>  8     8 00        0      0 526;107…           1     0     0     0     0       0
#>  9     9 00        0      0 111;1078           0     0     0     0     0       0
#> 10    10 00        0      0 377;1992           0     0     0     0     1       0
#> # … with 1,990 more rows, and 40 more variables: known_3 <int>,
#> #   n_visible_out <dbl>, known_visible_out <dbl>, hidden_visible_out <dbl>,
#> #   type_00_visible_out <dbl>, type_01_visible_out <dbl>,
#> #   type_10_visible_out <dbl>, type_11_visible_out <dbl>,
#> #   service_use_visible_out <dbl>, loc_1_visible_out <dbl>,
#> #   loc_2_visible_out <dbl>, loc_3_visible_out <dbl>, loc_4_visible_out <dbl>,
#> #   known_2_visible_out <dbl>, known_3_visible_out <dbl>,
#> #   known_visible_in <dbl>, hidden_visible_in <dbl>, type_00_visible_in <dbl>,
#> #   type_01_visible_in <dbl>, type_10_visible_in <dbl>,
#> #   type_11_visible_in <dbl>, service_use_visible_in <dbl>,
#> #   loc_1_visible_in <dbl>, loc_2_visible_in <dbl>, loc_3_visible_in <dbl>,
#> #   loc_4_visible_in <dbl>, known_2_visible_in <dbl>, known_3_visible_in <dbl>,
#> #   p_visible_known <dbl>, p_visible_hidden <dbl>, total <int>,
#> #   total_known <int>, total_hidden <int>, total_service_use <int>,
#> #   total_loc_1 <int>, total_loc_2 <int>, total_loc_3 <int>, total_loc_4 <int>,
#> #   total_known_2 <int>, total_known_3 <int>
```

### Show the network

``` r
g <-
  example_pop %$% {
    hiddenmeta:::retrieve_graph(links) %>%
      igraph::set_vertex_attr("name", value = name) %>%
      igraph::set_vertex_attr("type", value = type)
  }

igraph::V(g)$color <-
  plyr::mapvalues(igraph::V(g)$type,
                  from = unique(igraph::V(g)$type),
                  to = grDevices::palette.colors(n = length(unique(igraph::V(g)$type)), 
                                                 palette = "Set 3"))

plot(g,
     layout = igraph::layout_on_grid(g, dim = 2, width = 100),
     vertex.size = 3, vertex.dist = 2, vertex.label = NA, edge.width = 1,
     edge.arrow.size = .2, edge.curved = .2)

legend(x = -1, y = -1.2,
       legend = c("none", "known only", "hidden only", "both"),
       pt.bg = grDevices::palette.colors(n = length(unique(igraph::V(g)$type)), palette = "Set 3"),
       pch = 21, col = "#777777", pt.cex = 2, cex = 1.5, bty = "o", ncol = 2)
```

![](README-plot_pop_network-1.png)<!-- -->

## Step 3. Declare all relevant study sampling procedures

The sampling procedures are additive in a sense that each procedure
appends several columns relevant to the sampling procedure and
particular draw based on population simulation, but does not change the
study population data frame (unless you specify
`drop_nonsampled = TRUE`).

``` r
study_sample_rds <- 
  do.call(what = declare_sampling, 
          args = study_1$sample$rds)

set.seed(19872312)
draw_data(study_population + study_sample_rds)
#> # A tibble: 2,000 x 60
#>     name type  known hidden links    service_use loc_1 loc_2 loc_3 loc_4 known_2
#>    <int> <chr> <int>  <int> <chr>          <int> <int> <int> <int> <int>   <int>
#>  1     1 00        0      0 166;107…           0     0     0     0     0       0
#>  2     2 00        0      0 84;199;…           0     0     0     0     0       0
#>  3     3 00        0      0 143;167…           1     1     0     0     0       0
#>  4     4 00        0      0 272;300…           0     1     0     0     1       0
#>  5     5 00        0      0 29;536;…           0     0     1     0     0       0
#>  6     6 00        0      0 745;103…           1     0     1     0     0       0
#>  7     7 00        0      0 66;933;…           0     0     0     0     0       0
#>  8     8 00        0      0 526;107…           1     0     0     0     0       0
#>  9     9 00        0      0 111;1078           0     0     0     0     0       0
#> 10    10 00        0      0 377;1992           0     0     0     0     1       0
#> # … with 1,990 more rows, and 49 more variables: known_3 <int>,
#> #   n_visible_out <dbl>, known_visible_out <dbl>, hidden_visible_out <dbl>,
#> #   type_00_visible_out <dbl>, type_01_visible_out <dbl>,
#> #   type_10_visible_out <dbl>, type_11_visible_out <dbl>,
#> #   service_use_visible_out <dbl>, loc_1_visible_out <dbl>,
#> #   loc_2_visible_out <dbl>, loc_3_visible_out <dbl>, loc_4_visible_out <dbl>,
#> #   known_2_visible_out <dbl>, known_3_visible_out <dbl>,
#> #   known_visible_in <dbl>, hidden_visible_in <dbl>, type_00_visible_in <dbl>,
#> #   type_01_visible_in <dbl>, type_10_visible_in <dbl>,
#> #   type_11_visible_in <dbl>, service_use_visible_in <dbl>,
#> #   loc_1_visible_in <dbl>, loc_2_visible_in <dbl>, loc_3_visible_in <dbl>,
#> #   loc_4_visible_in <dbl>, known_2_visible_in <dbl>, known_3_visible_in <dbl>,
#> #   p_visible_known <dbl>, p_visible_hidden <dbl>, total <int>,
#> #   total_known <int>, total_hidden <int>, total_service_use <int>,
#> #   total_loc_1 <int>, total_loc_2 <int>, total_loc_3 <int>, total_loc_4 <int>,
#> #   total_known_2 <int>, total_known_3 <int>, rds <int>, rds_from <dbl>,
#> #   rds_t <dbl>, rds_wave <dbl>, rds_hidden <int>, rds_own_coupon <chr>,
#> #   rds_coupon_1 <chr>, rds_coupon_2 <chr>, rds_coupon_3 <chr>
```

``` r
study_sample_pps <- 
  do.call(what = declare_sampling, 
          args = study_1$sample$pps)

set.seed(19872312)
draw_data(study_population + study_sample_rds + study_sample_pps)
#> # A tibble: 2,000 x 68
#>     name type  known hidden links    service_use loc_1 loc_2 loc_3 loc_4 known_2
#>    <int> <chr> <int>  <int> <chr>          <int> <int> <int> <int> <int>   <int>
#>  1     1 00        0      0 166;107…           0     0     0     0     0       0
#>  2     2 00        0      0 84;199;…           0     0     0     0     0       0
#>  3     3 00        0      0 143;167…           1     1     0     0     0       0
#>  4     4 00        0      0 272;300…           0     1     0     0     1       0
#>  5     5 00        0      0 29;536;…           0     0     1     0     0       0
#>  6     6 00        0      0 745;103…           1     0     1     0     0       0
#>  7     7 00        0      0 66;933;…           0     0     0     0     0       0
#>  8     8 00        0      0 526;107…           1     0     0     0     0       0
#>  9     9 00        0      0 111;1078           0     0     0     0     0       0
#> 10    10 00        0      0 377;1992           0     0     0     0     1       0
#> # … with 1,990 more rows, and 57 more variables: known_3 <int>,
#> #   n_visible_out <dbl>, known_visible_out <dbl>, hidden_visible_out <dbl>,
#> #   type_00_visible_out <dbl>, type_01_visible_out <dbl>,
#> #   type_10_visible_out <dbl>, type_11_visible_out <dbl>,
#> #   service_use_visible_out <dbl>, loc_1_visible_out <dbl>,
#> #   loc_2_visible_out <dbl>, loc_3_visible_out <dbl>, loc_4_visible_out <dbl>,
#> #   known_2_visible_out <dbl>, known_3_visible_out <dbl>,
#> #   known_visible_in <dbl>, hidden_visible_in <dbl>, type_00_visible_in <dbl>,
#> #   type_01_visible_in <dbl>, type_10_visible_in <dbl>,
#> #   type_11_visible_in <dbl>, service_use_visible_in <dbl>,
#> #   loc_1_visible_in <dbl>, loc_2_visible_in <dbl>, loc_3_visible_in <dbl>,
#> #   loc_4_visible_in <dbl>, known_2_visible_in <dbl>, known_3_visible_in <dbl>,
#> #   p_visible_known <dbl>, p_visible_hidden <dbl>, total <int>,
#> #   total_known <int>, total_hidden <int>, total_service_use <int>,
#> #   total_loc_1 <int>, total_loc_2 <int>, total_loc_3 <int>, total_loc_4 <int>,
#> #   total_known_2 <int>, total_known_3 <int>, rds <int>, rds_from <dbl>,
#> #   rds_t <dbl>, rds_wave <dbl>, rds_hidden <int>, rds_own_coupon <chr>,
#> #   rds_coupon_1 <chr>, rds_coupon_2 <chr>, rds_coupon_3 <chr>,
#> #   pps_frame <dbl>, pps_strata <dbl>, pps_strata_prop <dbl>,
#> #   pps_cluster <int>, pps_cluster_prop <dbl>, pps <dbl>,
#> #   pps_sampled_cluster <dbl>, pps_weight <dbl>
```

``` r
study_sample_tls <- 
  do.call(what = declare_sampling, 
          args = study_1$sample$tls)

set.seed(19872312)
draw_data(study_population + 
            study_sample_rds + study_sample_pps + study_sample_tls)
#> # A tibble: 2,000 x 73
#>     name type  known hidden links    service_use loc_1 loc_2 loc_3 loc_4 known_2
#>    <int> <chr> <int>  <int> <chr>          <int> <int> <int> <int> <int>   <int>
#>  1     1 00        0      0 166;107…           0     0     0     0     0       0
#>  2     2 00        0      0 84;199;…           0     0     0     0     0       0
#>  3     3 00        0      0 143;167…           1     1     0     0     0       0
#>  4     4 00        0      0 272;300…           0     1     0     0     1       0
#>  5     5 00        0      0 29;536;…           0     0     1     0     0       0
#>  6     6 00        0      0 745;103…           1     0     1     0     0       0
#>  7     7 00        0      0 66;933;…           0     0     0     0     0       0
#>  8     8 00        0      0 526;107…           1     0     0     0     0       0
#>  9     9 00        0      0 111;1078           0     0     0     0     0       0
#> 10    10 00        0      0 377;1992           0     0     0     0     1       0
#> # … with 1,990 more rows, and 62 more variables: known_3 <int>,
#> #   n_visible_out <dbl>, known_visible_out <dbl>, hidden_visible_out <dbl>,
#> #   type_00_visible_out <dbl>, type_01_visible_out <dbl>,
#> #   type_10_visible_out <dbl>, type_11_visible_out <dbl>,
#> #   service_use_visible_out <dbl>, loc_1_visible_out <dbl>,
#> #   loc_2_visible_out <dbl>, loc_3_visible_out <dbl>, loc_4_visible_out <dbl>,
#> #   known_2_visible_out <dbl>, known_3_visible_out <dbl>,
#> #   known_visible_in <dbl>, hidden_visible_in <dbl>, type_00_visible_in <dbl>,
#> #   type_01_visible_in <dbl>, type_10_visible_in <dbl>,
#> #   type_11_visible_in <dbl>, service_use_visible_in <dbl>,
#> #   loc_1_visible_in <dbl>, loc_2_visible_in <dbl>, loc_3_visible_in <dbl>,
#> #   loc_4_visible_in <dbl>, known_2_visible_in <dbl>, known_3_visible_in <dbl>,
#> #   p_visible_known <dbl>, p_visible_hidden <dbl>, total <int>,
#> #   total_known <int>, total_hidden <int>, total_service_use <int>,
#> #   total_loc_1 <int>, total_loc_2 <int>, total_loc_3 <int>, total_loc_4 <int>,
#> #   total_known_2 <int>, total_known_3 <int>, rds <int>, rds_from <dbl>,
#> #   rds_t <dbl>, rds_wave <dbl>, rds_hidden <int>, rds_own_coupon <chr>,
#> #   rds_coupon_1 <chr>, rds_coupon_2 <chr>, rds_coupon_3 <chr>,
#> #   pps_frame <dbl>, pps_strata <dbl>, pps_strata_prop <dbl>,
#> #   pps_cluster <int>, pps_cluster_prop <dbl>, pps <dbl>,
#> #   pps_sampled_cluster <dbl>, pps_weight <dbl>, tls <dbl>,
#> #   tls_loc_present <chr>, tls_locs_sampled <chr>, tls_cluster <chr>,
#> #   tls_weight <dbl>
```

## Step 4. Declare study level estimands

``` r
study_estimands <- 
  do.call(what = declare_inquiry, 
          args = study_1$inquiries)

set.seed(19872312)
draw_inquiry(study_population + 
                # study_sample_rds + study_sample_pps + study_sample_tls + 
                study_estimands)
#>           inquiry_label estimand
#> 1           hidden_size 190.0000
#> 2           hidden_prev   0.0950
#> 3        degree_average   4.9310
#> 4 degree_hidden_average   0.7155
```

## Step 5. Declare estimators used in the study

``` r
est_sspse <- do.call(what = declare_estimator, 
                     args = study_1$estimators$rds$sspse)
est_chords <- do.call(what = declare_estimator, 
                      args = study_1$estimators$rds$chords)
est_multi <- do.call(what = declare_estimator, 
                     args = study_1$estimators$rds$multiplier)

est_ht_pps <- do.call(what = declare_estimator, 
                      args = study_1$estimators$pps$ht)
est_ht_tls <- do.call(what = declare_estimator, 
                      args = study_1$estimators$tls$ht)

est_nsum_tls <- do.call(what = declare_estimator,
                        args = study_1$estimators$tls$nsum)
est_nsum_pps <- do.call(what = declare_estimator, 
                        args = study_1$estimators$pps$nsum)

est_recap <- do.call(what = declare_estimator, 
                     args = study_1$estimators$all$recap)

est_recap2 <- do.call(what = declare_estimator, 
                      args = study_1$estimators$all$recap2)

set.seed(19872312)
draw_estimates(study_population +
                 study_sample_rds + study_sample_pps + study_sample_tls +
                 study_estimands +
                 est_sspse + est_chords + est_multi +
                 est_ht_tls + est_ht_pps +
                 est_nsum_tls + est_nsum_pps +
                 est_recap + est_recap2)
#>           estimator_label    estimate           se         inquiry_label
#> 1   hidden_size_rds_sspse 141.0000000 39.602863401           hidden_size
#> 2      hidden_size_chords 158.0000000 59.380276799           hidden_size
#> 3    degree_hidden_chords   4.9240506  0.807132243 degree_hidden_average
#> 4   hidden_size_rds_multi 180.5555556 25.649056775           hidden_size
#> 5      hidden_size_tls_ht  80.0000000 17.059572100           hidden_size
#> 6      hidden_prev_tls_ht   0.0400000  0.008529786           hidden_prev
#> 7      hidden_size_pps_ht 155.0000000 27.265348311           hidden_size
#> 8      hidden_prev_pps_ht   0.0775000  0.013632674           hidden_prev
#> 9    hidden_size_tls_nsum 433.3825108 55.589475936           hidden_size
#> 10        degree_tls_nsum   0.4857944  0.132946002        degree_average
#> 11   hidden_size_pps_nsum 170.1373183 28.171363409           hidden_size
#> 12        degree_pps_nsum  58.7286145 10.473314912        degree_average
#> 13  hidden_size_all_recap 198.1868365 11.900941786           hidden_size
#> 14 hidden_size_samp_recap 201.0000000 31.780497164           hidden_size
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

study_simulations <- 
  simulate_design(
    study_population +
      study_sample_rds + study_sample_pps + study_sample_tls +
      study_estimands +
      est_sspse + est_chords + est_multi +
      est_ht_tls + est_ht_pps +
      est_nsum_tls + est_nsum_pps +
      est_recap + est_recap2, 
    sims = 10)

diagnose_design(study_simulations, 
                diagnosands = study_diagnosands, sims = 2)
#> 
#> Research design diagnosis based on 10 simulations. Diagnosand estimates with bootstrapped standard errors in parentheses (100 replicates).
#> 
#>  Design Label         Inquiry Label        Estimator Label N Sims Mean Estimand
#>      design_1        degree_average        degree_pps_nsum     10          4.94
#>                                                                          (0.03)
#>      design_1        degree_average        degree_tls_nsum     10          4.94
#>                                                                          (0.03)
#>      design_1 degree_hidden_average   degree_hidden_chords     10          0.82
#>                                                                          (0.03)
#>      design_1           hidden_prev     hidden_prev_pps_ht     10          0.10
#>                                                                          (0.00)
#>      design_1           hidden_prev     hidden_prev_tls_ht     10          0.10
#>                                                                          (0.00)
#>      design_1           hidden_size  hidden_size_all_recap     10        203.50
#>                                                                          (3.32)
#>      design_1           hidden_size     hidden_size_chords     10        203.50
#>                                                                          (3.32)
#>      design_1           hidden_size     hidden_size_pps_ht     10        203.50
#>                                                                          (3.32)
#>      design_1           hidden_size   hidden_size_pps_nsum     10        203.50
#>                                                                          (3.32)
#>      design_1           hidden_size  hidden_size_rds_multi     10        203.50
#>                                                                          (3.32)
#>      design_1           hidden_size  hidden_size_rds_sspse     10        203.50
#>                                                                          (3.32)
#>      design_1           hidden_size hidden_size_samp_recap     10        203.50
#>                                                                          (3.32)
#>      design_1           hidden_size     hidden_size_tls_ht     10        203.50
#>                                                                          (3.32)
#>      design_1           hidden_size   hidden_size_tls_nsum     10        203.50
#>                                                                          (3.32)
#>  Mean Estimate SD Estimate Mean Se    Bias    RMSE
#>          44.70        7.49    6.21   39.77   40.41
#>         (2.32)      (1.06)  (0.50)  (2.34)  (2.31)
#>           0.46        0.07    0.11   -4.47    4.48
#>         (0.02)      (0.02)  (0.02)  (0.04)  (0.04)
#>           5.11        0.65    0.76    4.29    4.32
#>         (0.21)      (0.17)  (0.02)  (0.19)  (0.20)
#>           0.10        0.02    0.02    0.00    0.01
#>         (0.00)      (0.00)  (0.00)  (0.00)  (0.00)
#>           0.06        0.01    0.01   -0.04    0.04
#>         (0.00)      (0.00)  (0.00)  (0.00)  (0.01)
#>         199.21       16.21   11.48   -4.29   15.96
#>         (4.44)      (2.79)  (0.77)  (4.70)  (2.22)
#>         181.20       68.54   58.22  -22.30   67.80
#>        (18.86)     (16.49)  (2.27) (18.85) (18.05)
#>         208.00       30.48   31.04    4.50   23.73
#>         (9.49)      (4.40)  (0.80)  (7.59)  (2.72)
#>         234.51       45.74   30.81   31.01   47.22
#>        (14.31)      (9.23)  (1.18) (11.72) (11.33)
#>         204.39       15.86   36.45    0.89   11.04
#>         (5.06)      (2.89)  (1.53)  (3.59)  (2.54)
#>         161.85       20.94   64.25  -41.65   49.75
#>         (6.91)      (3.05)  (5.12)  (9.17)  (7.72)
#>         190.16       17.44   24.18  -13.34   20.67
#>         (4.43)      (2.79)  (1.42)  (4.50)  (4.28)
#>         124.00       23.07   20.96  -79.50   84.84
#>         (7.19)      (5.68)  (0.39)  (9.73) (10.18)
#>         562.51       84.41   63.66  359.01  367.78
#>        (23.19)     (22.14) (11.44) (22.75) (19.97)
```

## Meta-analysis Draft structure

1.  Get all estimates produced for each of the studies (depending on the
    design) in a single data frame that will include estimator and
    standard error for each estimator (or potentially for each
    sampling-estimator pair)

2.  Then we can *feed* this data frame into Stan model akin the one
    below

    1.  Use the data frame to produce observed estimator vectors for
        each estimator/sampling-estimator pair, `observed*`
    2.  Use the data frame produced by study designs to produce arrays
        of observed estimates and corresponding standard errors (with
        dummy value where missing), `est*` and `est*_sd`
    3.  Assume that deviation from true parameter of interest (`error`)
        is driven by use of estimator

3.  We then declare the relevant estimands and run diagnosis

``` r
stan_model_meta <- " 
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
  "

get_meta_estimands <- function(data) {
  
  data.frame(inquiry_label = c(paste0("prevalence_", 1:N)),
             estimand = c(data[,1]),
             stringsAsFactors = FALSE)
}

get_meta_estimators = function(data) {
  
  stan_data <- list(N = nrow(data),
                    K = (ncol(data)-1)/2)
  
  for (k in 1:stan_data$K) {
    stan_data[[ paste0("observed",k) ]] <- 
      which(!is.na(data[,(2 * k)]))
    stan_data[[ paste0("N",k) ]] <- 
      length(stan_data[[paste0("observed",k)]])
    stan_data[[ paste0("est",k) ]] <- 
      data[stan_data[[paste0("observed",k)]],(2 * k)]
    stan_data[[ paste0("est",k,"_sd") ]] <- 
      data[stan_data[[paste0("observed",k)]],(1 + 2 * k)]
  }
  
  fit <- 
    rstan::stan(fit = stan_model_meta, 
                data = stan_data, 
                iter = 4000) %>% 
    extract
  
  data.frame(estimator_label = c(paste0("prev_", 1:N)),
             estimate = c(apply(fit$alpha, 2, mean)),
             sd =   c(apply(fit$alpha, 2, sd)),
             inquiry_label = c(paste0("hidden_prev", 1:N)),
             big_Rhat = big_Rhat
             )
  
  }
```

## Meta declaration

### Additional study

``` r
study_2 <- 
  list(
    pop = 
      list(
        handler = get_study_population,
        network_handler = sim_block_network,
        network_handler_args = 
          list(N = 5000, K = 3, 
               prev_K = c(frame = .5, known = .2, hidden = .1), 
               rho_K = c(.05, .05, .05),
               p_edge_within = list(frame = c(0.05, 0.05), 
                                    known = c(0.1, 0.05), 
                                    hidden = c(0.2, 0.9)),
               p_edge_between = list(frame = 0.05, 
                                     known = 0.1, 
                                     hidden = 0.01),
               directed = FALSE),
        
        group_names = c("frame", "known", "hidden"),
        
        # probability of visibility (show-up) for each group
        p_visible = list(frame = 1, known = 1, hidden = .6),
        
        # probability of service utilization in hidden population
        # for service multiplier
        add_groups = list(p_service = 0.3, 
                          loc_1 = 0.3, loc_2 = 0.1, loc_3 = 0.2)
      ),
    sample = 
      list(
        rds = list(handler = sample_rds,
                   # RDS parameters
                   sampling_variable = "rds",
                   n_seed = 100,
                   n_coupons = 3,
                   target_type = "waves",
                   target_n_rds = 2),
        pps = list(handler = sample_pps,
                   sampling_variable = "pps",
                   # prop sampling parameters
                   sampling_frame = "frame",
                   strata = NULL,
                   cluster = NULL,
                   target_n_pps = 800)
      ),
    inquiries = list(handler = get_study_estimands,
                     known_pattern = "^known", 
                     hidden_pattern = "^hidden$"),
    estimators = 
      list(
        rds = 
          list(sspse = list(handler = get_study_est_sspse, 
                            label = "rds_sspse",
                            prior_mean = 450,
                            mcmc_params = list(interval = 5, burnin = 2000, samplesize = 500),
                            rds_prefix = "rds"),
               chords = list(handler = get_study_est_chords, 
                             type = "mle",
                             label = "rds_chords",
                             seed_condition = "rds_from == -999",
                             n_boot = 100,
                             rds_prefix = "rds")),
        pps = 
          list(
            ht = list(handler = get_study_est_ht, 
                      label = "pps_ht"),
            nsum = list(handler = get_study_est_nsum,
                        known = c("frame", "known"),
                        hidden = "hidden_visible_out",
                        survey_design = ~ pps_cluster + strata(pps_strata),
                        n_boot = 100,
                        label = "pps_nsum")),
        all = 
          list(recap = list(handler = get_study_est_recapture,
                            capture_vars = paste0("loc_", 1:2),
                            model = "Mt",
                            hidden_condition = "hidden == 1",
                            label = "all_recap"),
               recap2 = list(handler = get_study_est_recapture,
                             capture_vars = c("rds", "pps"),
                             model = "Mt",
                             hidden_condition = "hidden == 1",
                             label = "samp_recap"))
      )
  )
```

### Declare multiple studies

``` r
multi_population <- 
  declare_population(handler = get_multi_populations, 
                     pops_args = list(study_1 = study_1$pop,
                                      study_2 = study_2$pop))

multi_sampling <- 
  declare_sampling(handler = get_multi_samples, 
                   samples_args = list(study_1 = study_1$sample,
                                       study_2 = study_2$sample)) 

multi_inquiry <- 
  declare_inquiry(handler = get_multi_estimands, 
                  inquiries_args = list(study_1 = study_1$inquiries,
                                        study_2 = study_2$inquiries)) 

multi_estimators <-
  declare_estimator(handler = get_multi_estimates,
                    estimators_args = list(study_1 = study_1$estimators,
                                           study_2 = study_2$estimators))

multi_design <- multi_population + multi_sampling + multi_inquiry + multi_estimators


set.seed(19872312)
draw_data(multi_design)
#> # A tibble: 2 x 2
#>   study   population                
#>   <chr>   <list>                    
#> 1 study_1 <tibble[,73] [2,000 × 73]>
#> 2 study_2 <tibble[,69] [5,000 × 69]>
draw_inquiry(multi_design)
#>                   inquiry_label estimand
#> 1           study_1_hidden_size 190.0000
#> 2           study_1_hidden_prev   0.0950
#> 3        study_1_degree_average   4.8170
#> 4 study_1_degree_hidden_average   0.7115
#> 5           study_2_hidden_size 488.0000
#> 6           study_2_hidden_prev   0.0976
#> 7        study_2_degree_average   4.3252
#> 8 study_2_degree_hidden_average   0.2392
draw_estimates(multi_design)
#>             estimator_label    estimate           se
#> 1     hidden_size_rds_sspse 183.0000000 7.171208e+01
#> 2        hidden_size_chords 147.0000000 4.814045e+01
#> 3      degree_hidden_chords   6.1632653 6.233803e-01
#> 4     hidden_size_rds_multi 194.2857143 2.818115e+01
#> 5        hidden_size_tls_ht 110.0000000 2.296959e+01
#> 6        hidden_prev_tls_ht   0.0550000 1.148480e-02
#> 7      hidden_size_tls_nsum 482.3374887 6.767643e+01
#> 8           degree_tls_nsum   0.5324749 1.831385e-01
#> 9        hidden_size_pps_ht 195.0000000 3.014279e+01
#> 10       hidden_prev_pps_ht   0.0975000 1.507139e-02
#> 11     hidden_size_pps_nsum 265.6483899 3.360678e+01
#> 12          degree_pps_nsum  37.1267484 4.620739e+00
#> 13    hidden_size_all_recap 208.1144228 1.310946e+01
#> 14   hidden_size_samp_recap 138.3103448 1.126146e+01
#> 15    hidden_size_rds_sspse 333.0000000 1.005618e+02
#> 16   hidden_size_rds_chords  67.0000000 1.133623e+02
#> 17 degree_hidden_rds_chords   3.0000000 9.051586e-01
#> 18       hidden_size_pps_ht 256.8650000 2.710808e+01
#> 19       hidden_prev_pps_ht   0.0513730 5.421615e-03
#> 20     hidden_size_pps_nsum 175.5210595 2.548890e+01
#> 21          degree_pps_nsum  58.4656067 9.342514e+00
#> 22    hidden_size_all_recap 458.0000000 7.842576e+01
#> 23   hidden_size_samp_recap 485.9333333 6.365135e+01
#>                    inquiry_label
#> 1            study_1_hidden_size
#> 2            study_1_hidden_size
#> 3  study_1_degree_hidden_average
#> 4            study_1_hidden_size
#> 5            study_1_hidden_size
#> 6            study_1_hidden_prev
#> 7            study_1_hidden_size
#> 8         study_1_degree_average
#> 9            study_1_hidden_size
#> 10           study_1_hidden_prev
#> 11           study_1_hidden_size
#> 12        study_1_degree_average
#> 13           study_1_hidden_size
#> 14           study_1_hidden_size
#> 15           study_2_hidden_size
#> 16           study_2_hidden_size
#> 17 study_2_degree_hidden_average
#> 18           study_2_hidden_size
#> 19           study_2_hidden_prev
#> 20           study_2_hidden_size
#> 21        study_2_degree_average
#> 22           study_2_hidden_size
#> 23           study_2_hidden_size

multi_simulations <- 
  simulate_design(multi_population +
                    multi_sampling + 
                    multi_inquiry + 
                    multi_estimators, 
                  sims = 10)

diagnose_design(simulations_df = multi_simulations, 
                diagnosands = study_diagnosands)
#> 
#> Research design diagnosis based on 10 simulations. Diagnosand estimates with bootstrapped standard errors in parentheses (100 replicates).
#> 
#>  Design Label                 Inquiry Label          Estimator Label N Sims
#>      design_1        study_1_degree_average          degree_pps_nsum     10
#>                                                                            
#>      design_1        study_1_degree_average          degree_tls_nsum     10
#>                                                                            
#>      design_1 study_1_degree_hidden_average     degree_hidden_chords     10
#>                                                                            
#>      design_1           study_1_hidden_prev       hidden_prev_pps_ht     10
#>                                                                            
#>      design_1           study_1_hidden_prev       hidden_prev_tls_ht     10
#>                                                                            
#>      design_1           study_1_hidden_size    hidden_size_all_recap     10
#>                                                                            
#>      design_1           study_1_hidden_size       hidden_size_chords     10
#>                                                                            
#>      design_1           study_1_hidden_size       hidden_size_pps_ht     10
#>                                                                            
#>      design_1           study_1_hidden_size     hidden_size_pps_nsum     10
#>                                                                            
#>      design_1           study_1_hidden_size    hidden_size_rds_multi     10
#>                                                                            
#>      design_1           study_1_hidden_size    hidden_size_rds_sspse     10
#>                                                                            
#>      design_1           study_1_hidden_size   hidden_size_samp_recap     10
#>                                                                            
#>      design_1           study_1_hidden_size       hidden_size_tls_ht     10
#>                                                                            
#>      design_1           study_1_hidden_size     hidden_size_tls_nsum     10
#>                                                                            
#>      design_1        study_2_degree_average          degree_pps_nsum     10
#>                                                                            
#>      design_1 study_2_degree_hidden_average degree_hidden_rds_chords     10
#>                                                                            
#>      design_1           study_2_hidden_prev       hidden_prev_pps_ht     10
#>                                                                            
#>      design_1           study_2_hidden_size    hidden_size_all_recap     10
#>                                                                            
#>      design_1           study_2_hidden_size       hidden_size_pps_ht     10
#>                                                                            
#>      design_1           study_2_hidden_size     hidden_size_pps_nsum     10
#>                                                                            
#>      design_1           study_2_hidden_size   hidden_size_rds_chords     10
#>                                                                            
#>      design_1           study_2_hidden_size    hidden_size_rds_sspse     10
#>                                                                            
#>      design_1           study_2_hidden_size   hidden_size_samp_recap     10
#>                                                                            
#>  Mean Estimand Mean Estimate SD Estimate Mean Se     Bias     RMSE
#>           4.93         47.08        8.34    6.70    42.15    42.89
#>         (0.03)        (2.55)      (1.86)  (0.53)   (2.55)   (2.63)
#>           4.93          0.46        0.10    0.12    -4.47     4.47
#>         (0.03)        (0.03)      (0.01)  (0.02)   (0.02)   (0.02)
#>           0.77          5.13        0.51    0.70     4.36     4.38
#>         (0.02)        (0.15)      (0.11)  (0.05)   (0.14)   (0.13)
#>           0.10          0.10        0.02    0.01     0.00     0.02
#>         (0.00)        (0.01)      (0.00)  (0.00)   (0.01)   (0.00)
#>           0.10          0.06        0.01    0.01    -0.04     0.04
#>         (0.00)        (0.00)      (0.00)  (0.00)   (0.00)   (0.00)
#>         196.30        200.40       14.82   12.79     4.10    15.43
#>         (3.34)        (4.80)      (2.55)  (0.73)   (4.89)   (3.71)
#>         196.30        154.00       68.69   59.86   -42.30    74.18
#>         (3.34)       (20.69)      (9.83)  (3.47)  (19.46)   (9.37)
#>         196.30        199.00       37.77   29.34     2.70    35.33
#>         (3.34)       (11.62)      (5.46)  (0.82)  (11.35)   (4.52)
#>         196.30        219.69       40.04   29.70    23.39    43.29
#>         (3.34)       (12.48)      (7.38)  (1.38)  (11.97)   (9.76)
#>         196.30        210.77       31.89   35.46    14.47    30.39
#>         (3.34)       (10.18)      (4.75)  (3.30)   (9.08)   (8.01)
#>         196.30        169.35       27.42   61.36   -26.95    40.96
#>         (3.34)        (7.89)      (6.56)  (5.54)   (9.63)   (6.81)
#>         196.30        189.56       22.36   24.87    -6.74    22.31
#>         (3.34)        (6.38)      (5.02)  (2.02)   (6.71)   (3.24)
#>         196.30        122.50       22.39   22.10   -73.80    75.29
#>         (3.34)        (6.40)      (3.28)  (0.55)   (4.52)   (4.29)
#>         196.30        546.89       93.24   82.69   350.59   360.09
#>         (3.34)       (25.67)     (20.00) (13.05)  (23.83)  (25.75)
#>           4.24         49.50       10.52    7.25    45.26    46.35
#>         (0.01)        (2.81)      (2.53)  (0.88)   (2.81)   (3.07)
#>           0.25          1.65        0.81    0.97     1.41     1.60
#>         (0.00)        (0.19)      (0.14)  (0.10)   (0.19)   (0.22)
#>           0.10          0.06        0.00    0.01    -0.04     0.04
#>         (0.00)        (0.00)      (0.00)  (0.00)   (0.00)   (0.00)
#>         503.60        462.94       60.95   89.69   -40.66    71.81
#>         (3.91)       (16.40)     (11.43)  (5.08)  (16.95)  (14.61)
#>         503.60        288.67       19.83   28.18  -214.93   216.26
#>         (3.91)        (5.60)      (2.80)  (0.75)   (7.33)   (7.38)
#>         503.60        209.63       40.27   26.41  -293.97   296.63
#>         (3.91)       (10.89)      (6.25)  (1.13)  (11.34)  (11.52)
#>         503.60        324.80      200.03  151.50  -178.80   260.88
#>         (3.91)       (53.46)     (20.62)  (6.71)  (53.44)  (35.74)
#>         503.60       1374.55      876.77  364.01   870.95  1206.38
#>         (3.91)      (268.79)     (80.64) (88.12) (270.07) (214.06)
#>         503.60        496.45       81.39   61.88    -7.15    86.40
#>         (3.91)       (26.13)     (16.26)  (6.60)  (29.05)  (13.57)
```

``` r
# meta_population <- 
#   declare_estimand(handler = get_meta_population)
#   
# meta_estimands <- 
#   declare_estimand(handler = get_meta_estimands)
# 
# meta_estimators <- 
#   declare_estimator(handler = get_meta_estimators)
# 

meta_design <- meta_population + meta_inquiry + meta_estimators
 
```
