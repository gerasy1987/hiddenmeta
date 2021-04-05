
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
        add_groups = list(p_service = 0.3, 
                          loc_1 = 0.3, loc_2 = 0.1, loc_3 = 0.2, 
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
                   target_n_clusters = 2,
                   target_n_tls = 100,
                   cluster = paste0("loc_", 1:3)),
        pps = list(handler = sample_pps,
                   sampling_variable = "pps",
                   # prop sampling parameters
                   sampling_frame = NULL,
                   strata = NULL,
                   cluster = NULL,
                   target_n_pps = 400)
      ),
    inquiries = list(handler = get_study_estimands),
    estimators = 
      list(
        rds = 
          list(sspse = list(handler = get_study_est_sspse, 
                            label = "rds_sspse",
                            prior_mean = 200,
                            mcmc_params = list(interval = 5, burnin = 2000, samplesize = 500),
                            rds_prefix = "rds")),
        # tls = 
        #   list(nsum = list(handler = get_study_est_nsum,
        #                    prefix = "tls",
        #                    known = c("known", "known_2", "known_3"),
        #                    hidden = "hidden_visible_out",
        #                    label = "tls_nsum")),
        pps = 
          list(ht = list(handler = get_study_est_ht, 
                         label = "pps_ht"),
               nsum = list(handler = get_study_est_nsum,
                           known = c("known", "known_2", "known_3"),
                           hidden = "hidden_visible_out",
                           label = "pps_nsum"))
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
#> # A tibble: 2,000 x 47
#>     name type  known hidden links    p_service loc_1 loc_2 loc_3 known_2 known_3
#>    <int> <chr> <int>  <int> <list>       <int> <int> <int> <int>   <int>   <int>
#>  1     1 00        0      0 <igrph.…         0     1     0     0       0       1
#>  2     2 00        0      0 <igrph.…         0     0     0     1       0       0
#>  3     3 00        0      0 <igrph.…         1     0     1     0       0       0
#>  4     4 00        0      0 <igrph.…         0     0     0     0       1       0
#>  5     5 00        0      0 <igrph.…         0     1     0     1       0       1
#>  6     6 00        0      0 <igrph.…         1     0     0     0       0       0
#>  7     7 00        0      0 <igrph.…         1     0     0     0       0       0
#>  8     8 00        0      0 <igrph.…         1     0     0     0       0       0
#>  9     9 00        0      0 <igrph.…         0     1     0     0       0       0
#> 10    10 00        0      0 <igrph.…         0     0     0     0       1       0
#> # … with 1,990 more rows, and 36 more variables: n_visible_out <dbl>,
#> #   known_visible_out <dbl>, hidden_visible_out <dbl>,
#> #   type_00_visible_out <dbl>, type_01_visible_out <dbl>,
#> #   type_10_visible_out <dbl>, type_11_visible_out <dbl>,
#> #   p_service_visible_out <dbl>, loc_1_visible_out <dbl>,
#> #   loc_2_visible_out <dbl>, loc_3_visible_out <dbl>,
#> #   known_2_visible_out <dbl>, known_3_visible_out <dbl>,
#> #   known_visible_in <dbl>, hidden_visible_in <dbl>, type_00_visible_in <dbl>,
#> #   type_01_visible_in <dbl>, type_10_visible_in <dbl>,
#> #   type_11_visible_in <dbl>, p_service_visible_in <dbl>,
#> #   loc_1_visible_in <dbl>, loc_2_visible_in <dbl>, loc_3_visible_in <dbl>,
#> #   known_2_visible_in <dbl>, known_3_visible_in <dbl>, p_visible_known <dbl>,
#> #   p_visible_hidden <dbl>, total <int>, total_known <int>, total_hidden <int>,
#> #   total_p_service <int>, total_loc_1 <int>, total_loc_2 <int>,
#> #   total_loc_3 <int>, total_known_2 <int>, total_known_3 <int>
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
#> # A tibble: 2,000 x 56
#>     name type  known hidden links    p_service loc_1 loc_2 loc_3 known_2 known_3
#>    <int> <chr> <int>  <int> <list>       <int> <int> <int> <int>   <int>   <int>
#>  1     1 00        0      0 <igrph.…         0     1     0     0       0       1
#>  2     2 00        0      0 <igrph.…         0     0     0     1       0       0
#>  3     3 00        0      0 <igrph.…         1     0     1     0       0       0
#>  4     4 00        0      0 <igrph.…         0     0     0     0       1       0
#>  5     5 00        0      0 <igrph.…         0     1     0     1       0       1
#>  6     6 00        0      0 <igrph.…         1     0     0     0       0       0
#>  7     7 00        0      0 <igrph.…         1     0     0     0       0       0
#>  8     8 00        0      0 <igrph.…         1     0     0     0       0       0
#>  9     9 00        0      0 <igrph.…         0     1     0     0       0       0
#> 10    10 00        0      0 <igrph.…         0     0     0     0       1       0
#> # … with 1,990 more rows, and 45 more variables: n_visible_out <dbl>,
#> #   known_visible_out <dbl>, hidden_visible_out <dbl>,
#> #   type_00_visible_out <dbl>, type_01_visible_out <dbl>,
#> #   type_10_visible_out <dbl>, type_11_visible_out <dbl>,
#> #   p_service_visible_out <dbl>, loc_1_visible_out <dbl>,
#> #   loc_2_visible_out <dbl>, loc_3_visible_out <dbl>,
#> #   known_2_visible_out <dbl>, known_3_visible_out <dbl>,
#> #   known_visible_in <dbl>, hidden_visible_in <dbl>, type_00_visible_in <dbl>,
#> #   type_01_visible_in <dbl>, type_10_visible_in <dbl>,
#> #   type_11_visible_in <dbl>, p_service_visible_in <dbl>,
#> #   loc_1_visible_in <dbl>, loc_2_visible_in <dbl>, loc_3_visible_in <dbl>,
#> #   known_2_visible_in <dbl>, known_3_visible_in <dbl>, p_visible_known <dbl>,
#> #   p_visible_hidden <dbl>, total <int>, total_known <int>, total_hidden <int>,
#> #   total_p_service <int>, total_loc_1 <int>, total_loc_2 <int>,
#> #   total_loc_3 <int>, total_known_2 <int>, total_known_3 <int>, rds <int>,
#> #   rds_from <dbl>, rds_t <dbl>, rds_wave <dbl>, rds_hidden <int>,
#> #   rds_own_coupon <chr>, rds_coupon_1 <chr>, rds_coupon_2 <chr>,
#> #   rds_coupon_3 <chr>
```

``` r
study_sample_pps <- 
  do.call(what = declare_sampling, 
          args = study_1$sample$pps)

set.seed(19872312)
draw_data(study_population + study_sample_rds + study_sample_pps)
#> # A tibble: 2,000 x 62
#>     name type  known hidden links    p_service loc_1 loc_2 loc_3 known_2 known_3
#>    <int> <chr> <int>  <int> <list>       <int> <int> <int> <int>   <int>   <int>
#>  1     1 00        0      0 <igrph.…         0     1     0     0       0       1
#>  2     2 00        0      0 <igrph.…         0     0     0     1       0       0
#>  3     3 00        0      0 <igrph.…         1     0     1     0       0       0
#>  4     4 00        0      0 <igrph.…         0     0     0     0       1       0
#>  5     5 00        0      0 <igrph.…         0     1     0     1       0       1
#>  6     6 00        0      0 <igrph.…         1     0     0     0       0       0
#>  7     7 00        0      0 <igrph.…         1     0     0     0       0       0
#>  8     8 00        0      0 <igrph.…         1     0     0     0       0       0
#>  9     9 00        0      0 <igrph.…         0     1     0     0       0       0
#> 10    10 00        0      0 <igrph.…         0     0     0     0       1       0
#> # … with 1,990 more rows, and 51 more variables: n_visible_out <dbl>,
#> #   known_visible_out <dbl>, hidden_visible_out <dbl>,
#> #   type_00_visible_out <dbl>, type_01_visible_out <dbl>,
#> #   type_10_visible_out <dbl>, type_11_visible_out <dbl>,
#> #   p_service_visible_out <dbl>, loc_1_visible_out <dbl>,
#> #   loc_2_visible_out <dbl>, loc_3_visible_out <dbl>,
#> #   known_2_visible_out <dbl>, known_3_visible_out <dbl>,
#> #   known_visible_in <dbl>, hidden_visible_in <dbl>, type_00_visible_in <dbl>,
#> #   type_01_visible_in <dbl>, type_10_visible_in <dbl>,
#> #   type_11_visible_in <dbl>, p_service_visible_in <dbl>,
#> #   loc_1_visible_in <dbl>, loc_2_visible_in <dbl>, loc_3_visible_in <dbl>,
#> #   known_2_visible_in <dbl>, known_3_visible_in <dbl>, p_visible_known <dbl>,
#> #   p_visible_hidden <dbl>, total <int>, total_known <int>, total_hidden <int>,
#> #   total_p_service <int>, total_loc_1 <int>, total_loc_2 <int>,
#> #   total_loc_3 <int>, total_known_2 <int>, total_known_3 <int>, rds <int>,
#> #   rds_from <dbl>, rds_t <dbl>, rds_wave <dbl>, rds_hidden <int>,
#> #   rds_own_coupon <chr>, rds_coupon_1 <chr>, rds_coupon_2 <chr>,
#> #   rds_coupon_3 <chr>, pps <dbl>, pps_frame <dbl>, pps_cluster_id <int>,
#> #   pps_cluster_prop <dbl>, pps_weight <dbl>, pps_sampled_cluster <dbl>
```

``` r
study_sample_tls <- 
  do.call(what = declare_sampling, 
          args = study_1$sample$tls)

set.seed(19872312)
draw_data(study_population + 
            study_sample_rds + study_sample_pps + study_sample_tls)
#> # A tibble: 2,000 x 66
#>     name type  known hidden links    p_service loc_1 loc_2 loc_3 known_2 known_3
#>    <int> <chr> <int>  <int> <list>       <int> <int> <int> <int>   <int>   <int>
#>  1     1 00        0      0 <igrph.…         0     1     0     0       0       1
#>  2     2 00        0      0 <igrph.…         0     0     0     1       0       0
#>  3     3 00        0      0 <igrph.…         1     0     1     0       0       0
#>  4     4 00        0      0 <igrph.…         0     0     0     0       1       0
#>  5     5 00        0      0 <igrph.…         0     1     0     1       0       1
#>  6     6 00        0      0 <igrph.…         1     0     0     0       0       0
#>  7     7 00        0      0 <igrph.…         1     0     0     0       0       0
#>  8     8 00        0      0 <igrph.…         1     0     0     0       0       0
#>  9     9 00        0      0 <igrph.…         0     1     0     0       0       0
#> 10    10 00        0      0 <igrph.…         0     0     0     0       1       0
#> # … with 1,990 more rows, and 55 more variables: n_visible_out <dbl>,
#> #   known_visible_out <dbl>, hidden_visible_out <dbl>,
#> #   type_00_visible_out <dbl>, type_01_visible_out <dbl>,
#> #   type_10_visible_out <dbl>, type_11_visible_out <dbl>,
#> #   p_service_visible_out <dbl>, loc_1_visible_out <dbl>,
#> #   loc_2_visible_out <dbl>, loc_3_visible_out <dbl>,
#> #   known_2_visible_out <dbl>, known_3_visible_out <dbl>,
#> #   known_visible_in <dbl>, hidden_visible_in <dbl>, type_00_visible_in <dbl>,
#> #   type_01_visible_in <dbl>, type_10_visible_in <dbl>,
#> #   type_11_visible_in <dbl>, p_service_visible_in <dbl>,
#> #   loc_1_visible_in <dbl>, loc_2_visible_in <dbl>, loc_3_visible_in <dbl>,
#> #   known_2_visible_in <dbl>, known_3_visible_in <dbl>, p_visible_known <dbl>,
#> #   p_visible_hidden <dbl>, total <int>, total_known <int>, total_hidden <int>,
#> #   total_p_service <int>, total_loc_1 <int>, total_loc_2 <int>,
#> #   total_loc_3 <int>, total_known_2 <int>, total_known_3 <int>, rds <int>,
#> #   rds_from <dbl>, rds_t <dbl>, rds_wave <dbl>, rds_hidden <int>,
#> #   rds_own_coupon <chr>, rds_coupon_1 <chr>, rds_coupon_2 <chr>,
#> #   rds_coupon_3 <chr>, pps <dbl>, pps_frame <dbl>, pps_cluster_id <int>,
#> #   pps_cluster_prop <dbl>, pps_weight <dbl>, pps_sampled_cluster <dbl>,
#> #   tls <dbl>, tls_loc_sampled <chr>, tls_sampled_locs <chr>, tls_weight <dbl>
```

## Step 4. Declare study level estimands

``` r
study_estimands <- 
  do.call(what = declare_inquiry, 
          args = study_1$inquiries)

set.seed(19872312)
draw_inquiry(study_population + 
                study_sample_rds + study_sample_pps + study_sample_tls + 
                study_estimands)
#>           inquiry_label estimand
#> 1           hidden_size 190.0000
#> 2           hidden_prev   0.0950
#> 3        degree_average   4.9310
#> 4 degree_hidden_average   0.7155
```

## Step 5. Declare estimators used in the study

``` r
estimator_sspse <- do.call(what = declare_estimator, 
                           args = study_1$estimators$rds$sspse)
estimator_ht <- do.call(what = declare_estimator, 
                        args = study_1$estimators$pps$ht)
# estimator_nsum_tls <- do.call(what = declare_estimator, 
#                               args = study_1$estimators$tls$nsum)
estimator_nsum <- do.call(what = declare_estimator, 
                          args = study_1$estimators$pps$nsum)

set.seed(19872312)
draw_estimates(study_population +
                 study_sample_rds + study_sample_pps + study_sample_tls +
                 study_estimands +
                 estimator_sspse + estimator_ht + #estimator_nsum_tls + 
                 estimator_nsum)
#>           estimator_label  estimate           se  inquiry_label
#> 1   hidden_size_rds_sspse 221.50000 121.30762406    hidden_size
#> 2      hidden_prev_pps_ht   0.09750   0.01485044    hidden_prev
#> 3      hidden_size_pps_ht 195.00000  29.70088984    hidden_size
#> 4    hidden_size_pps_nsum 190.07131  13.08324953    hidden_size
#> 5 degree_average_pps_nsum   5.02443           NA degree_average
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
      study_sample_rds + study_sample_pps + 
      study_estimands + 
      estimator_sspse + estimator_ht + estimator_nsum, 
    sims = 10)

diagnose_design(study_simulations, 
                diagnosands = study_diagnosands, sims = 2)
#> 
#> Research design diagnosis based on 10 simulations. Diagnosand estimates with bootstrapped standard errors in parentheses (100 replicates).
#> 
#>  Design Label         Inquiry Label         Estimator Label N Sims
#>      design_1        degree_average degree_average_pps_nsum     10
#>                                                                   
#>      design_1 degree_hidden_average                    <NA>     10
#>                                                                   
#>      design_1           hidden_prev      hidden_prev_pps_ht     10
#>                                                                   
#>      design_1           hidden_size      hidden_size_pps_ht     10
#>                                                                   
#>      design_1           hidden_size    hidden_size_pps_nsum     10
#>                                                                   
#>      design_1           hidden_size   hidden_size_rds_sspse     10
#>                                                                   
#>  Mean Estimand Mean Estimate SD Estimate Mean Se   Bias   RMSE
#>           4.91          4.89        0.21      NA  -0.02   0.14
#>         (0.02)        (0.07)      (0.03)      NA (0.05) (0.02)
#>           0.78            NA          NA      NA     NA     NA
#>         (0.04)            NA          NA      NA     NA     NA
#>           0.10          0.10        0.01    0.01  -0.00   0.01
#>         (0.00)        (0.00)      (0.00)  (0.00) (0.00) (0.00)
#>         197.50        195.00       15.99   29.67  -2.50  15.73
#>         (5.19)        (5.09)      (3.09)  (0.35) (4.76) (2.59)
#>         197.50        211.86       15.11   13.91  14.36  18.59
#>         (5.19)        (4.99)      (3.14)  (0.16) (3.89) (3.84)
#>         197.50        165.05       16.69   65.42 -32.45  42.92
#>         (5.19)        (5.40)      (2.43)  (6.07) (9.77) (6.90)
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
                   target_n_pps = 1000)
      ),
    inquiries = list(handler = get_study_estimands),
    estimators = 
      list(
        rds = 
          list(sspse = list(handler = get_study_est_sspse, 
                            label = "rds_sspse",
                            prior_mean = 450,
                            mcmc_params = list(interval = 5, burnin = 2000, samplesize = 500),
                            rds_prefix = "rds")),
        pps = 
          list(
            ht = list(handler = get_study_est_ht, 
                      label = "pps_ht"),
            nsum = list(handler = get_study_est_nsum,
                        known = c("frame", "known"),
                        hidden = "hidden_visible_out",
                        label = "pps_nsum")
          )
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
#> # A tibble: 5 x 3
#>   study   sample population                
#>   <chr>   <chr>  <list>                    
#> 1 study_1 rds    <tibble[,56] [2,000 × 56]>
#> 2 study_1 tls    <tibble[,51] [2,000 × 51]>
#> 3 study_1 pps    <tibble[,53] [2,000 × 53]>
#> 4 study_2 rds    <tibble[,61] [5,000 × 61]>
#> 5 study_2 pps    <tibble[,58] [5,000 × 58]>
draw_inquiry(multi_design)
#>                   inquiry_label estimand
#> 1           study_1_hidden_size 215.0000
#> 2           study_1_hidden_prev   0.1075
#> 3        study_1_degree_average   4.9100
#> 4 study_1_degree_hidden_average   0.9235
#> 5           study_2_hidden_size 488.0000
#> 6           study_2_hidden_prev   0.0976
#> 7        study_2_degree_average   4.3332
#> 8 study_2_degree_hidden_average   0.2374
draw_estimates(multi_design)
#>            estimator_label   estimate           se          inquiry_label
#> 1    hidden_size_rds_sspse 160.500000 5.619659e+01    study_1_hidden_size
#> 2       hidden_prev_pps_ht   0.097500 1.485044e-02    study_1_hidden_prev
#> 3       hidden_size_pps_ht 195.000000 2.970089e+01    study_1_hidden_size
#> 4     hidden_size_pps_nsum 214.594203 1.363877e+01    study_1_hidden_size
#> 5  degree_average_pps_nsum   5.149254           NA study_1_degree_average
#> 6    hidden_size_rds_sspse 535.500000 6.421872e+02    study_2_hidden_size
#> 7       hidden_prev_pps_ht   0.102000 9.575369e-03    study_2_hidden_prev
#> 8       hidden_size_pps_ht 510.000000 4.787684e+01    study_2_hidden_size
#> 9     hidden_size_pps_nsum 170.831554 1.403949e+01    study_2_hidden_size
#> 10 degree_average_pps_nsum   4.185409           NA study_2_degree_average

multi_simulations <- 
  simulate_design(multi_population +
                    multi_sampling + 
                    multi_inquiry + 
                    multi_estimators, sims = 10)

diagnose_design(simulations_df = multi_simulations, 
                diagnosands = study_diagnosands)
#> 
#> Research design diagnosis based on 10 simulations. Diagnosand estimates with bootstrapped standard errors in parentheses (100 replicates).
#> 
#>  Design Label                 Inquiry Label         Estimator Label N Sims
#>      design_1        study_1_degree_average degree_average_pps_nsum     10
#>                                                                           
#>      design_1 study_1_degree_hidden_average                    <NA>     10
#>                                                                           
#>      design_1           study_1_hidden_prev      hidden_prev_pps_ht     10
#>                                                                           
#>      design_1           study_1_hidden_size      hidden_size_pps_ht     10
#>                                                                           
#>      design_1           study_1_hidden_size    hidden_size_pps_nsum     10
#>                                                                           
#>      design_1           study_1_hidden_size   hidden_size_rds_sspse     10
#>                                                                           
#>      design_1        study_2_degree_average degree_average_pps_nsum     10
#>                                                                           
#>      design_1 study_2_degree_hidden_average                    <NA>     10
#>                                                                           
#>      design_1           study_2_hidden_prev      hidden_prev_pps_ht     10
#>                                                                           
#>      design_1           study_2_hidden_size      hidden_size_pps_ht     10
#>                                                                           
#>      design_1           study_2_hidden_size    hidden_size_pps_nsum     10
#>                                                                           
#>      design_1           study_2_hidden_size   hidden_size_rds_sspse     10
#>                                                                           
#>  Mean Estimand Mean Estimate SD Estimate Mean Se     Bias     RMSE
#>           4.93          4.94        0.19      NA     0.01     0.18
#>         (0.02)        (0.06)      (0.04)      NA   (0.06)   (0.04)
#>           0.79            NA          NA      NA       NA       NA
#>         (0.02)            NA          NA      NA       NA       NA
#>           0.10          0.10        0.02    0.01    -0.00     0.02
#>         (0.00)        (0.01)      (0.00)  (0.00)   (0.00)   (0.00)
#>         200.80        195.00       36.67   29.56    -5.80    30.35
#>         (1.91)       (10.30)      (9.70)  (0.67)   (9.16)   (6.34)
#>         200.80        218.59       42.30   13.97    17.79    39.32
#>         (1.91)       (11.16)      (7.83)  (0.32)  (10.02)   (6.54)
#>         200.80        173.70       43.82   74.35   -27.10    53.39
#>         (1.91)       (11.87)      (9.41) (10.95)  (12.93)   (7.43)
#>           4.22          4.10        0.07      NA    -0.12     0.13
#>         (0.02)        (0.02)      (0.01)      NA   (0.01)   (0.01)
#>           0.24            NA          NA      NA       NA       NA
#>         (0.01)            NA          NA      NA       NA       NA
#>           0.10          0.10        0.01    0.01     0.00     0.01
#>         (0.00)        (0.00)      (0.00)  (0.00)   (0.00)   (0.00)
#>         499.40        505.00       45.34   47.61     5.60    36.25
#>         (6.53)       (13.16)      (7.39)  (0.56)  (10.99)   (5.83)
#>         499.40        182.01       35.46   14.56  -317.39   318.31
#>         (6.53)       (11.03)      (5.31)  (0.47)   (7.89)   (7.93)
#>         499.40       1761.60      702.87  298.24  1262.20  1430.97
#>         (6.53)      (221.56)    (197.94) (90.26) (224.41) (140.81)
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
