testthat::context("sampling")

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
    p_visible = list(known = .99, hidden = .7),

    # probability of service utilization in hidden population
    # for service multiplier
    add_groups = list(p_service = 0.3,
                      "paste0('loc_', 1:10) := lapply(rep(.2, times = 10), function(add) rbinom(.N, 1, 0.05 + hidden * add))"),

    # RDS parameters
    n_seed = 20,
    n_coupons = 3,
    target_type = "sample",
    target_n_rds = 100,

    # prop sampling parameters
    target_n_pps = 400,

    # TLS sampling parameters
    target_n_clusters = 5,
    target_cluster_type = "fixed",
    target_per_cluster = 15,
    clusters = paste0("loc_", 1:10)
  )

pop <- do.call(what = get_study_population, args = study_1[1:5])

testthat::test_that("RDS sampling with various parameters works", {

  pop_rds <- data.table::copy(pop)

  set.seed(23121987)
  rds_sample <-
    do.call(what = sample_rds,
          args = c(data = list(pop_rds),
                   hidden_var = "hidden",
                   sampling_variable = "rds1",
                   drop_nonsampled = FALSE,
                   study_1[6:9]))

  rds_sample_drop <-
    do.call(what = sample_rds,
          args = c(data = list(pop_rds),
                   hidden_var = "hidden",
                   sampling_variable = "RDS_dr",
                   drop_nonsampled = TRUE,
                   study_1[6:9]))

  rds_sample_other <-
    do.call(what = sample_rds,
            args = list(data = pop_rds,
                        hidden_var = "hidden",
                        sampling_variable = "tls",
                        drop_nonsampled = FALSE,
                        n_seed = 40,
                        n_coupons = 2,
                        n_waves = 3,
                        target_type = "waves",
                        target_n_rds = 70))

  # test that data was not subset
  testthat::expect_equal(nrow(rds_sample), study_1$network_handler_args$N)

  # test sampling column names
  testthat::expect_equal(
    names(rds_sample)[grep(pattern = "rds1", names(rds_sample))],
    paste0("rds1", c("", "_from", "_t", "_wave", "_hidden", "_own_coupon",
                     paste0("_coupon_", 1:study_1$n_coupons))))

  # test that data.table does not auto-update populate data
  testthat::expect_length(
    names(rds_sample_drop)[grep(pattern = "rds", names(rds_sample_drop))],
    n = 0)

  # test that drop_nonsampled works
  testthat::expect_equal(nrow(rds_sample_drop), study_1$target_n_rds)

  # test sampling column names
  testthat::expect_equal(
    names(rds_sample_drop)[grep(pattern = "RDS\\_dr", names(rds_sample_drop))],
    paste0("RDS_dr", c("", "_from", "_t", "_wave", "_hidden", "_own_coupon",
                       paste0("_coupon_", 1:study_1$n_coupons))))

  # test that specifying number of waves works
  testthat::expect_true(
    length(na.omit(unique(rds_sample_other$rds_wave_wave))) <= 3)
  testthat::expect_true(all(table(rds_sample_other$rds_wave_from)[-1] <= 2))

  # test that even naming rds sample tls works
  testthat::expect_equal(
    names(rds_sample_other)[grep(pattern = "tls", names(rds_sample_other))],
    paste0("tls", c("", "_from", "_t", "_wave", "_hidden", "_own_coupon",
                     paste0("_coupon_", 1:2))))

  # test number of columns
  testthat::expect_equal(ncol(rds_sample), ncol(rds_sample_drop))
  testthat::expect_equal(ncol(rds_sample), ncol(rds_sample_other) + 1)
  testthat::expect_true(ncol(pop) <= ncol(rds_sample))

})

testthat::test_that("PPS sampling with various parameters works", {

  pop_pps <- data.table::copy(pop)

  set.seed(23121987)
  pps_sample <-
    do.call(what = sample_pps,
            args = c(data = list(pop_pps),
                     sampling_variable = "pps",
                     drop_nonsampled = FALSE,
                     study_1[10]))

  pps_sample_drop <-
    do.call(what = sample_pps,
            args = c(data = list(pop_pps),
                     sampling_variable = "SPS",
                     drop_nonsampled = TRUE,
                     study_1[10]))

  pps_sample_other <-
    do.call(what = sample_pps,
            args = list(data = pop_pps,
                        sampling_variable = "rds",
                        drop_nonsampled = FALSE,
                        target_n_pps = 800))


  # test that data was not subset
  testthat::expect_equal(nrow(pps_sample), study_1$network_handler_args$N)

  # test sampling column names
  testthat::expect_equal(
    names(pps_sample)[grep(pattern = "^pps", names(pps_sample))],
    paste0("pps", c("", "_frame", "_cluster_prop",
                    "_strata_prop", "_sampled_cluster",
                    "_weight", "_cluster", "_strata")))

  # test that drop_nonsampled works
  testthat::expect_equal(nrow(pps_sample_drop), study_1$target_n_pps)

  # test sampling column names
  testthat::expect_equal(
    names(pps_sample_drop)[grep(pattern = "SPS", names(pps_sample_drop))],
    paste0("SPS", c("", "_frame", "_cluster_prop",
                    "_strata_prop", "_sampled_cluster",
                    "_weight", "_cluster", "_strata")))

  # test that changing sample size works
  testthat::expect_true(sum(pps_sample_other$rds == 1) == 800)

  # test that even naming pps sample rds works
  testthat::expect_equal(
    names(pps_sample_other)[grep(pattern = "rds", names(pps_sample_other))],
    paste0("rds", c("", "_frame", "_cluster_prop",
                    "_strata_prop", "_sampled_cluster",
                    "_weight", "_cluster", "_strata")))

  # test number of columns
  testthat::expect_equal(ncol(pps_sample), ncol(pps_sample_drop))
  testthat::expect_equal(dim(pps_sample_other), dim(pps_sample))
  testthat::expect_true(ncol(pop) <= ncol(pps_sample_other))

  # test that requested sample cannot be larger than population
  testthat::expect_error(
    do.call(what = sample_pps,
            args = list(data = pop,
                        sampling_variable = "rds",
                        drop_nonsampled = FALSE,
                        target_n_pps = (study_1$network_handler_args$N + 1))),
    regexp =
      "Requested PPS sample size is larger than population (or sampling frame) size",
    fixed = TRUE)

})

testthat::test_that("TLS sampling with various parameters works", {

  pop_tls <- data.table::copy(pop)

  set.seed(23121987)
  tls_sample <-
    do.call(what = sample_tls,
            args = c(data = list(pop_tls),
                     sampling_variable = "tls",
                     drop_nonsampled = FALSE,
                     study_1[11:14]))

  tls_sample_drop <-
    do.call(what = sample_tls,
            args = c(data = list(pop_tls),
                     sampling_variable = "PPS",
                     drop_nonsampled = TRUE,
                     study_1[11:14]))

  tls_sample_other <-
    do.call(what = sample_tls,
            args = list(data = pop_tls,
                        sampling_variable = "rds",
                        drop_nonsampled = FALSE,
                        target_n_clusters = 3,
                        target_cluster_type = "prop",
                        target_per_cluster = .5,
                        clusters = paste0("loc_", 1:9)))

  # test that data was not subset
  testthat::expect_equal(nrow(tls_sample), study_1$network_handler_args$N)

  # test sampling column names
  testthat::expect_equal(
    names(tls_sample)[grep(pattern = "tls", names(tls_sample))],
    paste0("tls", c("", "_loc_present","_loc_sampled", "_weight",
                    "_weight_visible","_locs_sampled")))

  # test that drop_nonsampled works
  testthat::expect_true(nrow(tls_sample_drop) <= study_1$network_handler_args$N)

  # test that data.table does not auto-update populate data
  testthat::expect_length(
    names(tls_sample_drop)[grep(pattern = "tls", names(tls_sample_drop))],
    n = 0)

  # test non-standard sampling column names
  testthat::expect_equal(
    names(tls_sample_drop)[grep(pattern = "PPS", names(tls_sample_drop))],
    paste0("PPS", c("", "_loc_present", "_loc_sampled",
                    "_weight","_weight_visible", "_locs_sampled")))

  # test that sampling reacts to changes in number of clusters
  testthat::expect_length(
    strsplit(na.omit(unique(tls_sample$tls_locs_sampled)), split = ";")[[1]],
    n = 5)
  testthat::expect_length(
    strsplit(na.omit(unique(tls_sample_other$rds_locs_sampled)), split = ";")[[1]],
    n = 3)

  # test that sampling columns created properly even if rds is used instead of tls
  testthat::expect_equal(
    names(tls_sample_other)[grep(pattern = "rds", names(tls_sample_other))],
    paste0("rds", c("", "_loc_present", "_loc_sampled",
                    "_weight","_weight_visible", "_locs_sampled")))

  # test number of columns
  testthat::expect_equal(ncol(tls_sample), ncol(tls_sample_drop))
  testthat::expect_equal(dim(tls_sample), dim(tls_sample_other))
  testthat::expect_true(ncol(pop) <= ncol(tls_sample))

  # test that target number of clusters is properly checked
  testthat::expect_error(
    do.call(what = sample_tls,
            args = list(data = pop,
                        sampling_variable = "rds",
                        drop_nonsampled = FALSE,
                        target_n_clusters = 20,
                        target_cluster_type = "prop",
                        target_per_cluster = .5,
                        clusters = paste0("loc_", 1:9))),
    regexp = "Number of requested locations for TLS sample exceeds number of available locations")

})
