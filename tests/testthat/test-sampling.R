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
                      loc_1 = 0.3, loc_2 = 0.1, loc_3 = 0.2),

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

pop <- do.call(what = get_study_population, args = study_1[1:5])

testthat::test_that("RDS sampling with various parameters works", {

  set.seed(23121987)
  do.call(what = sample_rds,
          args = c(data = list(pop),
                   hidden_var = "hidden",
                   sampling_variable = "rds1",
                   drop_nonsampled = FALSE,
                   study_1[6:9]))

  rds_sample_drop <-
    do.call(what = sample_rds,
          args = c(data = list(pop),
                   hidden_var = "hidden",
                   sampling_variable = "RDS_dr",
                   drop_nonsampled = TRUE,
                   study_1[6:9]))

  do.call(what = sample_rds,
          args = c(data = list(pop),
                   hidden_var = "hidden",
                   sampling_variable = "rds_wave",
                   drop_nonsampled = FALSE,
                   n_seed = 40,
                   n_coupons = 2,
                   n_waves = 3,
                   target_type = "waves",
                   target_n_rds = 70))

  testthat::expect_equal(nrow(pop), study_1$network_handler_args$N)
  testthat::expect_equal(
    names(pop)[grep(pattern = "rds1", names(pop))],
    c("rds1", "rds1_from", "rds1_t", "rds1_wave", "rds1_hidden", "rds1_own_coupon",
      paste0("rds1_coupon_", 1:study_1$n_coupons)))

  testthat::expect_equal(nrow(rds_sample_drop), study_1$target_n_rds)
  testthat::expect_equal(
    names(rds_sample_drop)[grep(pattern = "RDS\\_dr", names(rds_sample_drop))],
    c("RDS_dr", "RDS_dr_from", "RDS_dr_t",
      "RDS_dr_wave", "RDS_dr_hidden", "RDS_dr_own_coupon",
      paste0("RDS_dr_coupon_", 1:study_1$n_coupons)))

  testthat::expect_true(length(stats::na.omit(unique(pop$rds_wave_wave))) <= 3)
  testthat::expect_true(all(table(pop$rds_wave_from)[-1] <= 2))
  testthat::expect_equal(
    names(pop)[grep(pattern = "rds\\_wave", names(pop))],
    c("rds_wave", "rds_wave_from", "rds_wave_t",
      "rds_wave_wave", "rds_wave_hidden", "rds_wave_own_coupon",
      paste0("rds_wave_coupon_", 1:2)))

  testthat::expect_equal(ncol(pop), ncol(rds_sample_drop))
  testthat::expect_equal(ncol(rds_sample_drop), ncol(pop) - 6 - 2)

})

testthat::test_that("PPS sampling with various parameters works", {

  set.seed(23121987)
  pps_sample <-
    do.call(what = sample_pps,
            args = c(data = list(pop),
                     sampling_variable = "pps",
                     drop_nonsampled = FALSE,
                     study_1[10]))

  pps_sample_drop <-
    do.call(what = sample_pps,
            args = c(data = list(pop),
                     sampling_variable = "SPS",
                     drop_nonsampled = TRUE,
                     study_1[10]))

  pps_sample_other <-
    do.call(what = sample_pps,
            args = c(data = list(pop),
                     sampling_variable = "rds",
                     drop_nonsampled = FALSE,
                     target_n_pps = 800))

  testthat::expect_equal(nrow(pps_sample), study_1$network_handler_args$N)
  testthat::expect_equal(names(pps_sample)[grep(pattern = "pps", names(pps_sample))],
                         c("pps_share", "pps"))

  testthat::expect_equal(nrow(pps_sample_drop), study_1$target_n_pps)
  testthat::expect_true(all(!grepl(pattern = "pps", names(pps_sample_drop))))
  testthat::expect_equal(names(pps_sample_drop)[grep(pattern = "SPS", names(pps_sample_drop))],
                         c("SPS_share", "SPS"))

  testthat::expect_true(sum(pps_sample_other$rds == 1) == 800)
  testthat::expect_equal(names(pps_sample_other)[grep(pattern = "rds", names(pps_sample_other))],
                         c("rds_share", "rds"))

  testthat::expect_equal(ncol(pps_sample), ncol(pps_sample_drop))
  testthat::expect_equal(dim(pps_sample), dim(pps_sample_other))
  testthat::expect_equal(ncol(pps_sample), ncol(pop) + 2)

  testthat::expect_error(
    do.call(what = sample_pps,
            args = c(data = list(pop),
                     sampling_variable = "rds",
                     drop_nonsampled = FALSE,
                     target_n_pps = (study_1$network_handler_args$N + 1))),
    regexp = "Requested PPS sample size is larger than population size")

})

testthat::test_that("TLS sampling with various parameters works", {

  set.seed(23121987)
  tls_sample <-
    do.call(what = sample_tls,
            args = c(data = list(pop),
                     sampling_variable = "tls",
                     drop_nonsampled = FALSE,
                     study_1[11]))

  tls_sample_drop <-
    do.call(what = sample_tls,
            args = c(data = list(pop),
                     sampling_variable = "PPS",
                     drop_nonsampled = TRUE,
                     study_1[11]))

  tls_sample_other <-
    do.call(what = sample_tls,
            args = c(data = list(pop),
                     sampling_variable = "rds",
                     drop_nonsampled = FALSE,
                     target_n_tls = 1))

  testthat::expect_equal(nrow(tls_sample), study_1$network_handler_args$N)
  testthat::expect_equal(names(tls_sample)[grep(pattern = "tls", names(tls_sample))],
                         c("tls_loc_sampled", "tls"))

  testthat::expect_true(nrow(tls_sample_drop) <= study_1$network_handler_args$N)
  testthat::expect_true(all(!grepl(pattern = "tls", names(tls_sample_drop))))
  testthat::expect_equal(names(tls_sample_drop)[grep(pattern = "PPS", names(tls_sample_drop))],
                         c("PPS_loc_sampled", "PPS"))

  testthat::expect_true(sum(tls_sample_other$rds == 1, na.rm = T) <=
                          sum(tls_sample_other$rds_loc_sampled == "loc_1", na.rm = TRUE))
  testthat::expect_true(all(!grepl(pattern = "tls", names(tls_sample_other))))
  testthat::expect_length(
    object = unique(tls_sample_other$rds_loc_sampled[tls_sample_other$rds == 1]),
    n = 2
  )

  testthat::expect_equal(ncol(tls_sample), ncol(tls_sample_drop))
  testthat::expect_equal(dim(tls_sample), dim(tls_sample_other))
  testthat::expect_equal(ncol(tls_sample), ncol(pop) + 2)

  testthat::expect_error(
    do.call(what = sample_tls,
            args = c(data = list(pop),
                     sampling_variable = "rds",
                     drop_nonsampled = FALSE,
                     target_n_tls = 4)),
    regexp = "Number of requested locations for TLS sample exceeds number of available locations")

})
