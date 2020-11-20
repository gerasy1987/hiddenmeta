testthat::context("sampling")

testthat::test_that("RDS sampling with various parameters works", {

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
      add_groups = list(p_service = 0.3,
                        loc_1 = 0.3, loc_2 = 0.1, loc_3 = 0.2),

      # RDS parameters
      n_seed = 20,
      n_coupons = 3,
      target_type = "sample",
      target_n_rds = 100,

      # prop sampling parameters
      target_n_prop = 400,

      # TLS sampling parameters
      target_n_tls = 1
    )

  pop <- do.call(what = get_study_population, args = study_1[1:8])

  set.seed(23121987)
  rds_sample <-
    do.call(what = sample_rds,
            args = c(data = list(pop),
                     sampling_variable = "rds",
                     drop_nonsampled = FALSE,
                     study_1[9:12]))

  rds_sample_drop <-
    do.call(what = sample_rds,
            args = c(data = list(pop),
                     sampling_variable = "DRS",
                     drop_nonsampled = TRUE,
                     study_1[9:12]))

  rds_sample_wave <-
    do.call(what = sample_rds,
            args = c(data = list(pop),
                     sampling_variable = "rds",
                     drop_nonsampled = FALSE,
                     n_seed = 40,
                     n_coupons = 2,
                     target_type = "waves",
                     target_n_rds = 2))

  testthat::expect_equal(nrow(rds_sample), study_1$N)
  testthat::expect_equal(names(rds_sample)[grep(pattern = "rds", names(rds_sample))],
                         c("rds", "rds_from", "rds_t", "rds_wave", "rds_hidden"))

  testthat::expect_equal(nrow(rds_sample_drop), study_1$target_n_rds)
  testthat::expect_true(all(!grepl(pattern = "rds", names(rds_sample_drop))))
  testthat::expect_equal(names(rds_sample_drop)[grep(pattern = "DRS", names(rds_sample_drop))],
                         c("DRS", "DRS_from", "DRS_t", "DRS_wave", "DRS_hidden"))

  testthat::expect_true(length(unique(rds_sample_wave$rds_wave)) <= 3)
  testthat::expect_true(all(table(rds_sample_wave$rds_from) <= 2))
  testthat::expect_equal(names(rds_sample_wave)[grep(pattern = "rds", names(rds_sample_wave))],
                         c("rds", "rds_from", "rds_t", "rds_wave", "rds_hidden"))

  testthat::expect_equal(ncol(rds_sample), ncol(rds_sample_drop))
  testthat::expect_equal(dim(rds_sample), dim(rds_sample_wave))
  testthat::expect_equal(ncol(rds_sample), ncol(pop) + 5)

})
