testthat::context("population")

testthat::test_that("generation of population with two groups works", {

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

  pop <- do.call(what = get_study_population, args = study_1[1:5])

  testthat::expect_equal(nrow(pop), study_1$network_handler_args$N)

  testthat::expect_true(all(c("known","hidden","links") %in% names(pop)))

  testthat::expect_false("known_1" %in% names(pop))

})

testthat::test_that("generation of population with more than two groups works", {

  ## STUDY 2
  study_2 <-
    list(
      network_handler = sim_block_network,
      network_handler_args =
        list(N = 3000,
             K = 3,
             prev_K = c(known_1 = .2, known_2 = .3, hidden = .1),
             rho_K = c(.05),
             p_edge_within = list(known_1 = c(0.05, 0.05),
                                  known_2 = c(0.05, 0.05),
                                  hidden = c(0.05, 0.8)),
             p_edge_between = list(known_1 = 0.05,
                                   known_2 = 0.05,
                                   hidden = 0.01),
             directed = FALSE),

      group_names = c("known_1", "known_2", "hidden"),
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

  pop <- do.call(what = get_study_population, args = study_2[1:5])

  testthat::expect_equal(nrow(pop), study_2$network_handler_args$N)

  testthat::expect_equal(length(grep(pattern = "^hidden|^known_", names(pop))),
                         2 * study_2$network_handler_args$K)

  testthat::expect_equal(length(unique(pop$type)),
                         length(grep(pattern = "^type_visible", names(pop))))

})
