testthat::context("population")

testthat::test_that("generation of population with two groups works", {

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

  testthat::expect_equal(
    dim(do.call(what = get_study_population, args = study_1[1:8])),
    c(study_1$N, 18))

})
