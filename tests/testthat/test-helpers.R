testthat::context("helpers")

testthat::test_that("generation of group sizes works", {

  N <- 1000
  K <- 3
  prev_K <- c(known1 = .3, known2 = .3, hidden = .1)
  rho_K <- c(.05,.05,.01)

  type_names <-
    lapply(1:K, function(x) 0:1) %>%
    expand.grid() %>%
    apply(X = ., MARGIN = 1, FUN = paste0, collapse = "")

  testthat::expect_equal(length(gen_group_sizes(N, prev_K, rho_K)), length(type_names))
  testthat::expect_equal(names(gen_group_sizes(N, prev_K, rho_K)), type_names)

  testthat::expect_error(gen_group_sizes(N, prev_K, rho_K = c(.5, .1)))
  testthat::expect_error(gen_group_sizes(N, prev_K = c(.3,.1), rho_K))

})

testthat::test_that("generation of block probability matrix works", {

  N <- 1000
  K <- 3
  p_edge_within <- list(known1 = c(0.1,0.5), known2 = c(0.1,0.5), hidden = c(0.1,0.8))
  p_edge_between <- list(known1 = 0.1, known2 = 0.2, hidden = 0.3)

  type_names <-
    lapply(1:K, function(x) 0:1) %>%
    expand.grid() %>%
    apply(X = ., MARGIN = 1, FUN = paste0, collapse = "")

  block_mat <- gen_block_matrix(p_edge_within, p_edge_between)

  testthat::expect_equal(dim(block_mat), rep(length(type_names), times = 2))
  testthat::expect_equal(rownames(block_mat), type_names)
  testthat::expect_equal(colnames(block_mat), type_names)

  testthat::expect_error(gen_block_matrix(p_edge_within, p_edge_between = c(.1, .8)))
  testthat::expect_error(gen_block_matrix(p_edge_within = list(known1 = c(0.1,0.5),
                                                               known2 = c(0.1,0.5),
                                                               hidden = c(0.8)),
                                          p_edge_between))

})

