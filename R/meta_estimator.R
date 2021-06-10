#' Use Stan model to estimate bias and estimands
#'
#' @param data pass-through meta population or meta sample data frame
#' @param sampling_variable name of variable storing meta analysis sampling information
#' @param which_estimand name of study level estimand for meta analysis
#'
#' @return
#' @export
#'
#' @import dplyr
#' @importFrom purrr map2_int
get_meta_estimators <-
  function(data,
           sampling_variable = "meta",
           which_estimand = "hidden_size") {

    stan_data <-
      data %>%
      dplyr::filter(across(all_of(sampling_variable), ~ . == 1),
                    inquiry %in% which_estimand)

    samp_ests <- unique(stan_data[,c("sample", "estimator")])
    studies <- unique(stan_data$study)


    N <- length(studies)
    K <- nrow(samp_ests)

    # get ids of unique samp-est pairs in observed data
    samp_est_ids <-
      mapply(
        x = samp_ests$sample, y = samp_ests$estimator,
        function(x,y) {
          lapply(stan_data$sample, paste0, collapse = "_") ==
            paste0(x, collapse = "_") & stan_data$estimator == y
        },
        SIMPLIFY = FALSE
      )

    stan_data <-
      samp_est_ids %>%
      {
        c(lapply(X = ., FUN = sum),
          lapply(X = ., FUN = function(x) which(studies %in% stan_data$study[x])),
          lapply(X = ., FUN = function(x) stan_data$estimate[x]),
          lapply(X = ., FUN = function(x) stan_data$se[x]))
      } %>%
      `names<-`(
        .,
        apply(X = expand.grid(1:K, c("n", "observed", "est", "est_se")),
              MARGIN = 1,
              FUN = function(x) paste0(x[c(2,1)], collapse = ""))
      ) %>%
      c(N = N, K = K, .)

  # parse_model <- stanc(model_code = get_meta_stan(K), model_name = "meta_stan")

  fit <-
    rstan::stan(
      # file = "R/meta.stan",
      model_code = get_meta_stan(stan_data),
      data = stan_data,
      iter = 1000, chains = 4, seed = 19871223, cores = parallel::detectCores() - 1,
      verbose = TRUE)

  # data.frame(estimator_label = c(paste0("prev_", 1:N)),
  #            estimate = c(apply(fit$alpha, 2, mean)),
  #            sd =   c(apply(fit$alpha, 2, sd)),
  #            inquiry_label = c(paste0("hidden_prev", 1:N)),
  #            big_Rhat = big_Rhat
  # )

}
