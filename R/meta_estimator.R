#' Use Stan model to estimate bias and estimands
#'
#' @param data pass-through meta population or meta sample data frame
#' @param sampling_variable name of variable storing meta analysis sampling information
#' @param which_estimand name of study level estimand for meta analysis
#' @param benchmark named list of length 2 giving benchmark sampling-estimator pair (only accepts one value across studies for now)
#'
#' @return
#' @export
#'
#' @import dplyr
#' @import rstan
get_meta_estimates <- function(
  data,
  sampling_variable = "meta",
  which_estimand = "hidden_size",
  benchmark = list(sample = "pps", estimator = "ht")) {

  .stan_data <-
    data %>%
    dplyr::filter(across(all_of(sampling_variable), ~ . == 1),
                  inquiry %in% which_estimand)

  .samp_ests <- unique(.stan_data[,c("sample", "estimator")])
  .benchmark_sample_est <-
    which(lapply(.samp_ests$sample, paste0, collapse = "_") ==
            paste0(benchmark$sample, collapse = "_") &
            .samp_ests$estimator == benchmark$estimator)
  .samp_ests <-
    rbind(.samp_ests[.benchmark_sample_est,],
          .samp_ests[-.benchmark_sample_est,])
  .samp_est_names <-
    paste0(sapply(.samp_ests$sample, paste0, collapse = "_"), "_", .samp_ests$estimator)
  .studies <- unique(.stan_data$study)

  .N <- length(.studies)
  .K <- nrow(.samp_ests)

  # get ids of unique samp-est pairs in observed data
  .samp_est_ids <-
    mapply(
      x = .samp_ests$sample, y = .samp_ests$estimator,
      function(x, y) {
        lapply(.stan_data$sample, paste0, collapse = "_") ==
          paste0(x, collapse = "_") & .stan_data$estimator == y
      },
      SIMPLIFY = FALSE
    )

  .stan_data <-
    .samp_est_ids %>%
    {
      c(lapply(X = ., FUN = sum),
        lapply(X = ., FUN = function(x) which(.studies %in% .stan_data$study[x])),
        lapply(X = ., FUN = function(x) .stan_data$estimate[x]),
        lapply(X = ., FUN = function(x) .stan_data$se[x]))
    } %>%
    `names<-`(
      .,
      apply(X = expand.grid(1:.K, c("n", "observed", "est", "est_se")),
            MARGIN = 1,
            FUN = function(x) paste0(x[c(2,1)], collapse = ""))
    ) %>%
    c(N = .N, K = .K, .)

  .fit <-
    suppressMessages(
      rstan::stan(model_code = get_meta_stan(.stan_data),
                  data = .stan_data,
                  iter = 2000, chains = 8, seed = 19871223, cores = parallel::detectCores() - 1,
                  verbose = FALSE) %>%
        rstan::extract(.)
    )

  .biases <-
    .fit$rel_bias %>%
    cbind(1, .) %>%
    plyr::alply(
      .margins = 1,
      .fun = function(x) {
        c(matrix(x, nrow = .K, ncol = .K, byrow = TRUE) /
            matrix(x, nrow = .K, ncol = .K, byrow = FALSE))
      }) %>%
    do.call(rbind, .) %>%
    plyr::aaply(.margins = 2,
                .fun = function(x) c(est = mean(x),
                                     se = sd(x))) %>%
    plyr::alply(.margins = 2,
                function(x) {
                  matrix(data = x,
                         nrow = .K, ncol = .K,
                         byrow = FALSE,
                         dimnames = list(.samp_est_names,
                                         .samp_est_names))
                })


  .study_ests <-
    .fit$alpha %>%
    plyr::aaply(.margins = 2,
                .fun = function(x) c(est = mean(x),
                                     se = sd(x)))

  return(
    data.frame(estimator_label = c(paste0("rel_bias_", .samp_est_names, "_meta"),
                                   paste0(.studies, "_meta")),
               estimate = c(unname(.biases[[1]][1,]), unname(.study_ests[,1])),
               se =   c(unname(.biases[[2]][1,]), unname(.study_ests[,2])),
               inquiry_label = c(paste0("rel_bias_", .samp_est_names),
                                 paste0(.studies, "_", which_estimand))
    )
  )

}
