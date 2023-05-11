#' Use Stan model to estimate bias and estimands
#'
#' @param data pass-through meta population or meta sample data frame
#' @param sample_label name of variable storing meta analysis sampling information
#' @param which_estimand name of study level estimand for meta analysis
#' @param benchmark named list of length 2 giving benchmark sampling-estimator pair (only accepts one value across studies for now)
#' @param stan_handler function that takes stan_data as input and produces compilable stan model object
#' @param hidden_prior list of two hyperpriors, on means and standard errors of each included sampling-estimator pairs. Names of list objects should be "mean" and "se". If one number provided for a hyperprior it gets expanded to all sampling-estimator pairs
#' @param rel_bias_prior list of two hyperpriors, on means and standard errors of relative bias. Names of list objects should be "mean" and "se".
#' @param control_params list of additional parameters to pass to Stan fit function. These can include number of iterations, chains, thinning, seed and number of cores to use
#'
#' @return Data frame of meta level estimates and pertaining estimand names
#'
#'
#' @import tidyselect
#' @import rstan
#' @importFrom magrittr `%>%` `%$%`
#' @importFrom dplyr mutate filter select group_by ungroup summarize pull arrange rename_with left_join bind_rows if_all if_any
#' @importFrom parallel detectCores
#' @importFrom plyr alply
get_meta_estimates_old <- function(
  data,
  sample_label = "meta",
  which_estimand = "hidden_size",
  benchmark = list(sample = "pps", estimator = "ht"),
  stan_handler = get_meta_stan3,
  hidden_prior = NULL,
  rel_bias_prior = list(mean = 1, se = 10),
  control_params = list(
    iter = 8000, chains = 8, thin = 10,
    seed = 872312,
    cores = 1)
) {

  # allow for NULL sample label to keep all data
  if (is.null(sample_label)) {
    .stan_data <- data
  } else if (is.character(sample_label)) {
    .stan_data <-
      data %>%
      dplyr::filter(dplyr::if_any(sample_label, ~ . == 1),
                    inquiry %in% which_estimand)
  }

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

  if (!is.null(hidden_prior) & is.list(hidden_prior)) {
    .alpha_prior <- do.call(cbind, hidden_prior)

    if (nrow(.alpha_prior) == 1) {
      .alpha_prior <- .alpha_prior[rep(1, times = .N),]
    } else if (nrow(.alpha_prior) == .N) {
      .alpha_prior <- .alpha_prior[.studies,]
    } else {
      stop("There is mismatch in length of hyperpriors on mean and std. error of target parameters. Length of each element should be equal either to length of number of studies included or to 1")
    }
  } else if (is.null(hidden_prior)) {
    .alpha_prior <-
      .stan_data %>%
      dplyr::group_by(study) %>%
      dplyr::summarize(mean = sum(estimate / (se * sum(1/se))), se = 2 * mean(se / (se * sum(1/se)))) %>%
      { .[order(match(.$study,.studies)),] } %>%
      dplyr::select(-study) %>%
      as.matrix()
  } else {
    stop("Hyperpriors on target parameters are provided in wrong format (should be NULL or list of two numeric objects)")
  }

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
            FUN = function(x) paste0(x[2], as.integer(x[1])))
    ) %>%
    c(list(N = .N, K = .K,
           alpha_mean = unname(.alpha_prior[,"mean"]), alpha_se = unname(.alpha_prior[,"se"]),
           rel_bias_mean = rel_bias_prior$mean, rel_bias_se = rel_bias_prior$se),
      .)

  .stan_model <-
    rstan::stan_model(model_code = stan_handler(.stan_data))

  .fit <-
    do.call(rstan::sampling, c(list(object = .stan_model,
                                  data = .stan_data,
                                  verbose = FALSE,
                                  refresh = 0),
                               control_params)) %>%
    rstan::extract(.)

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
    data.frame(estimator = c(paste0("rel_bias_", .samp_est_names, "_meta"),
                                   paste0(.studies, "_meta")),
               estimate = c(unname(.biases[[1]][1,]), unname(.study_ests[,1])),
               se =   c(unname(.biases[[2]][1,]), unname(.study_ests[,2])),
               inquiry = c(paste0("rel_bias_", .samp_est_names),
                           paste0(.studies, "_", which_estimand)),
               # note assumes first item is benchmark
               prior = c(1, rep(rel_bias_prior[[1]], .K-1), .alpha_prior[,1]) %>% unlist,
               prior_sd = c(0, rep(rel_bias_prior[[2]], .K-1), .alpha_prior[,2]) %>% unlist,
               stringsAsFactors = FALSE
    )
  )

}


#' Use Stan model to estimate bias and estimands
#'
#' @param data pass-through meta population or meta sample data frame
#' @param sample_label name of variable storing meta analysis sampling information
#' @param benchmark named list of length 2 giving benchmark sampling-estimator pair (only accepts one value across studies for now)
#' @param hidden_prior list of two hyperpriors, on means and standard errors of each included sampling-estimator pairs. Names of list objects should be "mean" and "se". If one number provided for a hyperprior it gets expanded to all sampling-estimator pairs
#' @param rel_bias_prior list of two hyperpriors, on means and standard errors of relative bias. Names of list objects should be "mean" and "se".
#' @param stan_handler function that takes stan_data as input and produces compilable stan model object
#' @param control_params list of additional parameters to pass to Stan fit function. These can include number of iterations, chains, thinning, seed and number of cores to use
#' @param return_stanfit logical value identifying whether the Stan model fit object should be returned in addition to estimates data frame
#'
#' @return Data frame of meta level estimates and pertaining estimand names
#'
#' @export
#'
#' @import tidyselect
#' @import rstan
#' @importFrom dplyr mutate filter select group_by ungroup summarize pull arrange rename_with left_join bind_rows if_all if_any
#' @importFrom parallel detectCores
#' @importFrom plyr alply aaply
get_meta_estimates <- function(
    data,
    sample_label = "meta",
    benchmark = list(sample = "pps", estimator = "ht"),
    hidden_prior = list(mean = .3, sd = .2),
    bias_prior = list(mean = 1, sd = 3),
    stan_handler = get_meta_stan4,
    control_params = list(
      iter = 8000, chains = 8, thin = 10,
      seed = 872312,
      cores = 1),
    return_stanfit = FALSE
) {

  # allow for NULL sample label to keep all data
  if (is.null(sample_label)) {
    .stan_data <- data
  } else if (is.character(sample_label)) {
    .stan_data <-
      data |> dplyr::filter(dplyr::if_any(dplyr::all_of(sample_label), ~ . == 1))
  }

  # make order of studies
  if (is.null(names(hidden_prior$mean)) & length(hidden_prior$mean) == 1) {
    .study_order <- unique(paste0(.stan_data$study, "_", .stan_data$inquiry))
  } else if (is.null(names(hidden_prior$mean)) & length(hidden_prior$mean) != 1) {
    warning("Study level estimand names are inferred from the sample data")

    .study_order <- unique(paste0(.stan_data$study, "_", .stan_data$inquiry))
  } else {
    .study_order <- names(hidden_prior$mean)
  }


  .stan_data <-
    .stan_data |>
    dplyr::mutate(
      study_inquiry = paste0(study, "_", inquiry),
      study_inquiry_id =
        as.integer(
          plyr::mapvalues(study_inquiry,
                          from = .study_order,
                          to = seq_along(.study_order))),
      study_design =
        paste0(inquiry, "_",
               sapply(sample, paste0, collapse = "_"), "_",
               estimator),
      study_design_id = as.integer(as.factor(study_design))
    )

  # prepare design mapping and bias priors
  .designs <-
    unique(.stan_data[,c("inquiry", "sample", "estimator",
                         "study_design", "study_design_id")]) |>
    dplyr::mutate(
      benchmark =
        as.integer((lapply(sample, paste0, collapse = "_") %in% benchmark$sample) &
        (estimator == benchmark$estimator)),
      bias_prior_mean = bias_prior$mean,
      bias_prior_sd = dplyr::if_else(benchmark == 1, 0.001, bias_prior$sd)
    ) |>
    dplyr::arrange(study_design_id)

  # dims of priors
  .N <- length(.study_order)
  .K <- nrow(.designs)

  # prepare priors for study level inquiries
  # if no priors are provided use provided estimates
  if (!is.null(hidden_prior) & is.list(hidden_prior)) {
    .alpha_prior <-
      do.call(cbind, hidden_prior) |>
      `rownames<-`(NULL)

    if (nrow(.alpha_prior) == 1)
      .alpha_prior <- .alpha_prior[rep(1, times = .N),]
  } else if (is.null(hidden_prior)) {
    .alpha_prior <-
      .stan_data |>
      dplyr::group_by(study_inquiry) |>
      dplyr::summarize(mean = sum(estimate / (se * sum(1/se))),
                       sd = 2 * mean(se / (se * sum(1/se)))) |>
      dplyr::arrange(match(study_inquiry, .study_order)) |>
      dplyr::select(-study_inquiry) |>
      as.matrix()
  } else {
    stop("Priors on study level inquiries are provided in wrong format. Should be NULL or list of two numeric objects of equal length giving mean and SD of priors.")
  }

  .stan_data <-
    list(
      N = .N,
      K = .K,
      n_ests = nrow(.stan_data),
      alpha_prior_mean = .alpha_prior[,"mean"],
      alpha_prior_sd = .alpha_prior[,"sd"],
      bias_prior_mean = .designs$bias_prior_mean,
      bias_prior_sd = .designs$bias_prior_sd,
      study = .stan_data$study_inquiry_id,
      design = .stan_data$study_design_id,
      ests = .stan_data$estimate,
      ses = .stan_data$se
    )


  .stan_model <-
    rstan::stan_model(model_code = stan_handler(.stan_data))

  .fit <-
    do.call(rstan::sampling, c(list(object = .stan_model,
                                    data = .stan_data,
                                    verbose = TRUE,
                                    refresh = 0),
                               control_params)) |>
    rstan::extract()

  .biases <-
    .fit$bias |>
    plyr::aaply(.margins = 2,
                .fun = function(x) c(est = mean(x),
                                     se = sd(x))) |>
    (\(.) .[.designs$benchmark != 1,])()

  .study_ests <-
    plyr::aaply(.data = .fit$alpha,
                .margins = 2,
                .fun = function(x) c(est = mean(x),
                                     se = sd(x)))

  .out <-
    data.frame(
      estimator =
        c(paste0("rel_bias_", .designs$study_design[.designs$benchmark != 1], "_meta"),
          paste0(.study_order, "_meta")),
      estimate = unname(c(.biases[,"est"],.study_ests[,"est"])),
      se =   unname(c(.biases[,"se"], .study_ests[,"se"])),
      inquiry = c(paste0("rel_bias_",
                         .designs$study_design[.designs$benchmark != 1]),
                  paste0(.study_order)),
      # note assumes first item is benchmark
      prior = c(.designs$bias_prior_mean[.designs$benchmark != 1],
                .alpha_prior[,"mean"]),
      prior_sd = c(.designs$bias_prior_sd[.designs$benchmark != 1],
                   .alpha_prior[,"sd"]),
      stringsAsFactors = FALSE
    )

  if (!return_stanfit) {
    return(.out)
  } else {
    return(
      list(estimates = .out,
           fit = .fit)
    )
  }

}
