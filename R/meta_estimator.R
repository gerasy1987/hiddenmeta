#' Use Stan model to estimate bias and estimands
#'
#' @param data pass-through meta population or meta sample data frame
#' @param sample_label name of variable storing meta analysis sampling information
#' @param benchmark named list of length 2 giving benchmark sampling-estimator pair (only accepts one value across studies for now)
#' @param hidden_prior list of two hyperpriors, on means and standard errors of each included sampling-estimator pairs. Names of list objects should be "mean" and "se". If one number provided for a hyperprior it gets expanded to all sampling-estimator pairs
#' @param rel_bias_prior list of two hyperpriors, on means and standard errors of relative bias. Names of list objects should be "mean" and "se".
#' @param stan_model function that takes \code{stan_data} as input and produces compilable stan model object. Alternatively this can me character string that specifies name of model pre-compiled by the package, either \code{"stan_prev"} or \code{"stan_size"}
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
    stan_model = "stan_prev",
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


  if (is.function(stan_model)) {

    .stan_compile <-
      rstan::stan_model(model_code = stan_model(.stan_data))

    .fit <-
      do.call(rstan::sampling, c(list(object = .stan_compile,
                                      data = .stan_data,
                                      verbose = TRUE,
                                      refresh = 0),
                                 control_params)) |>
      rstan::extract()

  } else if (is.character(stan_model)) {

    .fit <-
      do.call(rstan::sampling, c(list(object = stanmodels[[stan_model]],
                                      data = .stan_data,
                                      verbose = TRUE,
                                      refresh = 0),
                                 control_params))

  }

  .fit_extract <- rstan::extract(.fit)

  .biases <-
    .fit_extract$bias |>
    plyr::aaply(.margins = 2,
                .fun = function(x) c(est = mean(x),
                                     se = sd(x))) |>
    (\(.) .[.designs$benchmark != 1,])()

  .fit_extract$bias <-
    .fit_extract$bias[,.designs$benchmark != 1] |>
    `colnames<-`(.designs$study_design[.designs$benchmark != 1])

  .fit_extract$alpha <-
    .fit_extract$alpha |>
    `colnames<-`(.study_order)

  .study_ests <-
    plyr::aaply(.data = .fit_extract$alpha,
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
    return( list(estimates = .out, fit = .fit) )
  }

}
