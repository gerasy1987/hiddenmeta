#' Sequential Sampling (SS) prevalence estimator by Gile (2011)
#'
#' @param data pass-through population data frame.
#' @param sampling_frame character string giving name of variable with sampling frame indicator for RDS sample. Defaults to "hidden" for the RDS samples from hidden population.
#' @param hidden_var character string specifying names of the hidden group variable name (associated probability of visibility should be named. \code{p_visible_[hidden_var]}). Defaults to "target" for the simulations
#' @param total_var name of variable containing size of population for prevalence estimation. Defaults to column named "total" in \code{data}.
#' @param n_coupons The number of recruitment coupons distributed to each enrolled subject (i.e. the maximum number of recruitees for any subject). By default it is taken by the attribute or data, else the maximum recorded number of coupons.
#' @param prefix character prefix used for RDS sample variables.
#' @param label character string describing the estimator.
#'
#' @return Data frame of SS prevalence estimates for a single study
#'
#' @references {Gile, Krista J. 2011 Improved Inference for Respondent-Driven Sampling Data with Application to HIV Prevalence Estimation, Journal of the American Statistical Association, 106, 135-146.}
#' @references {Gile, Krista J., Handcock, Mark S., 2010 Respondent-driven Sampling: An Assessment of Current Methodology, Sociological Methodology, 40, 285-327.}
#'
#' @export
#'
#' @import tidyselect
#' @importFrom dplyr mutate filter select group_by ungroup summarize pull arrange if_all
#' @importFrom RDS as.rds.data.frame RDS.SS.estimates
#' @importFrom purrr quietly
get_study_est_ss <-
  function(data,
           sampling_frame = "hidden",
           hidden_var = "target",
           total_var = "total",
           n_coupons = 3,
           prefix = "rds",
           label = "ss") {

    .quiet_rds_ss <- purrr::quietly(RDS::RDS.SS.estimates)

    .data_mod <-
      data |>
      dplyr::filter(dplyr::if_all(all_of(prefix), ~ . == 1))

    .total <- unique(.data_mod[[total_var]])

    if (length(.total) != 1) stop("Specified 'total_var' has more than one unique value in the sample.")

    .fit_rds_ss <-
      .data_mod |>
      dplyr::select(name,
                    all_of(c(hidden_var, paste0(sampling_frame, "_visible_out"))),
                    starts_with(prefix)) |>
      RDS::as.rds.data.frame(id = "name",
                             recruiter.id = paste0(prefix, "_from"),
                             network.size = paste0(sampling_frame, "_visible_out"),
                             time = paste0(prefix, "_t"),
                             max.coupons = n_coupons) |>
      .quiet_rds_ss(outcome.variable = hidden_var, N = .total)  |>
      (\(.) .$result$interval[c(1,5)])()

    return(
      ifelse(sampling_frame == "all",
             yes = paste0(hidden_var, "_prev"),
             no = paste0(hidden_var, "_prev_in_", sampling_frame)) |>
        (\(.)
         data.frame(estimator = paste0(., "_", label),
                    estimate = c(unname(.fit_rds_ss[1])),
                    se =   c(unname(.fit_rds_ss[2])),
                    inquiry = .)
        )()
    )
  }

#' Horvitz-Thompson prevalence estimator
#'
#' @param data pass-through population data frame
#' @param hidden_var variable containing hidden group membership indicator
#' @param weight_var variable containing sampling weights
#' @param total_var name of variable containing size of population for prevalence estimation
#' @param survey_design a formula describing the design of the survey (for bootstrap)
#' @param n_boot number of bootstrap resamples
#' @param parallel_boot logical, whether to compute bootstrap samples in parallel using \code{foreach} package
#' @param prefix character prefix used for sampling variables (has to include \code{[prefix]_weights})
#' @param label character string describing the estimator
#'
#' @return Data frame of HT estimates for a single study
#' @export
#'
#' @references {Salganik, Matthew J. "Variance estimation, design effects, and sample size calculations for respondent-driven sampling." Journal of Urban Health 83, no. 1 (2006): 98. \url{https://doi.org/10.1007/s11524-006-9106-x}.}
#'
#' @import tidyselect data.table
#' @importFrom dplyr mutate filter select group_by ungroup summarize pull arrange
#' @importFrom plyr laply
#' @importFrom doParallel registerDoParallel
get_study_est_ht <- function(data,
                             hidden_var = "hidden",
                             weight_var = "pps_weight",
                             total_var = "total",
                             survey_design = ~ pps_cluster + strata(pps_strata),
                             n_boot = 100,
                             parallel_boot = FALSE,
                             prefix = "pps",
                             label = "ht") {

  if (parallel_boot) {
    requireNamespace(c("doParallel", "parallel"))
    doParallel::registerDoParallel(cores = parallel::detectCores() - 1)
  }

  if (!data.table::is.data.table(data))
    data <- data.table::as.data.table(data)

  .data_mod <- data[get(prefix) == 1,]

  .total <- unique(.data_mod[[total_var]])

  if (length(.total) != 1) stop("Specified 'total_var' has more than one unique value in the sample.")

  if (any(is.na(.data_mod[[weight_var]])))
    stop("There are missing values in sampling weights provided")
  if (any(is.na(.data_mod[[hidden_var]])))
    stop("There are missing values in hidden population indicator provided")
  if (any(sort(unique(.data_mod[[hidden_var]])) != c(0,1)))
    stop("Hidden variable should be a binary indicator")

  .est_ht <-
    crossprod(.data_mod[[hidden_var]],
              .data_mod[[weight_var]])

  .est_ht_boot <-
    hiddenmeta:::get_rescaled_boot(data = .data_mod,
                                   survey_design = survey_design,
                                   n_boot = n_boot) |>
    plyr::laply(function(samp) {
      .data_mod  |>
        dplyr::mutate(rescale_wgt = samp$weight.scale)  |>
        dplyr::filter(rescale_wgt != 0)  |>
        (\(.) crossprod(.[[hidden_var]],
                        .$rescale_wgt * .[[weight_var]])
        )()
    }, .parallel = parallel_boot)

  return(
    data.frame(estimator = paste0(c("hidden_size", "hidden_prev"), "_", label),
               estimate = c(.est_ht, .est_ht/.total),
               se =  c(sd(.est_ht_boot), sd(.est_ht_boot/.total)),
               inquiry = c("hidden_size", "hidden_prev"))
  )
}

#' NSUM estimator
#'
#' @param data pass-through population data frame
#' @param hidden_var character vector containing names of hidden groups
#' @param known_vars character vector containing names of known groups
#' @param total integer giving total size of population
#' @param prefix character prefix used for PPS sample variables
#' @param label character string describing the estimator
#' @param weight_var character string giving name of the sampling weights variable
#' @param survey_design a formula describing the design of the survey
#' @param n_boot integer giving number of bootstrap re-samples
#' @param parallel_boot logical, whether to compute bootstrap samples in parallel using \code{foreach} package. Defaults to \code{FALSE}
#'
#' @return Data frame of NSUM estimates for a single study with PPS sample
#' @export
#'
#' @references Zheng, Tian, Matthew J. Salganik, and Andrew Gelman. "How many people do you know in prison? Using overdispersion in count data to estimate social structure in networks." Journal of the American Statistical Association 101.474 (2006): 409-423.
#' @references Laga, Ian, Le Bao, and Xiaoyue Niu. "A Correlated Network Scale-up Model: Finding the Connection Between Subpopulations." Journal of the American Statistical Association just-accepted (2023): 1-18.
#' @references Salganik, Matthew J. "Variance estimation, design effects, and sample size calculations for respondent-driven sampling." Journal of Urban Health 83, no. 1 (2006): 98. \url{https://doi.org/10.1007/s11524-006-9106-x}
#'
#' @import surveybootstrap
#' @importFrom magrittr `%>%` `%$%`
#' @importFrom networkscaleup killworth
#' @importFrom plyr llply
get_study_est_nsum <- function(
    data,
    hidden_var = "hidden_visible_out",
    known_vars = paste0(c("known", paste0("known_", 2:10)), "_visible_out"),
    total = 1000,
    prefix = "pps",
    label = "pps_nsum",
    weight_var = "pps_weight",
    survey_design = ~ pps_cluster + strata(pps_strata),
    n_boot = 1000,
    parallel_boot = FALSE) {

  if (parallel_boot) {
    requireNamespace(c("doParallel", "parallel"))
    doParallel::registerDoParallel(cores = parallel::detectCores() - 1)
  }
  .data_mod <- data[get(prefix) == 1,]

  if (!is.null(weight_var))
    .weights <- .data_mod[[weight_var]]
  else
    .weights <- 1

  .known_pops <-
    unlist(.data_mod[, lapply(.SD, unique),
                     .SDcols = paste0("total_", gsub("_visible_out",
                                                     replacement = "",
                                                     x = known_vars))])

  .data_est <-
    .data_mod[
      , c(known_vars, hidden_var), with = FALSE
    ][
      , lapply(.SD, function(x) x * .weights)
    ]

  .fit_nsum_est <-
    .data_est %>%
    networkscaleup::killworth(.,
                              known_sizes = .known_pops,
                              known_ind = 1:(ncol(.) - length(hidden_var)),
                              N = total,
                              model = "MLE") %>%
    .$sizes


  .fit_nsum_boot_se <-
    get_rescaled_boot(data = .data_mod,
                      survey_design = survey_design,
                      n_boot = n_boot) %>%
    plyr::llply(
      .data = .,
      .fun = function(wgt) {
        .data_est[, lapply(.SD, function(x) x * wgt[,2])] %>%
          networkscaleup::killworth(.,
                                    known_sizes = .known_pops,
                                    known_ind = 1:(ncol(.) - length(hidden_var)),
                                    N = total,
                                    model = "MLE") %>%
          .$sizes

      },
      .parallel = parallel_boot) %>%
    do.call("c", .)

  return(
    data.frame(
      estimator = paste0(c("hidden_size_"), label),
      estimate = .fit_nsum_est,
      se = sd(.fit_nsum_boot_se),
      inquiry = c("hidden_size"))
  )
}


#' Service/Object multiplier estimator
#'
#' @param data pass-through population data frame
#' @param service_var name of variable that indicates service/object use by respondent
#' @param total_service numeric value that indicates number of hidden population members who used the service. Defaults to truth from simulated dataset
#' @param seed_condition character string containing condition to define seeds. Defaults to "rds_from == -999" that applies to simulated RDS samples
#' @param n_boot number of bootstrap resamples
#' @param parallel_boot logical, whether to compute bootstrap samples in parallel using \code{foreach} package
#' @param prefix character prefix used for RDS sample variables
#' @param label character string describing the estimator
#'
#' @details Function currently requires variable "hidden_visible_out" to be present in the data supplied and represent the hidden population out-report
#'
#' @return Data frame of service/object multiplier population size estimates for single study
#' @export
#'
#' @references
#' Salganik, Matthew J. "Variance estimation, design effects, and sample size calculations for respondent-driven sampling." Journal of Urban Health 83, no. 1 (2006): 98. \url{https://doi.org/10.1007/s11524-006-9106-x}
#'
#' @import tidyselect data.table
#' @importFrom magrittr `%>%` `%$%`
#' @importFrom dplyr mutate filter select group_by ungroup summarize pull arrange bind_rows
#' @importFrom plyr laply
#' @importFrom doParallel registerDoParallel
get_study_est_multiplier <- function(data,
                                     service_var = "service_use",
                                     total_service = sum(data$service_use[data$hidden == 1]),
                                     seed_condition = "rds_from == -999",
                                     n_boot = 100,
                                     parallel_boot = FALSE,
                                     prefix = "rds",
                                     label = "multiplier") {

  if (parallel_boot) {
    requireNamespace(c("doParallel", "parallel"))
    doParallel::registerDoParallel(cores = parallel::detectCores() - 1)
  }

  if (!data.table::is.data.table(data))
    data <- data.table::as.data.table(data)

  .data_mod <- data[get(prefix) == 1,]

  .est_out <-
    total_service/mean(.data_mod[[service_var]])


  .est_boot <-
    get_rds_boot(data = .data_mod,
                 seed_condition = seed_condition,
                 in_coupon = paste0(prefix, "_own_coupon"),
                 out_coupon = paste0(prefix, "_coupon_"),
                 trait_var = "hidden_visible_out",
                 other_vars = c("hidden_visible_out", service_var, "name"),
                 n_boot = n_boot) %>%
    plyr::laply(., function(samp) {
      total_service/mean(samp[[service_var]])
    })

  return(
    data.frame(estimator = paste0("hidden_size_", label),
               estimate = .est_out,
               se =  sd(.est_boot),
               inquiry = "hidden_size")
  )
}

#' Generalized NSUM estimator
#'
#' @param data pass-through population data frame
#' @param label character string describing the estimator
#'
#' @keywords internal
#'
#' @return Data frame of HT estimates for single study
#'
get_study_est_gnsum <- function(data, label = "gnsum") {

  # res$sample.y.F.H <- with(frame.df,
  #                          sum((y.FH + y.notFH)*
  #                                sampling.weight))
  #
  # res$sample.dbar.F.U <- with(frame.df,
  #                             sum(d.degree*sampling.weight)/
  #                               sum(sampling.weight))
  #
  # res$sample.d.F.U <- with(frame.df,
  #                          sum(d.degree*sampling.weight))
  #
  # res$sample.dbar.F.F <- with(frame.df,
  #                             sum((d.FH + d.FnotH)*
  #                                   sampling.weight)/
  #                               sum(sampling.weight))
  #
  # res$sample.vbar.H.F <- with(hidden.df,
  #                             sum((v.FnotH + v.FH)*
  #                                   sampling.weight)/
  #                               sum(sampling.weight))
  #
  # res$sample.basic <- with(res,
  #                          (sample.y.F.H / sample.d.F.U)*
  #                            this.N)
  #
  # res$sample.adapted <- with(res,
  #                            sample.y.F.H / sample.dbar.F.F)
  #
  # res$sample.generalized <- with(res,
  #                                sample.y.F.H / sample.vbar.H.F)


  return(
    data.frame(estimator = paste0("hidden_prev_", label),
               estimate = NA,
               se =  NA,
               inquiry = "hidden_prev")
  )
}

#' Mark-recapture estimator for closed population
#'
#' @param data pass-through population data frame that contains capture indicators
#' @param capture_vars character vector giving names of variables with capture indicators
#' @param capture_parse character string giving expression to evaluation of which produces character vector giving names of variable with capture indicators. Defaults to \code{NULL}. This is useful when capture variables are stored in one column (e.g. if using TLS sampled locations for recapture indicators)
#' @param sample_condition character string with condition if the capture-recapture conducted on subsample of population (e.g. tls sample only)
#' @param model character string giving capture-recapture Log-Linear model to estimate
#' @param hidden_variable character string giving indicator of hidden population membership
#' @param label character string describing the estimator
#'
#' @return Data frame of Mark-recapture estimates for a single study
#' @export
#'
#' @references Louis-Paul Rivest, Sophie Baillargeon. “The Rcapture package.” (2019). \url{https://cran.r-project.org/package=Rcapture}.
#'
#' @import tidyselect data.table
#' @importFrom dplyr mutate filter select group_by ungroup summarize pull arrange bind_rows
#' @importFrom Rcapture closedp.bc periodhist
get_study_est_recapture <- function(
  data,
  capture_vars = NULL,
  capture_parse = NULL,
  sample_condition = NULL,
  model = "Mt",
  hidden_variable = "hidden",
  label = "recapture"
) {

  if (!data.table::is.data.table(data))
    data <- data.table::as.data.table(data)

  if (!is.null(sample_condition)) {
    data <- data[eval(parse(text = sample_condition)),]
  }

  if (!is.null(capture_parse)) {
    capture_vars <- data[, eval(parse(text = capture_parse))]
  }

  .est_out <-
    data[
      apply(sapply(capture_vars, function(x) get(x) == 1), 1, any),
    ][
      apply(sapply(hidden_variable, function(x) get(x) == 1), 1, any),
    ]

  if (nrow(.est_out) == 0) {
    warning("There were no hidden population member recaptures in the sample!")

    return(
      data.frame(estimator = paste0("hidden_size_", label),
                 estimate = NA_real_,
                 se =  NA_real_,
                 inquiry = "hidden_size")
    )

  } else {

    if (length(capture_vars) > 10) {

      .pool <- rep(length(capture_vars) %/% 5, 5)
      .pool[seq_len(length(capture_vars) %% 5)] <- (length(capture_vars) %/% 5) + 1

      .est_out <-
        .est_out  |>
        dplyr::select(all_of(capture_vars))  |>
        Rcapture::periodhist(vt = .pool) |>
        Rcapture::closedp.bc(dfreq = TRUE,
                             dtype = "hist",
                             # t = length(capture_vars),
                             m = model)
    } else {

      .est_out <-
        .est_out |>
        dplyr::select(all_of(capture_vars)) |>
        Rcapture::closedp.bc(dfreq = FALSE,
                             dtype = "hist",
                             # t = length(capture_vars),
                             m = model)

    }

    return(
      data.frame(estimator = paste0("hidden_size_", label),
                 estimate = .est_out$results[1],
                 se =  .est_out$results[2],
                 inquiry = "hidden_size")
    )

  }
}


#' Multiple Systems estimator of population size for Sparse Capture Data by Chan et al.
#'
#' @param data pass-through population data frame that contains capture indicators
#' @param capture_vars character vector giving names of variables with capture indicators
#' @param capture_parse character string giving expression to evaluation of which produces character vector giving names of variable with capture indicators. Defaults to \code{NULL}. This is useful when capture variables are stored in one column (e.g. if using TLS sampled locations for recapture indicators)
#' @param sample_condition character string with condition if the capture-recapture conducted on subsample of population (e.g. tls sample only)
#' @param method character string giving the estimation method to be passed to \code{estimatepopulation.0}. See \code{?SparseMSE::estimatepopulation.0} for more details. Defaults to \code{"stepwise"}
#' @param n_boot number of bootstrap resamples. Defaults to \code{NULL} which runs parametric estimation
#' @param hidden_variable character string giving indicator of hidden population membership
#' @param total_var name of variable containing size of population for prevalence estimation
#' @param label character string describing the estimator
#'
#' @return Data frame of Multiple System estimates for a single study
#'
#' @references Chan, Lax, Bernard W. Silverman, and Kyle Vincent. "Multiple systems estimation for sparse capture data: Inferential challenges when there are nonoverlapping lists." Journal of the American Statistical Association 116.535 (2021): 1297-1306.
#' @references Chan, Lax, Bernard W. Silverman, and Kyle Vincent. “The SparseMSE package.” (2022). \url{https://cran.r-project.org/package=SparseMSE}.
#'
#' @export
#'
#' @import data.table
#' @importFrom SparseMSE estimatepopulation.0
#' @importFrom Rcapture periodhist
get_study_est_mse <- function(
    data,
    capture_vars = NULL,
    capture_parse = NULL,
    sample_condition = NULL,
    method = "stepwise",
    n_boot = NULL,
    hidden_variable = "hidden",
    total_var = "total",
    label = "mse"
) {

  if (!data.table::is.data.table(data))
    data <- data.table::as.data.table(data)

  if (!is.null(sample_condition)) {
    data <- data[eval(parse(text = sample_condition)),]
  }

  .total <- unique(data[[total_var]])

  if (length(.total) != 1) stop("Specified 'total_var' has more than one unique value in the sample.")

  if (!is.null(capture_parse)) {
    capture_vars <- data[, eval(parse(text = capture_parse))]
  }

  .est_out <-
    data[
      apply(sapply(capture_vars, function(x) get(x) == 1), 1, any),
    ][
      apply(sapply(hidden_variable, function(x) get(x) == 1), 1, any),
    ]



  if (nrow(.est_out) == 0) {
    warning("There were no hidden population member recaptures in the sample!")

    return(
      data.frame(estimator = paste0("hidden_size_", label),
                 estimate = NA_real_,
                 se =  NA_real_,
                 inquiry = "hidden_size")
    )

  } else {

    if (length(capture_vars) > 10) {
      .pool <- rep(length(capture_vars) %/% 5, 5)
      .pool[seq_len(length(capture_vars) %% 5)] <- (length(capture_vars) %/% 5) + 1

      if (!is.null(n_boot)) {
        .est_out <-
          .est_out |>
          dplyr::select(all_of(capture_vars)) |>
          Rcapture::periodhist(vt = .pool) |>
          SparseMSE::estimatepopulation(method = method, nboot = n_boot)
      } else {
        .est_out <-
          .est_out |>
          dplyr::select(all_of(capture_vars)) |>
          Rcapture::periodhist(vt = .pool) |>
          SparseMSE::estimatepopulation.0(method = method)
      }

    } else if (length(capture_vars) <= 10 & !is.null(n_boot)) {
      if (!is.null(n_boot)) {
        .est_out <-
          SparseMSE::estimatepopulation(.est_out[, .(n = .N), by = capture_vars],
                                        method = method, nboot = n_boot)
      } else {
        .est_out <-
          SparseMSE::estimatepopulation(.est_out[, .(n = .N), by = capture_vars],
                                        method = method)
      }
    }


    if (!is.null(n_boot)) {
      return(
        data.frame(estimator = paste0(c("hidden_size", "hidden_prev"), "_", label),
                   estimate = c(unname(.est_out$popest["point est."]),
                                unname(.est_out$popest["point est."])/.total),
                   se =  c(sd(.est_out$bootreps), sd(.est_out$bootreps/.total)),
                   inquiry = c("hidden_size", "hidden_prev"))
      )
    } else {
      return(
        data.frame(estimator = paste0("hidden_size_", label),
                   estimate = unname(.est_out$estimate["point est."]),
                   # calculate SE as SD of log-normal distribution
                   # with logmean and logsd from MSEfit object
                   se =  unname(
                     sqrt((exp(summary(.est_out$MSEfit$fit)$cov.unscaled[1, 1]) - 1) *
                            exp(2*.est_out$MSEfit$fit$coefficients[1] +
                                  summary(.est_out$MSEfit$fit)$cov.unscaled[1, 1]))),
                   inquiry = "hidden_size")
      )
    }
  }
}