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
  function(
    data,
    sampling_frame = "hidden",
    hidden_var = "target",
    total_var = "total",
    n_coupons = 3,
    prefix = "rds",
    label = "ss"
  ) {
    .quiet_rds_ss <- purrr::quietly(RDS::RDS.SS.estimates)

    .data_mod <-
      data |>
      dplyr::filter(dplyr::if_all(all_of(prefix), ~ . == 1))

    .total <- unique(.data_mod[[total_var]])

    if (length(.total) != 1) {
      stop(
        "Specified 'total_var' has more than one unique value in the sample."
      )
    }

    .fit_rds_ss <-
      .data_mod |>
      dplyr::select(
        name,
        all_of(c(hidden_var, paste0(sampling_frame, "_visible_out"))),
        starts_with(prefix)
      ) |>
      RDS::as.rds.data.frame(
        id = "name",
        recruiter.id = paste0(prefix, "_from"),
        network.size = paste0(sampling_frame, "_visible_out"),
        time = paste0(prefix, "_t"),
        max.coupons = n_coupons
      ) |>
      .quiet_rds_ss(outcome.variable = hidden_var, N = .total) |>
      (\(.) .$result$interval[c(1, 5)])()

    return(
      ifelse(
        sampling_frame == "all",
        yes = paste0(hidden_var, "_prev"),
        no = paste0(hidden_var, "_prev_in_", sampling_frame)
      ) |>
        (\(.) {
          data.frame(
            estimator = paste0(., "_", label),
            estimate = c(unname(.fit_rds_ss[1])),
            se = c(unname(.fit_rds_ss[2])),
            inquiry = .
          )
        })()
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
get_study_est_ht <- function(
  data,
  hidden_var = "hidden",
  weight_var = "pps_weight",
  total_var = "total",
  survey_design = ~ pps_cluster + strata(pps_strata),
  n_boot = 100,
  parallel_boot = FALSE,
  prefix = "pps",
  label = "ht"
) {
  if (parallel_boot) {
    requireNamespace(c("doParallel", "parallel"))
    doParallel::registerDoParallel(cores = parallel::detectCores() - 1)
  }

  if (!data.table::is.data.table(data)) {
    data <- data.table::as.data.table(data)
  }

  .data_mod <- data[get(prefix) == 1, ]

  .total <- unique(.data_mod[[total_var]])

  if (length(.total) != 1) {
    stop("Specified 'total_var' has more than one unique value in the sample.")
  }

  if (any(is.na(.data_mod[[weight_var]]))) {
    stop("There are missing values in sampling weights provided")
  }
  if (any(is.na(.data_mod[[hidden_var]]))) {
    stop("There are missing values in hidden population indicator provided")
  }
  if (any(sort(unique(.data_mod[[hidden_var]])) != c(0, 1))) {
    stop("Hidden variable should be a binary indicator")
  }

  .est_ht <-
    crossprod(.data_mod[[hidden_var]], .data_mod[[weight_var]])

  .est_ht_boot <-
    hiddenmeta:::get_rescaled_boot(
      data = .data_mod,
      survey_design = survey_design,
      n_boot = n_boot
    ) |>
    plyr::laply(
      function(samp) {
        .data_mod |>
          dplyr::mutate(rescale_wgt = samp$weight.scale) |>
          dplyr::filter(rescale_wgt != 0) |>
          (\(.) crossprod(.[[hidden_var]], .$rescale_wgt * .[[weight_var]]))()
      },
      .parallel = parallel_boot
    )

  return(
    data.frame(
      estimator = paste0(c("hidden_size", "hidden_prev"), "_", label),
      estimate = c(.est_ht, .est_ht / .total),
      se = c(sd(.est_ht_boot), sd(.est_ht_boot / .total)),
      inquiry = c("hidden_size", "hidden_prev")
    )
  )
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
  if (!data.table::is.data.table(data)) {
    data <- data.table::as.data.table(data)
  }

  if (!is.null(sample_condition)) {
    data <- data[eval(parse(text = sample_condition)), ]
  }

  .total <- unique(data[[total_var]])

  if (length(.total) != 1) {
    stop("Specified 'total_var' has more than one unique value in the sample.")
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
      data.frame(
        estimator = paste0("hidden_size_", label),
        estimate = NA_real_,
        se = NA_real_,
        inquiry = "hidden_size"
      )
    )
  } else {
    if (length(capture_vars) > 10) {
      .pool <- rep(length(capture_vars) %/% 5, 5)
      .pool[seq_len(length(capture_vars) %% 5)] <- (length(capture_vars) %/%
        5) +
        1

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
          SparseMSE::estimatepopulation(
            .est_out[, .(n = .N), by = capture_vars],
            method = method,
            nboot = n_boot
          )
      } else {
        .est_out <-
          SparseMSE::estimatepopulation(
            .est_out[, .(n = .N), by = capture_vars],
            method = method
          )
      }
    }

    if (!is.null(n_boot)) {
      return(
        data.frame(
          estimator = paste0(c("hidden_size", "hidden_prev"), "_", label),
          estimate = c(
            unname(.est_out$popest["point est."]),
            unname(.est_out$popest["point est."]) / .total
          ),
          se = c(sd(.est_out$bootreps), sd(.est_out$bootreps / .total)),
          inquiry = c("hidden_size", "hidden_prev")
        )
      )
    } else {
      return(
        data.frame(
          estimator = paste0("hidden_size_", label),
          estimate = unname(.est_out$estimate["point est."]),
          # calculate SE as SD of log-normal distribution
          # with logmean and logsd from MSEfit object
          se = unname(
            sqrt(
              (exp(summary(.est_out$MSEfit$fit)$cov.unscaled[1, 1]) - 1) *
                exp(
                  2 *
                    .est_out$MSEfit$fit$coefficients[1] +
                    summary(.est_out$MSEfit$fit)$cov.unscaled[1, 1]
                )
            )
          ),
          inquiry = "hidden_size"
        )
      )
    }
  }
}
