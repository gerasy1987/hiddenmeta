#' SS-PSE population size estimator by Handcock, Gile and Mar
#'
#' @param data pass-through population data frame
#' @param prior_mean the mean of the prior distribution on the population size for SS-PSE estimation
#' @param n_coupons The number of recruitment coupons distributed to each enrolled subject (i.e. the maximum number of recruitees for any subject). By default it is taken by the attribute or data, else the maximum recorded number of coupons.
#' @param mcmc_params named list of parameters passed to \code{sspse::posteriorsize} for MCMC sampling,
#' @param additional_params named list of additional parameter passed to \code{sspse::posteriorsize} . If empty \code{sspse::posteriorsize} uses default parameters.
#' @param total integer giving the total size of population
#' @param prefix character prefix used for RDS sample variables
#' @param label character string describing the estimator
#'
#' @return Data frame of SS-PSE estimates for a single study
#'
#' @references {Handcock, Mark S., Krista J. Gile, and Corinne M. Mar. “Estimating Hidden Population Size Using Respondent-Driven Sampling Data.” Electronic Journal of Statistics 8, no. 1 (2014): 1491–1521. \url{https://doi.org/10.1214/14-EJS923}.}
#'
#' @export
#'
#' @import tidyselect
#' @importFrom magrittr `%>%` `%$%`
#' @importFrom dplyr mutate filter select group_by ungroup summarize pull arrange if_all
#' @importFrom sspse posteriorsize
#' @importFrom purrr quietly
#' @importFrom stats quantile
get_study_est_sspse <- function(
  data,
  prior_mean = .1 * nrow(data),
  n_coupons = 3,
  total = 2000,
  mcmc_params = list(interval = 5, warmup = 2000, samplesize = 500),
  additional_params = list(),
  prefix = "rds",
  label = "sspse"
) {
  .quiet_sspse <- purrr::quietly(sspse::posteriorsize)

  network_sizes <- data %>%
    dplyr::filter(dplyr::if_all(all_of(prefix), ~ . == 1)) %>%
    dplyr::pull("hidden_visible_out")

  # Check 1: Empty sample
  if (length(network_sizes) == 0) {
    warning("No RDS sample found. Returning NA.")
    return(
      data.frame(
        estimator = paste0("hidden_size_", label),
        estimate = NA_real_,
        se = NA_real_,
        inquiry = c("hidden_size")
      )
    )
  }

  # Check 2: Too sparse (>90% zeros or max <= 1)
  prop_zero <- mean(network_sizes == 0)
  if (prop_zero > 0.90 || max(network_sizes) <= 1) {
    warning(paste0(
      "Network sizes too sparse (",
      round(prop_zero * 100, 1),
      "% zeros, max=",
      max(network_sizes),
      "). SS-PSE estimate unreliable. Returning NA."
    ))
    return(
      data.frame(
        estimator = paste0("hidden_size_", label),
        estimate = NA_real_,
        se = NA_real_,
        inquiry = c("hidden_size")
      )
    )
  }

  # Check 3: Verify we have variation
  if (length(unique(network_sizes)) < 3) {
    warning("Insufficient variation in network sizes. Returning NA.")
    return(
      data.frame(
        estimator = paste0("hidden_size_", label),
        estimate = NA_real_,
        se = NA_real_,
        inquiry = c("hidden_size")
      )
    )
  }

  # Try estimation with error handling
  .fit_sspse <- tryCatch(
    {
      do.call(
        .quiet_sspse,
        c(
          list(
            s = network_sizes,
            interval = mcmc_params$interval,
            samplesize = mcmc_params$samplesize,
            warmup = mcmc_params$warmup,
            mean.prior.size = prior_mean,
            verbose = FALSE,
            K = round(stats::quantile(network_sizes, 0.80)),
            max.coupons = n_coupons
          ),
          additional_params
        )
      )
    },
    error = function(e) {
      warning("SS-PSE estimation failed: ", e$message)
      return(NULL)
    }
  )

  # Check if estimation failed
  if (is.null(.fit_sspse) || is.null(.fit_sspse$result)) {
    return(
      data.frame(
        estimator = paste0("hidden_size_", label),
        estimate = NA_real_,
        se = NA_real_,
        inquiry = c("hidden_size")
      )
    )
  }

  return(
    data.frame(
      estimator = paste0("hidden_size_", label),
      estimate = c(unname(.fit_sspse$result$N["Median AP"])),
      se = c(sd(.fit_sspse$result$sample[, "N"])),
      inquiry = c("hidden_size")
    )
  )
}

#' Chords population size estimator by Berchenko, Rosenblatt and Frost
#'
#' @param data pass-through population data frame
#' @param type a character vector with the type of estimation. Can be one of \code{mle}, \code{integrated}, or \code{jeffreys}. See \code{?chords::Estimate.b.k} and the original paper from the references for details
#' @param seed_condition character string containing condition to define seeds. Defaults to "rds_from == -999" that applies to simulated RDS samples
#' @param n_boot number of bootstrap resamples
#' @param parallel_boot logical, whether to compute bootstrap samples in parallel using \code{foreach} package
#' @param prefix character string prefix used for RDS sample variables
#' @param label character string describing the estimator
#'
#' @return Data frame of Chords estimates for a single study with RDS sample
#'
#' @references {Berchenko, Yakir, Jonathan D. Rosenblatt, and Simon D. W. Frost. “Modeling and Analyzing Respondent-Driven Sampling as a Counting Process.” Biometrics 73, no. 4 (2017): 1189–98. \url{https://doi.org/10.1111/biom.12678}.}
#' @references {Salganik, Matthew J. "Variance estimation, design effects, and sample size calculations for respondent-driven sampling." Journal of Urban Health 83, no. 1 (2006): 98. \url{https://doi.org/10.1007/s11524-006-9106-x}.}
#'
#' @export
#'
#' @import tidyselect
#' @importFrom magrittr `%>%` `%$%`
#' @importFrom dplyr mutate filter select group_by ungroup summarize pull arrange rename_with left_join bind_rows if_all
#' @importFrom chords initializeRdsObject Estimate.b.k makeJackControl
#' @importFrom purrr quietly
#' @importFrom doParallel registerDoParallel
#' @importFrom stats weighted.mean
get_study_est_chords <- function(
  data,
  type = c("mle", "integrated", "jeffreys"),
  seed_condition = "rds_from == -999",
  prefix = "rds",
  n_boot = 100,
  parallel_boot = FALSE,
  label = "chords"
) {
  if (parallel_boot) {
    requireNamespace(c("doParallel", "parallel"))
    doParallel::registerDoParallel(cores = parallel::detectCores() - 1)
  }

  type <- match.arg(type)
  .K <- log2(length(grep(pattern = "^type_.*_out$", names(data))))
  .pattern <- paste0(
    "^type_",
    paste0(rep("[0-9]", .K - 1), collapse = ""),
    "1_visible_out$"
  )

  .data_mod <-
    data %>%
    dplyr::filter(dplyr::if_all(all_of(prefix), ~ . == 1)) %>%
    dplyr::mutate(
      NS1 = apply(
        .[, grep(pattern = .pattern, x = names(data)), with = FALSE],
        1,
        sum
      ),
      refCoupNum = get(paste0(prefix, "_own_coupon")),
      interviewDt = get(paste0(prefix, "_t"))
    ) %>%
    dplyr::rename_with(
      .cols = starts_with(paste0(prefix, "_coupon_")),
      ~ gsub(pattern = paste0(prefix, "\\_coupon\\_"), replacement = "coup", .)
    )

  # if (type == "leave-d-out") {
  #   .jack_control <- chords::makeJackControl(1, 1e2)
  #
  #   .fit_chords <-
  #     chords::Estimate.b.k(rds.object = chords::initializeRdsObject(.data_mod),
  #                          type = type,
  #                          jack.control = .jack_control) %>%
  #     {c(est = sum(.$estimates$Nk.estimates[.$estimates$Nk.estimates < Inf]),
  #        degree_hidden =
  #          stats::weighted.mean(
  #            x = as.numeric(names(.$estimates$Nk.estimates))[.$estimates$Nk.estimates < Inf],
  #            w = .$estimates$Nk.estimates[.$estimates$Nk.estimates < Inf]))}
  #
  # } else {

  .fit_chords <-
    suppressWarnings(
      suppressMessages(
        chords::Estimate.b.k(
          rds.object = chords::initializeRdsObject(.data_mod),
          type = type
        ) %>%
          {
            c(
              est = sum(.$estimates$Nk.estimates[
                .$estimates$Nk.estimates < Inf
              ]),
              degree_hidden = stats::weighted.mean(
                x = as.numeric(names(.$estimates$Nk.estimates))[
                  .$estimates$Nk.estimates < Inf
                ],
                w = .$estimates$Nk.estimates[.$estimates$Nk.estimates < Inf]
              )
            )
          }
      )
    )

  # }

  .fit_chords_boot <-
    get_rds_boot(
      data = .data_mod,
      seed_condition = seed_condition,
      in_coupon = "refCoupNum",
      out_coupon = "coup",
      trait_var = "NS1",
      other_vars = c("NS1", "interviewDt", "hidden_visible_out", "name"),
      n_boot = n_boot
    ) %>%
    plyr::laply(
      .,
      function(samp) {
        suppressWarnings(
          suppressMessages(
            chords::Estimate.b.k(
              rds.object = chords::initializeRdsObject(samp),
              type = type
            ) %>%
              {
                c(
                  est = sum(.$estimates$Nk.estimates[
                    .$estimates$Nk.estimates < Inf
                  ]),
                  degree_hidden = stats::weighted.mean(
                    x = as.numeric(names(.$estimates$Nk.estimates))[
                      .$estimates$Nk.estimates < Inf
                    ],
                    w = .$estimates$Nk.estimates[.$estimates$Nk.estimates < Inf]
                  )
                )
              }
          )
        )
      },
      .parallel = parallel_boot
    )

  return(
    # data.frame(estimator = paste0(c("hidden_size_", "degree_hidden_"), label),
    #            estimate = c(.fit_chords["est"], .fit_chords["degree_hidden"]),
    #            se =  apply(.fit_chords_boot, 2, sd, na.rm = TRUE),
    #            inquiry = c("hidden_size", "degree_hidden"))
    data.frame(
      estimator = paste0(c("hidden_size_"), label),
      estimate = c(.fit_chords["est"]),
      se = sd(.fit_chords_boot[, 1], na.rm = TRUE),
      inquiry = c("hidden_size")
    )
  )
}

#' Estimator of population size based on Link-Tracing Sample by Vincent and Thompson
#'
#' @param data pass through sample
#' @param total integer giving the total size of the population
#' @param strata string specifying column name of strata vector
#' @param gibbs_params named list of parameters passed to Gibbs sampler
#' @param priors named list of prior specification for population size, stratum membership and links. p_n is an integer specifying the power law prior for population size (0 = flat). p_l is a positive rational numeric vector of length n_strata specifying the dirichlet prior for stratum membership (0.1 = non-informative). p_b is an integer specifying the beta distribution prior for links (1 = non-informative).
#' @param progress logical indicating whether to display progress bar. Defaults to \code{FALSE}
#' @param prefix character string giving name of the column with RDS+ sampling indicator
#' @param label character string giving label for the estimator
#' @return Data frame of link tracing estimates for single study
#'
#' @export
#'
#' @references Vincent, Kyle, and Steve Thompson. "Estimating the size and distribution of networked populations with snowball sampling." Journal of Survey Statistics and Methodology 10.2 (2022): 397-418.
#'
#' @import data.table
#' @importFrom magrittr `%>%` `%$%`
get_study_est_linktrace <- function(
  data,
  total = 2000,
  strata,
  gibbs_params = list(
    n_samples = 50L,
    chain_samples = 250L,
    chain_burnin = 50L
  ),
  priors = list(p_n = 0L, p_l = 0.1, p_b = 1L),
  progress = FALSE,
  prefix = "lts",
  label = "link"
) {
  varname_map <-
    names(data)[which(grepl(paste0("^", prefix, "\\_"), names(data)))]

  names(varname_map) <- gsub(
    pattern = paste0("^", prefix),
    replacement = "",
    x = varname_map
  )

  data.table::setnames(
    data,
    old = varname_map,
    new = names(varname_map)
  )

  # switch node names to sample only
  data <-
    data[
      get(prefix) == 1,
    ][,
      c("_from", "name", "strata_id") := list(
        plyr::mapvalues(
          `_from`,
          from = c(-999, as.numeric(name)),
          to = c(NA, 1:.N),
          warn_missing = F
        ),
        plyr::mapvalues(
          as.numeric(name),
          from = as.numeric(name),
          to = 1:.N,
          warn_missing = F
        ),
        as.numeric(factor(get(strata)))
      )
    ]

  .samp_graph <-
    igraph::graph_from_data_frame(
      data[!is.na(name) & !is.na(`_from`), .(`_from`, name)],
      vertices = data$name,
      directed = FALSE
    )

  # add network data for RDS sample
  data <-
    igraph::as_adj_list(.samp_graph) %>%
    lapply(function(i) as.numeric(i$name)) %>%
    data.table::data.table(name = data$name, links_list = .) %>%
    data[., on = "name"]

  y_samp <- as.matrix(igraph::as_adjacency_matrix(.samp_graph))

  n_strata <- length(unique(data$strata_id))

  if (length(priors$p_l) == 1) {
    priors$p_l <- rep(priors$p_l, n_strata)
  }

  if (n_strata != length(priors$p_l)) {
    stop(
      "mismatch between number of strata and number of priors specified for strata in p_l"
    )
  }

  res <- hiddenmeta:::lt_gibbs_cpp(
    links_list = data$links_list,
    wave = data$`_wave`,
    name = data$name,
    y_samp = y_samp,
    strata = data$strata_id,
    n_strata = n_strata,
    n_waves = max(data$`_wave`),
    total = total,
    chain_samples = gibbs_params$chain_samples,
    chain_burnin = gibbs_params$chain_burnin,
    prior_n = priors$p_n,
    prior_l = priors$p_l,
    prior_b = priors$p_b,
    n_0 = total,
    l_0 = rep(1 / n_strata, n_strata),
    b_0 = matrix(rep(0.1, n_strata * n_strata), n_strata, n_strata),
    n_samples = gibbs_params$n_samples,
    progress = progress
  )

  colnames(res$L) <- unique(data[[strata]])

  return(
    data.frame(
      estimator = paste0("hidden_size_", label),
      estimate = mean(res$N),
      se = sd(res$N),
      inquiry = "hidden_size"
    )
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
  parallel_boot = FALSE
) {
  if (parallel_boot) {
    requireNamespace(c("doParallel", "parallel"))
    doParallel::registerDoParallel(cores = parallel::detectCores() - 1)
  }
  .data_mod <- data[get(prefix) == 1, ]

  if (!is.null(weight_var)) {
    .weights <- .data_mod[[weight_var]]
  } else {
    .weights <- 1
  }

  .known_pops <-
    unlist(.data_mod[,
      lapply(.SD, unique),
      .SDcols = paste0(
        "total_",
        gsub("_visible_out", replacement = "", x = known_vars)
      )
    ])

  .data_est <-
    .data_mod[,
      c(known_vars, hidden_var),
      with = FALSE
    ][,
      lapply(.SD, function(x) x * .weights)
    ]

  .fit_nsum_est <-
    .data_est %>%
    networkscaleup::killworth(
      .,
      known_sizes = .known_pops,
      known_ind = 1:(ncol(.) - length(hidden_var)),
      N = total,
      model = "MLE"
    ) %>%
    .$sizes

  .fit_nsum_boot_se <-
    get_rescaled_boot(
      data = .data_mod,
      survey_design = survey_design,
      n_boot = n_boot
    ) %>%
    plyr::llply(
      .data = .,
      .fun = function(wgt) {
        .data_est[, lapply(.SD, function(x) x * wgt[, 2])] %>%
          networkscaleup::killworth(
            .,
            known_sizes = .known_pops,
            known_ind = 1:(ncol(.) - length(hidden_var)),
            N = total,
            model = "MLE"
          ) %>%
          .$sizes
      },
      .parallel = parallel_boot
    ) %>%
    do.call("c", .)

  return(
    data.frame(
      estimator = paste0(c("hidden_size_"), label),
      estimate = .fit_nsum_est,
      se = sd(.fit_nsum_boot_se),
      inquiry = c("hidden_size")
    )
  )
}

#' Generalized NSUM estimator (Feehan & Salganik 2016)
#'
#' @param data pass-through population data frame
#' @param hidden_var character vector containing names of hidden groups
#' @param known_vars character vector containing names of known groups (for reference)
#' @param total integer giving total size of population
#' @param prefix character prefix used for sample variables
#' @param label character string describing the estimator
#' @param weight_var character string giving name of the sampling weights variable
#' @param survey_design a formula describing the design of the survey
#' @param n_boot integer giving number of bootstrap re-samples
#' @param parallel_boot logical, whether to compute bootstrap samples in parallel
#'
#' @return Data frame of generalized NSUM estimates
#' @export
#'
#' @references Feehan, Dennis M., and Matthew J. Salganik. "Generalizing the network scale-up method: a new estimator for the size of hidden populations." Sociological Methodology 46.1 (2016): 153-186.
#'
#' @import surveybootstrap
#' @importFrom magrittr `%>%` `%$%`
#' @importFrom plyr llply
get_study_est_gnsum <- function(
  data,
  hidden_var = "hidden_visible_out",
  known_vars = paste0(c("known", paste0("known_", 2:10)), "_visible_out"),
  total = 1000,
  prefix = "pps",
  label = "gnsum",
  weight_var = "pps_weight",
  survey_design = ~ pps_cluster + strata(pps_strata),
  n_boot = 1000,
  parallel_boot = FALSE
) {
  if (parallel_boot) {
    requireNamespace(c("doParallel", "parallel"))
    doParallel::registerDoParallel(cores = parallel::detectCores() - 1)
  }

  # ---- Frame sample ----
  .data_F <- data[get(prefix) == 1, ]

  # Weights
  if (!is.null(weight_var)) {
    .w <- .data_F[[weight_var]]
  } else {
    .w <- rep(1, nrow(.data_F))
  }

  # ---- STEP 1: y_{F,H} = weighted total out-reports from frame ----
  y_F_H <- sum(.data_F[[hidden_var]] * .w, na.rm = TRUE)

  # ---- STEP 2: v̄_{H,F} = mean in-degree (visibility) among hidden members ----
  visibility_var <- gsub("_out$", "_in", hidden_var)

  .data_H <- data[hidden == 1, ]

  if (!visibility_var %in% names(.data_H)) {
    stop(
      "Visibility variable '",
      visibility_var,
      "' not found. ",
      "Simplified gNSUM requires a population-level hidden visibility variable."
    )
  }

  if (nrow(.data_H) == 0) {
    warning("No hidden population rows found; returning NA.")
    v_bar_H_F <- NA
  } else {
    v_bar_H_F <- mean(.data_H[[visibility_var]], na.rm = TRUE)
  }

  # ---- STEP 3: simplified generalized NSUM estimate ----
  if (is.na(v_bar_H_F) || v_bar_H_F == 0) {
    N_H_est <- NA
    prev_est <- NA
    warning("Cannot compute simplified gNSUM: visibility is zero or NA.")
  } else {
    N_H_est <- y_F_H / v_bar_H_F
    prev_est <- N_H_est / total
  }

  # ---- STEP 4: Bootstrap SE (resampled frame weights only) ----
  if (is.na(N_H_est)) {
    size_se <- NA
    prev_se <- NA
  } else {
    boot_vals <-
      get_rescaled_boot(
        data = .data_F,
        survey_design = survey_design,
        n_boot = n_boot
      ) %>%
      plyr::llply(
        .fun = function(wgt) {
          y_F_H_boot <- sum(.data_F[[hidden_var]] * wgt[, 2], na.rm = TRUE)
          N_H_boot <- y_F_H_boot / v_bar_H_F # visibility held fixed (simulation truth)
          N_H_boot
        },
        .parallel = parallel_boot
      ) %>%
      unlist()

    size_se <- sd(boot_vals, na.rm = TRUE)
    prev_se <- size_se / total
  }

  # ---- Return ----
  data.frame(
    estimator = paste0(c("hidden_size_", "hidden_prev_"), label),
    estimate = c(N_H_est, prev_est),
    se = c(size_se, prev_se),
    inquiry = c("hidden_size", "hidden_prev")
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
get_study_est_multiplier <- function(
  data,
  service_var = "service_use",
  total_service = sum(data$service_use[data$hidden == 1]),
  seed_condition = "rds_from == -999",
  n_boot = 100,
  parallel_boot = FALSE,
  prefix = "rds",
  label = "multiplier"
) {
  if (parallel_boot) {
    requireNamespace(c("doParallel", "parallel"))
    doParallel::registerDoParallel(cores = parallel::detectCores() - 1)
  }

  if (!data.table::is.data.table(data)) {
    data <- data.table::as.data.table(data)
  }

  .data_mod <- data[get(prefix) == 1, ]

  .est_out <-
    total_service / mean(.data_mod[[service_var]])

  .est_boot <-
    get_rds_boot(
      data = .data_mod,
      seed_condition = seed_condition,
      in_coupon = paste0(prefix, "_own_coupon"),
      out_coupon = paste0(prefix, "_coupon_"),
      trait_var = "hidden_visible_out",
      other_vars = c("hidden_visible_out", service_var, "name"),
      n_boot = n_boot
    ) %>%
    plyr::laply(., function(samp) {
      total_service / mean(samp[[service_var]])
    })

  return(
    data.frame(
      estimator = paste0("hidden_size_", label),
      estimate = .est_out,
      se = sd(.est_boot),
      inquiry = "hidden_size"
    )
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
  if (!data.table::is.data.table(data)) {
    data <- data.table::as.data.table(data)
  }

  if (!is.null(sample_condition)) {
    data <- data[eval(parse(text = sample_condition)), ]
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

      .est_out <-
        .est_out |>
        dplyr::select(all_of(capture_vars)) |>
        Rcapture::periodhist(vt = .pool) |>
        Rcapture::closedp.bc(
          dfreq = TRUE,
          dtype = "hist",
          # t = length(capture_vars),
          m = model
        )
    } else {
      .est_out <-
        .est_out |>
        dplyr::select(all_of(capture_vars)) |>
        Rcapture::closedp.bc(
          dfreq = FALSE,
          dtype = "hist",
          # t = length(capture_vars),
          m = model
        )
    }

    return(
      data.frame(
        estimator = paste0("hidden_size_", label),
        estimate = .est_out$results[1],
        se = .est_out$results[2],
        inquiry = "hidden_size"
      )
    )
  }
}
