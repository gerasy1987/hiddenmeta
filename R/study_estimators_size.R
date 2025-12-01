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
  mcmc_params = list(interval = 5, burnin = 2000, samplesize = 500),
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

  # Check 2: Too sparse (>75% zeros or max <= 1)
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

  # Calculate K with minimum of 2 (need at least 2 for cutabove to work)
  K_value <- max(2, round(stats::quantile(network_sizes, 0.80)))

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
            burnin = mcmc_params$burnin,
            mean.prior.size = prior_mean,
            verbose = FALSE,
            K = K_value,
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
  # Create a copy to avoid modifying original data
  data <- data.table::copy(data)

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

  # Filter to sample only first, then process
  data <- data[get(prefix) == 1, ]

  # Check if strata column exists and has no NAs
  if (!strata %in% names(data)) {
    stop(paste0("Strata column '", strata, "' not found in data"))
  }

  if (any(is.na(data[[strata]]))) {
    warning(paste0(
      "NAs found in strata column '",
      strata,
      "'. Removing affected rows."
    ))
    data <- data[!is.na(get(strata)), ]
  }

  # Create strata_id and ensure it's valid
  data[, strata_id := as.numeric(factor(get(strata)))]

  # Ensure name column exists and is valid
  if (!"name" %in% names(data)) {
    stop("'name' column not found in data")
  }

  # Create sequential IDs for the sample
  data[, sample_id := 1:.N]

  # Map the _from variable
  if ("_from" %in% names(data)) {
    # Create mapping from original names to sample IDs
    name_to_id <- setNames(data$sample_id, data$name)

    data[,
      `_from` := ifelse(
        `_from` == -999,
        NA,
        name_to_id[as.character(`_from`)]
      )
    ]
  } else {
    stop("'_from' column not found after renaming")
  }

  # Build graph only from valid edges
  valid_edges <- data[
    !is.na(`_from`) & !is.na(sample_id),
    .(`_from`, sample_id)
  ]

  if (nrow(valid_edges) == 0) {
    warning("No valid edges found in link-tracing sample")
    return(data.frame(
      estimator = paste0("hidden_size_", label),
      estimate = NA,
      se = NA,
      inquiry = "hidden_size"
    ))
  }

  .samp_graph <-
    igraph::graph_from_data_frame(
      valid_edges,
      vertices = data$sample_id,
      directed = FALSE
    )

  # Get adjacency list
  adj_list <- igraph::as_adj_list(.samp_graph)

  # Convert to numeric list of links
  data[,
    links_list := lapply(sample_id, function(id) {
      if (id %in% names(adj_list)) {
        as.numeric(names(adj_list[[as.character(id)]]))
      } else {
        numeric(0)
      }
    })
  ]

  # Get adjacency matrix
  y_samp <- as.matrix(igraph::as_adjacency_matrix(.samp_graph))

  # Set up strata information
  n_strata <- length(unique(data$strata_id))

  if (length(priors$p_l) == 1) {
    priors$p_l <- rep(priors$p_l, n_strata)
  }

  if (n_strata != length(priors$p_l)) {
    stop(
      "mismatch between number of strata and number of priors specified for strata in p_l"
    )
  }

  # Check for wave variable
  if (!"_wave" %in% names(data)) {
    stop("'_wave' column not found after renaming")
  }

  if (any(is.na(data$`_wave`))) {
    stop("NAs found in wave variable")
  }

  # DIAGNOSTIC OUTPUT - Print everything being passed to C++
  cat("\n=== DIAGNOSTIC INFO ===\n")
  cat("Sample size (n):", nrow(data), "\n")
  cat("Number of strata:", n_strata, "\n")
  cat("Strata IDs:", paste(unique(data$strata_id), collapse = ", "), "\n")
  cat(
    "Strata distribution:",
    paste(table(data$strata_id), collapse = ", "),
    "\n"
  )
  cat("Wave range:", min(data$`_wave`), "to", max(data$`_wave`), "\n")
  cat("Number of waves:", max(data$`_wave`), "\n")
  cat("Total population:", total, "\n")
  cat("Chain samples:", gibbs_params$chain_samples, "\n")
  cat("Chain burnin:", gibbs_params$chain_burnin, "\n")
  cat("Number of samples:", gibbs_params$n_samples, "\n")
  cat("\nPriors:\n")
  cat("  p_n:", priors$p_n, "\n")
  cat("  p_l:", paste(priors$p_l, collapse = ", "), "\n")
  cat("  p_b:", priors$p_b, "\n")
  cat("\nInitial values:\n")
  cat("  n_0:", total, "\n")
  cat("  l_0:", paste(rep(1 / n_strata, n_strata), collapse = ", "), "\n")
  cat("  b_0 dimensions:", n_strata, "x", n_strata, "\n")
  cat("\nData structure checks:\n")
  cat(
    "  Any NA in links_list:",
    any(sapply(data$links_list, function(x) any(is.na(x)))),
    "\n"
  )
  cat("  Any NA in wave:", any(is.na(data$`_wave`)), "\n")
  cat("  Any NA in sample_id:", any(is.na(data$sample_id)), "\n")
  cat("  Any NA in strata_id:", any(is.na(data$strata_id)), "\n")
  cat("  Dimensions of y_samp:", paste(dim(y_samp), collapse = " x "), "\n")
  cat("  Any NA in y_samp:", any(is.na(y_samp)), "\n")

  # Check if any links point to non-existent nodes
  all_links <- unlist(data$links_list)
  if (length(all_links) > 0) {
    cat("  Link range:", min(all_links), "to", max(all_links), "\n")
    cat(
      "  Links outside sample range:",
      any(all_links < 1 | all_links > nrow(data)),
      "\n"
    )
  }

  # Check the ratio of sample to total
  cat("\nSample to total ratio:", round(nrow(data) / total, 3), "\n")

  cat("======================\n\n")

  # Save inputs for inspection
  debug_data <- list(
    links_list = data$links_list,
    wave = data$`_wave`,
    name = data$sample_id,
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

  # Save to temp file for inspection
  temp_file <- tempfile(fileext = ".rds")
  saveRDS(debug_data, temp_file)
  cat("Debug data saved to:", temp_file, "\n\n")

  # Call C++ function
  res <- lt_gibbs_cpp(
    links_list = data$links_list,
    wave = data$`_wave`,
    name = data$sample_id,
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
    data.frame(
      estimator = paste0("hidden_prev_", label),
      estimate = NA,
      se = NA,
      inquiry = "hidden_prev"
    )
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
