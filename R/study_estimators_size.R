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
get_study_est_sspse <- function(data,
                                prior_mean = .1 * nrow(data),
                                n_coupons = 3,
                                total = 2000,
                                mcmc_params = list(interval = 5,
                                                   burnin = 2000,
                                                   samplesize = 500),
                                additional_params = list(),
                                prefix = "rds",
                                label = "sspse") {

  .quiet_sspse <- purrr::quietly(sspse::posteriorsize)

  network_sizes <- data %>%
    dplyr::filter(dplyr::if_all(all_of(prefix), ~ . == 1)) %>%
    dplyr::pull("hidden_visible_out")

  # Check 1: Empty sample
  if (length(network_sizes) == 0) {
    warning("No RDS sample found. Returning NA.")
    return(
      data.frame(estimator = paste0("hidden_size_", label),
                 estimate = NA_real_,
                 se = NA_real_,
                 inquiry = c("hidden_size"))
    )
  }

  # Check 2: Too sparse (>75% zeros or max <= 1)
  prop_zero <- mean(network_sizes == 0)
  if (prop_zero > 0.90 || max(network_sizes) <= 1) {
    warning(paste0("Network sizes too sparse (", round(prop_zero*100, 1),
                   "% zeros, max=", max(network_sizes),
                   "). SS-PSE estimate unreliable. Returning NA."))
    return(
      data.frame(estimator = paste0("hidden_size_", label),
                 estimate = NA_real_,
                 se = NA_real_,
                 inquiry = c("hidden_size"))
    )
  }

  # Calculate K with minimum of 2 (need at least 2 for cutabove to work)
  K_value <- max(2, round(stats::quantile(network_sizes, 0.80)))

  # Check 3: Verify we have variation
  if (length(unique(network_sizes)) < 3) {
    warning("Insufficient variation in network sizes. Returning NA.")
    return(
      data.frame(estimator = paste0("hidden_size_", label),
                 estimate = NA_real_,
                 se = NA_real_,
                 inquiry = c("hidden_size"))
    )
  }

  # Try estimation with error handling
  .fit_sspse <- tryCatch({
    do.call(
      .quiet_sspse,
      c(list(s = network_sizes,
             interval = mcmc_params$interval,
             samplesize = mcmc_params$samplesize,
             burnin = mcmc_params$burnin,
             mean.prior.size = prior_mean,
             verbose = FALSE,
             K = K_value,
             max.coupons = n_coupons),
        additional_params)
    )
  }, error = function(e) {
    warning("SS-PSE estimation failed: ", e$message)
    return(NULL)
  })

  # Check if estimation failed
  if (is.null(.fit_sspse) || is.null(.fit_sspse$result)) {
    return(
      data.frame(estimator = paste0("hidden_size_", label),
                 estimate = NA_real_,
                 se = NA_real_,
                 inquiry = c("hidden_size"))
    )
  }

  return(
    data.frame(estimator = paste0("hidden_size_", label),
               estimate = c(unname(.fit_sspse$result$N["Median AP"])),
               se = c(sd(.fit_sspse$result$sample[,"N"])),
               inquiry = c("hidden_size"))
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
get_study_est_chords <- function(data,
                                 type = c("mle", "integrated", "jeffreys"),
                                 seed_condition = "rds_from == -999",
                                 prefix = "rds",
                                 n_boot = 100,
                                 parallel_boot = FALSE,
                                 label = "chords") {

  if (parallel_boot) {
    requireNamespace(c("doParallel", "parallel"))
    doParallel::registerDoParallel(cores = parallel::detectCores() - 1)
  }

  type <- match.arg(type)
  .K <- log2(length(grep(pattern = "^type_.*_out$", names(data))))
  .pattern <- paste0("^type_", paste0(rep("[0-9]", .K - 1), collapse = ""), "1_visible_out$")

  .data_mod <-
    data %>%
    dplyr::filter(dplyr::if_all(all_of(prefix), ~ . == 1)) %>%
    dplyr::mutate(
      NS1 = apply(.[, grep(pattern = .pattern, x = names(data)), with = FALSE], 1, sum),
      refCoupNum = get(paste0(prefix, "_own_coupon")),
      interviewDt = get(paste0(prefix, "_t"))) %>%
    dplyr::rename_with(
        .cols = starts_with(paste0(prefix, "_coupon_")),
        ~ gsub(pattern = paste0(prefix, "\\_coupon\\_"), replacement = "coup", .))

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
        chords::Estimate.b.k(rds.object = chords::initializeRdsObject(.data_mod),
                             type = type) %>%
          {c(est = sum(.$estimates$Nk.estimates[.$estimates$Nk.estimates < Inf]),
             degree_hidden =
               stats::weighted.mean(
                 x = as.numeric(names(.$estimates$Nk.estimates))[.$estimates$Nk.estimates < Inf],
                 w = .$estimates$Nk.estimates[.$estimates$Nk.estimates < Inf]))}
      )
    )

  # }

  .fit_chords_boot <-
    get_rds_boot(data = .data_mod,
                 seed_condition = seed_condition,
                 in_coupon = "refCoupNum",
                 out_coupon = "coup",
                 trait_var = "NS1",
                 other_vars = c("NS1", "interviewDt", "hidden_visible_out", "name"),
                 n_boot = n_boot) %>%
    plyr::laply(., function(samp) {
      suppressWarnings(
        suppressMessages(
          chords::Estimate.b.k(rds.object = chords::initializeRdsObject(samp),
                               type = type) %>%
            {c(est = sum(.$estimates$Nk.estimates[.$estimates$Nk.estimates < Inf]),
               degree_hidden =
                 stats::weighted.mean(
                   x = as.numeric(names(.$estimates$Nk.estimates))[.$estimates$Nk.estimates < Inf],
                   w = .$estimates$Nk.estimates[.$estimates$Nk.estimates < Inf]))}
        )
      )
    },
    .parallel = parallel_boot)

  return(
    # data.frame(estimator = paste0(c("hidden_size_", "degree_hidden_"), label),
    #            estimate = c(.fit_chords["est"], .fit_chords["degree_hidden"]),
    #            se =  apply(.fit_chords_boot, 2, sd, na.rm = TRUE),
    #            inquiry = c("hidden_size", "degree_hidden"))
    data.frame(estimator = paste0(c("hidden_size_"), label),
               estimate = c(.fit_chords["est"]),
               se =  sd(.fit_chords_boot[,1], na.rm = TRUE),
               inquiry = c("hidden_size"))
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
    gibbs_params = list(n_samples = 50L, chain_samples = 250L, chain_burnin = 50L),
    priors = list(p_n = 0L, p_l = 0.1, p_b = 1L),
    progress = FALSE,
    prefix = "lts",
    label = "link"
){

  varname_map <-
    names(data)[which(grepl(paste0("^", prefix, "\\_"), names(data)))]

  names(varname_map) <- gsub(pattern = paste0("^", prefix),
                             replacement = "", x = varname_map)

  data.table::setnames(
    data,
    old = varname_map,
    new = names(varname_map)
  )

  # switch node names to sample only
  data <-
    data[
      get(prefix) == 1,
    ][
      , c("_from", "name", "strata_id") :=
        list(plyr::mapvalues(`_from`,
                             from = c(-999, as.numeric(name)), to = c(NA, 1:.N),
                             warn_missing = F),
             plyr::mapvalues(as.numeric(name),
                             from = as.numeric(name), to = 1:.N,
                             warn_missing = F),
             as.numeric(factor(get(strata))))
    ]

  .samp_graph <-
    igraph::graph_from_data_frame(
      data[!is.na(name) & !is.na(`_from`), .(`_from`, name)],
      vertices = data$name,
      directed = FALSE)

  # add network data for RDS sample
  data <-
    igraph::as_adj_list(.samp_graph) %>%
    lapply(function(i) as.numeric(i$name)) %>%
    data.table::data.table(name = data$name, links_list = .) %>%
    data[., on = "name"]

  y_samp <- as.matrix(igraph::as_adjacency_matrix(.samp_graph))

  n_strata <- length(unique(data$strata_id))

  if(length(priors$p_l) == 1){
    priors$p_l <- rep(priors$p_l, n_strata)
  }

  if(n_strata != length(priors$p_l)){
    stop("mismatch between number of strata and number of priors specified for strata in p_l")
  }

  res <- lt_gibbs_cpp(links_list = data$links_list,
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
                      l_0 = rep(1/n_strata, n_strata),
                      b_0 = matrix(rep(0.1, n_strata * n_strata), n_strata, n_strata),
                      n_samples = gibbs_params$n_samples,
                      progress = progress)

  colnames(res$L) <- unique(data[[strata]])

  return(
    data.frame(estimator = paste0("hidden_size_", label),
               estimate = mean(res$N),
               se = sd(res$N),
               inquiry = "hidden_size")
  )

}
