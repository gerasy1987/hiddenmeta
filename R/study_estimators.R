#' SS-PSE population size estimator by Handcock, Gile and Mar
#'
#' @param data pass-through population data frame
#' @param prior_mean the mean of the prior distribution on the population size for SS-PSE estimation
#' @param n_coupons The number of recruitment coupons distributed to each enrolled subject (i.e. the maximum number of recruitees for any subject). By default it is taken by the attribute or data, else the maximum recorded number of coupons.
#' @param mcmc_params named list of parameters passed to \code{sspse::posteriorsize} for MCMC sampling
#' @param rds_prefix character prefix used for RDS sample variables
#' @param label character string describing the estimator
#'
#' @return Data frame of SS-PSE estimates for a single study
#'
#' @references
#' Handcock, Mark S., Krista J. Gile, and Corinne M. Mar. “Estimating Hidden Population Size Using Respondent-Driven Sampling Data.” Electronic Journal of Statistics 8, no. 1 (2014): 1491–1521. \url{https://doi.org/10.1214/14-EJS923}.
#' @export
#'
#' @import dplyr
#' @importFrom sspse posteriorsize
#' @importFrom RDS as.rds.data.frame
#' @importFrom purrr quietly
get_study_est_sspse <- function(data,
                                prior_mean = .1 * nrow(data),
                                n_coupons = 3,
                                mcmc_params = list(interval = 5, burnin = 2000, samplesize = 500),
                                rds_prefix = "rds",
                                label = "sspse") {

  .quiet_sspse <- purrr::quietly(sspse::posteriorsize)

  .fit_sspse <-
    data %>%
    dplyr::filter(dplyr::if_all(dplyr::all_of(rds_prefix), ~ . == 1)) %>%
    dplyr::select(name, hidden_visible_out, starts_with(rds_prefix), total) %>%
    RDS::as.rds.data.frame(df = .,
                           id = "name",
                           recruiter.id = paste0(rds_prefix, "_from"),
                           network.size = "hidden_visible_out",
                           time = paste0(rds_prefix, "_t"),
                           population.size = unique(.$total),
                           max.coupons = n_coupons) %>%
    {
      .quiet_sspse(s = .,
                   interval = mcmc_params$interval,
                   samplesize = mcmc_params$samplesize,
                   burnin = mcmc_params$burnin,
                   mean.prior.size = prior_mean,
                   verbose = FALSE,
                   # visibility = TRUE,
                   max.coupons = n_coupons
      )
    }

  return(
    data.frame(estimator_label = paste0("hidden_size_", label),
               estimate = c(unname(.fit_sspse$result$N["Median AP"])),
               se =   c(sd(.fit_sspse$result$sample[,"N"])),
               inquiry_label = c("hidden_size"))
  )
}

#' Horvitz-Thompson prevalence estimator using
#'
#' @param data pass-through population data frame
#' @param hidden_var variable containing hidden group membership indicator
#' @param weight_var variable containing sampling weights
#' @param total_var variable containing size of population for prevalence estimation
#' @param survey_design a formula describing the design of the survey (for bootstrap)
#' @param n_boot number of bootstrap resamples
#' @param parallel_boot logical, whether to compute bootstrap samples in parallel using \code{foreach} package
#' @param prefix character prefix used for sampling variables (has to include \code{[prefix]_weights})
#' @param label character string describing the estimator
#'
#' @return Data frame of HT estimates for a single study
#' @export
#'
#' @references
#' Salganik, Matthew J. "Variance estimation, design effects, and sample size calculations for respondent-driven sampling." Journal of Urban Health 83, no. 1 (2006): 98. \url{https://doi.org/10.1007/s11524-006-9106-x}
#'
#' @import dplyr
get_study_est_ht <- function(data,
                             hidden_var = "hidden",
                             weight_var = "pps_weight",
                             total_var = "total",
                             survey_design = ~ pps_cluster + strata(pps_strata),
                             n_boot = 100,
                             parallel_boot = TRUE,
                             prefix = "pps",
                             label = "ht") {

  if (parallel_boot) {
    requireNamespace(c("doParallel", "parallel"))
    doParallel::registerDoParallel(cores = parallel::detectCores() - 1)
  }

  .data_mod <-
    data %>%
    dplyr::filter(dplyr::if_all(dplyr::all_of(prefix), ~ . == 1))

  if (any(is.na(.data_mod[[weight_var]])))
    stop("there are missing values in sampling weights provided")
  if (any(is.na(.data_mod[[hidden_var]])))
    stop("there are missing values in hidden population indicator provided")

  .est_ht <-
    crossprod(.data_mod[[hidden_var]],
              .data_mod[[weight_var]])

  .est_ht_boot <-
    get_rescaled_boot(data = .data_mod,
                      survey_design = survey_design,
                      n_boot = n_boot) %>%
    plyr::laply(., function(samp) {
      .data_mod %>%
        dplyr::mutate(rescale_wgt = samp$weight.scale) %>%
        dplyr::filter(rescale_wgt != 0) %>%
        {
          crossprod(.[[hidden_var]],
                    .$rescale_wgt * .[[weight_var]])
        }
    }, .parallel = parallel_boot)

  # .fit_ht <-
  #   data %>%
  #   dplyr::filter(dplyr::if_all(dplyr::all_of(prefix), ~ . == 1)) %>%
  #   estimatr::horvitz_thompson(hidden ~ 1, weights = get(paste0(prefix, "_weight")), data = .) %>%
  #   {c(est = unname(.$coefficients), se = unname(.$std.error))}

 return(
   data.frame(estimator_label = paste0(c("hidden_size", "hidden_prev"), "_", label),
              estimate = c(.est_ht, .est_ht/unique(.data_mod[[total_var]])),
              se =  c(sd(.est_ht_boot), sd(.est_ht_boot/unique(.data_mod[[total_var]]))),
              inquiry_label = c("hidden_size", "hidden_prev"))
 )
}


#' Chords population size estimatior by Berchenko, Rosenblatt and Frost
#'
#' @param data pass-through population data frame
#' @param type a character vector with the type of estimation. Can be one of \code{mle}, \code{integrated}, or \code{jeffreys}. See \code{?chords::Estimate.b.k} and the original paper from the references for details
#' @param seed_condition character string containing condition to define seeds. Defaults to "rds_from == -999" that applies to simulated RDS samples
#' @param n_boot number of bootstrap resamples
#' @param parallel_boot logical, whether to compute bootstrap samples in parallel using \code{foreach} package
#' @param rds_prefix character prefix used for RDS sample variables
#' @param label character string describing the estimator
#'
#' @return Data frame of Chords estimates for a single study with RDS sample
#'
#' @references
#' Berchenko, Yakir, Jonathan D. Rosenblatt, and Simon D. W. Frost. “Modeling and Analyzing Respondent-Driven Sampling as a Counting Process.” Biometrics 73, no. 4 (2017): 1189–98. \url{https://doi.org/10.1111/biom.12678}.
#' Salganik, Matthew J. "Variance estimation, design effects, and sample size calculations for respondent-driven sampling." Journal of Urban Health 83, no. 1 (2006): 98. \url{https://doi.org/10.1007/s11524-006-9106-x}
#'
#' @export
#'
#' @import dplyr
#' @importFrom chords initializeRdsObject Estimate.b.k makeJackControl
#' @importFrom purrr quietly
get_study_est_chords <- function(data,
                                 type = c("mle", "integrated", "jeffreys"),
                                 seed_condition = "rds_from == -999",
                                 rds_prefix = "rds",
                                 n_boot = 100,
                                 parallel_boot = TRUE,
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
    dplyr::filter(dplyr::across(dplyr::all_of(rds_prefix), ~ . == 1)) %>%
    dplyr::mutate(
      NS1 = apply(.[,grep(pattern = .pattern, x = names(data))], 1, sum),
      refCoupNum = get(paste0(rds_prefix, "_own_coupon")),
      interviewDt = get(paste0(rds_prefix, "_t"))) %>%
    dplyr::rename_with(
        .cols = dplyr::starts_with(paste0(rds_prefix, "_coupon_")),
        ~ gsub(pattern = paste0(rds_prefix, "\\_coupon\\_"), replacement = "coup", .))

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
    chords::Estimate.b.k(rds.object = chords::initializeRdsObject(.data_mod),
                         type = type) %>%
    {c(est = sum(.$estimates$Nk.estimates[.$estimates$Nk.estimates < Inf]),
       degree_hidden =
         stats::weighted.mean(
           x = as.numeric(names(.$estimates$Nk.estimates))[.$estimates$Nk.estimates < Inf],
           w = .$estimates$Nk.estimates[.$estimates$Nk.estimates < Inf]))}

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
      suppressMessages(
        chords::Estimate.b.k(rds.object = chords::initializeRdsObject(samp),
                             type = type) %>%
          {c(est = sum(.$estimates$Nk.estimates[.$estimates$Nk.estimates < Inf]),
             degree_hidden =
               stats::weighted.mean(
                 x = as.numeric(names(.$estimates$Nk.estimates))[.$estimates$Nk.estimates < Inf],
                 w = .$estimates$Nk.estimates[.$estimates$Nk.estimates < Inf]))}
      )
    },
    .parallel = parallel_boot)

  return(
    data.frame(estimator_label = paste0(c("hidden_size_", "degree_hidden_"), label),
               estimate = c(.fit_chords["est"], .fit_chords["degree_hidden"]),
               se =  apply(.fit_chords_boot, 2, sd, na.rm = TRUE),
               inquiry_label = c("hidden_size", "degree_hidden_average"))
  )
}


#' NSUM estimatior
#'
#' @param data pass-through population data frame
#' @param known character vector containing names of known groups
#' @param hidden character vector containing names of hidden groups
#' @param survey_design a formula describing the design of the survey
#' @param degree_ratio numeric value between 0 and 1 representing degree ratio
#' @param transmission_rate numeric value between 0 and 1 representing information transmission rate
#' @param n_boot number of bootstrap resamples
#' @param parallel_boot logical, whether to compute bootstrap samples in parallel using \code{foreach} package
#' @param prefix character prefix used for PPS sample variables
#' @param label character string describing the estimator
#'
#' @return Data frame of NSUM estimates for a single study with PPS sample
#' @export
#'
#' @references
#' Dennis M. Feehan, Matthew J. Salganik. “The networkreporting package.” (2014). \url{https://cran.r-project.org/package=networkreporting}.
#' Dennis M. Feehan, Matthew J. Salganik. “The surveybootstrap package.” (2016). \url{https://cran.r-project.org/package=surveybootstrap}.
#' Salganik, Matthew J. "Variance estimation, design effects, and sample size calculations for respondent-driven sampling." Journal of Urban Health 83, no. 1 (2006): 98. \url{https://doi.org/10.1007/s11524-006-9106-x}
#'
#' @import dplyr
#' @importFrom networkreporting kp.degree.estimator nsum.estimator
#' @importFrom surveybootstrap bootstrap.estimates rescaled.bootstrap.sample
#' @importFrom purrr quietly
get_study_est_nsum <- function(data,
                               known,
                               hidden,
                               survey_design = ~ pps_cluster + strata(pps_strata),
                               degree_ratio = 1,
                               transmission_rate = 1,
                               n_boot = 1000,
                               parallel_boot = TRUE,
                               prefix = "pps",
                               label = "nsum") {

  if (parallel_boot) {
    requireNamespace(c("doParallel", "parallel"))
    doParallel::registerDoParallel(cores = parallel::detectCores() - 1)
  }

  .data_mod <-
    data %>%
    dplyr::filter(dplyr::if_all(dplyr::all_of(prefix), ~ . == 1))

  .known_pops <-
    .data_mod %>%
    dplyr::select(total, dplyr::all_of(paste0("total_", known))) %>%
    dplyr::rename_with(
      .cols = dplyr::starts_with("total_"), ~ paste0(gsub("^total\\_", "", .), "_visible_out")) %>%
    apply(., MARGIN = 2, unique)

  .data_mod$d_est <-
    purrr::quietly(networkreporting::kp.individual.estimator)(
      resp.data = .data_mod,
      alter.popn.size = .known_pops[1],
      known.populations = names(.known_pops[2:length(.known_pops)]),
      total.kp.size = sum(.known_pops[2:length(.known_pops)]))$result$dbar.Fcell.F

  .fit_nsum <-
    purrr::quietly(networkreporting::nsum.estimator)(
      survey.data = .data_mod,
      d.hat.vals = "d_est",
      total.popn.size = .known_pops[1],
      killworth.se = FALSE,
      y.vals = hidden,
      weights = paste0(prefix, "_weight"),
      missing = "ignore",
      deg.ratio = degree_ratio,
      tx.rate = transmission_rate)$result

  .fit_nsum_boot <-
    get_rescaled_boot(data = .data_mod,
                      survey_design = survey_design,
                      n_boot = n_boot) %>%
    plyr::llply(
      .data = .,
      .fun = function(wgt) {
        .data_mod %>%
          dplyr::mutate(index = 1:n()) %>%
          dplyr::left_join(., wgt, by = "index") %>%
          dplyr::mutate(rescaled_weight = weight.scale * get(paste0(prefix, "_weight"))) %>%
          networkreporting::nsum.estimator(
            survey.data = .,
            d.hat.vals = "d_est",
            total.popn.size = .known_pops[1],
            killworth.se = FALSE,
            y.vals = hidden,
            weights = "rescaled_weight",
            missing = "ignore",
            deg.ratio = degree_ratio,
            tx.rate = transmission_rate)
      },
      .parallel = parallel_boot) %>%
    bind_rows()

  return(
    data.frame(estimator_label = paste0(c("hidden_size_", "degree_"), label),
               estimate = c(unname(.fit_nsum$estimate), unname(.fit_nsum$sum.d.hat/.fit_nsum$estimate)),
               se = c(sd(.fit_nsum_boot$estimate), sd(.fit_nsum_boot$sum.d.hat/.fit_nsum_boot$estimate)),
               inquiry_label = c("hidden_size", "degree_average"))
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
#' @param rds_prefix character prefix used for RDS sample variables
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
#' @import dplyr
get_study_est_multiplier <- function(data,
                                     service_var = "service_use",
                                     total_service = sum(data$service_use[data$hidden == 1]),
                                     seed_condition = "rds_from == -999",
                                     n_boot = 100,
                                     parallel_boot = TRUE,
                                     rds_prefix = "rds",
                                     label = "multiplier") {

  if (parallel_boot) {
    requireNamespace(c("doParallel", "parallel"))
    doParallel::registerDoParallel(cores = parallel::detectCores() - 1)
  }

  .data_mod <-
    data %>%
    dplyr::filter(dplyr::across(dplyr::all_of(rds_prefix), ~ . == 1))

  .est_out <-
    total_service/mean(.data_mod[[service_var]])


  .est_boot <-
    get_rds_boot(data = .data_mod,
                 seed_condition = seed_condition,
                 in_coupon = paste0(rds_prefix, "_own_coupon"),
                 out_coupon = paste0(rds_prefix, "_coupon_"),
                 trait_var = "hidden_visible_out",
                 other_vars = c("hidden_visible_out", service_var, "name"),
                 n_boot = n_boot) %>%
    plyr::laply(., function(samp) {
      total_service/mean(samp[[service_var]])
    })

  return(
    data.frame(estimator_label = paste0("hidden_size_", label),
               estimate = .est_out,
               se =  sd(.est_boot),
               inquiry_label = "hidden_size")
  )
}

#' Generalized NSUM estimatior
#'
#' @param data pass-through population data frame
#' @param label character string describing the estimator
#'
#' @return Data frame of HT estimates for single study
#' @export
#'
#' @import dplyr
#' @importFrom estimatr lm_robust
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
    data.frame(estimator_label = paste0("hidden_prev_", label),
               estimate = fit_ht["est"],
               se =  fit_ht["se"],
               inquiry_label = "hidden_prev")
  )
}



#' Capture-recapture estimator for closed population
#'
#' @param data pass-through population data frame that contains capture indicators
#' @param capture_vars character vector giving names of variables with capture indicators
#' @param model character string giving capture-recapture Log-Linear model to estimate
#' @param hidden_condition character string giving condition to identify members of hidden population
#' @param label character string describing the estimator
#'
#' @return Data frame of HT estimates for single study
#' @export
#'
#' @references
#' Louis-Paul Rivest, Sophie Baillargeon. “The Rcapture package.” (2019). \url{https://cran.r-project.org/package=Rcapture}.
#'
#' @import dplyr
#' @importFrom Rcapture closedp.bc
get_study_est_recapture <- function(data,
                                    capture_vars = paste0("loc_", 1:4),
                                    model = "Mt",
                                    hidden_condition = "hidden == 1",
                                    label = "recapture") {

  .est_out <-
    data %>%
    dplyr::filter(if_any(all_of(capture_vars), ~ . == 1), eval(parse(text = hidden_condition))) %>%
    dplyr::select(all_of(capture_vars)) %>%
    Rcapture::closedp.bc(X = .,
                         dfreq = FALSE,
                         dtype = "hist",
                         t = length(capture_vars),
                         m = model)

  return(
    data.frame(estimator_label = paste0("hidden_size_", label),
               estimate = .est_out$results[1],
               se =  .est_out$results[2],
               inquiry_label = "hidden_size")
  )
}



