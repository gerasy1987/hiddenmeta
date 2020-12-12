#' SS-PSE population size estimator by Handcock, Gile and Mar
#'
#' @param data pass-through population data frame
#' @param prior_median prior median of hidden population size for SS-PSE estimation
#' @param rds_prefix character prefix used for RDS sample variables
#'
#' @return Data frame of SS-PSE estimates for a single study
#'
#' @references
#' Handcock, Mark S., Krista J. Gile, and Corinne M. Mar. “Estimating Hidden Population Size Using Respondent-Driven Sampling Data.” Electronic Journal of Statistics 8, no. 1 (2014): 1491–1521. \url{https://doi.org/10.1214/14-EJS923}.
#' @export
#'
#' @import dplyr
#' @importFrom sspse posteriorsize
get_study_est_sspse <- function(data, prior_median = 150, rds_prefix = "rds") {

  .quiet_sspse <- quietly(sspse::posteriorsize)

  .fit_sspse <-
    data %>%
    dplyr::filter_at(dplyr::vars(dplyr::all_of(rds_prefix)), ~ . == 1) %>%
    {
      .quiet_sspse(s = .$hidden_visible[order(.[,paste0(rds_prefix, "_t")])],
                   interval = 10,
                   median.prior.size = prior_median,
                   verbose = FALSE,
                   max.coupons = 3
      )
    }

  data.frame(estimator_label = c("hidden_size_sspse"),
             estimate = c(unname(.fit_sspse$result$N["Median AP"])),
             se =   c(sd(.fit_sspse$result$sample[,"N"])),
             estimand_label = c("hidden_size")
  )
}

#' Horvitz-Thompson prevalence estimator using weighted regression
#'
#' @param data pass-through population data frame
#' @param pps_prefix character prefix used for RDS sample variables
#'
#' @return Data frame of HT estimates for a single study
#' @export
#'
#' @import dplyr
#' @importFrom estimatr lm_robust
get_study_est_ht <- function(data, pps_prefix = "pps") {

  .fit_ht <-
    data %>%
    dplyr::filter_at(dplyr::vars(dplyr::all_of(pps_prefix)), ~ . == 1) %>%
    estimatr::lm_robust(hidden ~ 1, weights = 1/get(paste0(pps_prefix, "_share")), data = .) %>%
    {c(est = unname(.$coefficients), se = unname(.$std.error))}

  # quiet_HT <- quietly(sampling::HTestimator)
  # quiet_varHT <- quietly(sampling::varHT)
  #
  # fit_ht <-
  #   data %>%
  #   dplyr::filter_at(dplyr::vars(dplyr::all_of(pps_prefix)), ~ . == 1) %>%
  #   {
  #     c(est = quiet_HT(y = .$hidden_visible, pik = unlist(.[,paste0(pps_prefix, "_share")]))$result,
  #       se = sqrt(quiet_varHT(
  #         y = .$hidden_visible,
  #         pikl = sampling::UPmidzunopi2(unlist(.[,paste0(pps_prefix, "_prob")])), method = 1)$result))
  #
  #   }

 return(
   data.frame(estimator_label = "hidden_prev_ht",
              estimate = .fit_ht["est"],
              se =  .fit_ht["se"],
              estimand_label = "hidden_prev")
 )
}


#' Chords population size estimatior by Berchenko, Rosenblatt and Frost
#'
#' @param data pass-through population data frame
#' @param type a character vector with the type of estimation. Can be one of \code{mle}, \code{integrated}, \code{jeffreys} or \code{leave-d-out}. See \code{?chords::Estimate.b.k} and the original paper from the references for details
#' @param rds_prefix character prefix used for RDS sample variables
#'
#' @return Data frame of Chords estimates for a single study with RDS sample
#'
#' @references Berchenko, Yakir, Jonathan D. Rosenblatt, and Simon D. W. Frost. “Modeling and Analyzing Respondent-Driven Sampling as a Counting Process.” Biometrics 73, no. 4 (2017): 1189–98. \url{https://doi.org/10.1111/biom.12678}.
#'
#' @export
#'
#' @import dplyr
#' @importFrom chords initializeRdsObject Estimate.b.k makeJackControl
get_study_est_chords <- function(data,
                                 type = c("mle", "integrated", "jeffreys", "leave-d-out"),
                                 rds_prefix = "rds") {

  type <- match.arg(type)
  .K <- log2(length(grep(pattern = "^type_visible_", names(data))))
  .pattern <- paste0("^type_visible_", paste0(rep("[0-9]", .K - 1), collapse = ""), "1$")

  .data_mod <-
    data %>%
    dplyr::filter_at(dplyr::vars(dplyr::all_of(rds_prefix)), ~ . == 1) %>%
    dplyr::mutate(
      NS1 = apply(.[,grep(pattern = .pattern, x = names(data))], 1, sum),
      refCoupNum = get(paste0(rds_prefix, "_own_coupon")),
      interviewDt = get(paste0(rds_prefix, "_t"))) %>%
    dplyr::rename_at(vars(starts_with(paste0(rds_prefix, "_coupon_"))),
                     ~ gsub(pattern = paste0(rds_prefix, "\\_coupon\\_"), replacement = "coup", .))

  if (type == "leave-d-out") {
    .jack_control <- chords::makeJackControl(1, 1e2)

    .fit_chords <-
      chords::Estimate.b.k(rds.object = chords::initializeRdsObject(.data_mod),
                           type = type,
                           jack.control = .jack_control) %>%
      {c(est = sum(.$estimates$Nk.estimates[.$estimates$Nk.estimates < Inf]),
         degree_hidden =
           stats::weighted.mean(
             x = as.numeric(names(.$estimates$Nk.estimates))[.$estimates$Nk.estimates < Inf],
             w = .$estimates$Nk.estimates[.$estimates$Nk.estimates < Inf]))}

  } else {

    .fit_chords <-
      chords::Estimate.b.k(rds.object = chords::initializeRdsObject(.data_mod),
                           type = type) %>%
      {c(est = sum(.$estimates$Nk.estimates[.$estimates$Nk.estimates < Inf]),
         degree_hidden =
           stats::weighted.mean(
             x = as.numeric(names(.$estimates$Nk.estimates))[.$estimates$Nk.estimates < Inf],
             w = .$estimates$Nk.estimates[.$estimates$Nk.estimates < Inf]))}

  }

  return(
    data.frame(estimator_label = c("hidden_size_chords", "degree_hidden_chords"),
               estimate = c(.fit_chords["est"], .fit_chords["degree_hidden"]),
               estimand_label = c("hidden_size", "degree_hidden_average"))
  )
}


#' NSUM estimatior
#'
#' @param data pass-through population data frame
#' @param pps_prefix character prefix used for PPS sample variables
#' @param known_prefix character prefix for known population variables
#' @param hidden_prefix character prefix for hidden population variable
#'
#' @return Data frame of NSUM estimates for a single study with PPS sample
#' @export
#'
#' @references Dennis M. Feehan, Matthew J. Salganik. “The networkreporting package.” (2014). \url{http://cran.r-project.org/package=networkreporting}.
#'
#' @import dplyr
#' @importFrom networkreporting kp.degree.estimator nsum.estimator
get_study_est_nsum <- function(data, pps_prefix = "pps",
                               known_prefix = "known", hidden_prefix = "hidden_visible") {

  .quiet_nsum <- quietly(networkreporting::nsum.estimator)
  .quiet_degree_nsum <- quietly(networkreporting::kp.degree.estimator)

  .data_mod <-
    data %>%
    dplyr::filter_at(dplyr::vars(dplyr::all_of(pps_prefix)), ~ . == 1)

  .known_pops <-
    .data_mod %>%
    dplyr::select(total, matches(paste0("^total_", known_prefix))) %>%
    dplyr::rename_at(vars(starts_with("total_")), ~ paste0(gsub("^total\\_", "", .), "_visible")) %>%
    apply(., MARGIN = 2, unique)

  .data_mod$d_est <-
    .quiet_degree_nsum(survey.data = .data_mod,
                       known.popns = .known_pops[2:length(.known_pops)],
                       total.popn.size = .known_pops[1],
                       missing = "complete.obs")$result

  .fit_nsum <-
    .quiet_nsum(survey.data = .data_mod,
                d.hat.vals = .data_mod$d_est,
                total.popn.size = .known_pops[1],
                y.vals = hidden_prefix,
                missing = "complete.obs")$result

  # idu.est <-
  #   surveybootstrap::bootstrap.estimates(## this describes the sampling design of the
  #     ## survey; here, the PSUs are given by the
  #     ## variable cluster, and the strata are given
  #     ## by the variable region
  #     survey.design = ~ known + known_2 + known_3,
  #     ## the number of bootstrap resamples to obtain
  #     ## (NOTE: in practice, you should use more than 100.
  #     ##  this keeps building the package relatively fast)
  #     num.reps = 100,
  #     ## this is the name of the function
  #     ## we want to use to produce an estimate
  #     ## from each bootstrapped dataset
  #     estimator.fn = "nsum.estimator",
  #     bootstrap.fn = "rescaled.bootstrap.sample",
  #     ## our dataset
  #     survey.data = dat,
  #     ## other parameters we need to pass
  #     ## to the nsum.estimator function
  #     d.hat.vals = dat$d.hat,
  #     total.popn.size = known_pops[1],
  #     y.vals = "hidden_visible",
  #     missing = "complete.obs")

  return(
    data.frame(estimator_label = c("degree_average_nsum", "hidden_size_nsum"),
               estimate = c(unname(.fit_nsum$estimate), mean(.data_mod$d_est, na.rm = TRUE)),
               estimand_label = c("hidden_size", "degree_average"))
  )
}


#' Generalized NSUM estimatior
#'
#' @param data pass-through population data frame
#'
#' @return Data frame of HT estimates for single study
#' @export
#'
#' @import dplyr
#' @importFrom estimatr lm_robust
get_study_est_gnsum <- function(data) {

  return(
    data.frame(estimator_label = "hidden_prev_ht",
               estimate = fit_ht["est"],
               se =  fit_ht["se"],
               estimand_label = "hidden_prev")
  )
}


#' Service multiplier estimator
#'
#' @param data pass-through population data frame
#'
#' @return Data frame of HT estimates for single study
#' @export
#'
#' @import dplyr
#' @importFrom estimatr lm_robust
get_study_est_service <- function(data) {

  return(
    data.frame(estimator_label = "hidden_prev_ht",
               estimate = fit_ht["est"],
               se =  fit_ht["se"],
               estimand_label = "hidden_prev")
  )
}
