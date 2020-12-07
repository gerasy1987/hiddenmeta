#' SS-PSE estimator by Handcock, Gile and Mar
#'
#' @param data pass-through population data frame
#' @param prior_median prior median of hidden population size for SS-PSE estimation
#' @param rds_prefix character prefix used for RDS sample variable
#'
#' @return Data frame of SS-PSE estimates for single study
#'
#' @references
#' Handcock, Mark S., Krista J. Gile, and Corinne M. Mar. “Estimating Hidden Population Size Using Respondent-Driven Sampling Data.” Electronic Journal of Statistics 8, no. 1 (2014): 1491–1521. \url{https://doi.org/10.1214/14-EJS923}.
#' @export
#'
#' @import dplyr
#' @importFrom sspse posteriorsize
get_study_est_sspse <- function(data, prior_median = 150, rds_prefix = "rds") {



  quiet_sspse <- quietly(sspse::posteriorsize)

  fit_sspse <-
    data %>%
    dplyr::filter_at(dplyr::vars(dplyr::all_of(rds_prefix)), ~ . == 1) %>%
    {
      quiet_sspse(s = .$hidden_visible[order(.[,paste0(rds_prefix, "_t")])],
                  interval = 10,
                  median.prior.size = prior_median,
                  verbose = FALSE,
                  max.coupons = 3
      )
    }


  data.frame(estimator_label = c("hidden_size_sspse"),
             estimate = c(unname(fit_sspse$result$N["Median AP"])),
             se =   c(sd(fit_sspse$result$sample[,"N"])),
             estimand_label = c("hidden_size")
  )
}

#' Horvitz-Thompson estimatior using weighted regression
#'
#' @param data pass-through population data frame
#' @param pps_prefix character prefix used for RDS sample variable
#'
#' @return Data frame of HT estimates for single study
#' @export
#'
#' @import dplyr
#' @importFrom estimatr lm_robust
get_study_est_ht <- function(data, pps_prefix = "pps") {

  fit_ht <-
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
              estimate = fit_ht["est"],
              se =  fit_ht["se"],
              estimand_label = "hidden_prev")
 )
}


#' Chords estimatior
#'
#' @param data pass-through population data frame
#'
#' @return Data frame of HT estimates for single study
#'
#' @references Berchenko, Yakir, Jonathan D. Rosenblatt, and Simon D. W. Frost. “Modeling and Analyzing Respondent-Driven Sampling as a Counting Process.” Biometrics 73, no. 4 (2017): 1189–98. \url{https://doi.org/10.1111/biom.12678}.
#'
#' @export
#'
#' @import dplyr
#' @importFrom estimatr lm_robust
get_study_est_chords <- function(data) {

  return(
    data.frame(estimator_label = "hidden_prev_ht",
               estimate = fit_ht["est"],
               se =  fit_ht["se"],
               estimand_label = "hidden_prev")
  )
}


#' NSUM estimatior
#'
#' @param data pass-through population data frame
#'
#' @return Data frame of HT estimates for single study
#' @export
#'
#' @import dplyr
#' @importFrom estimatr lm_robust
get_study_est_nsum <- function(data) {


  return(
    data.frame(estimator_label = "hidden_prev_ht",
               estimate = fit_ht["est"],
               se =  fit_ht["se"],
               estimand_label = "hidden_prev")
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
