#' SS-PSE estimator
#'
#' @param data pass-through population data frame
#' @param prior_median prior median of hidden population size for SS-PSE estimation
#' @param rds_prefix character prefix used for RDS sample variable
#'
#' @return Data frame of SS-PSE estimates for single study
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

#' Horvitz-Thompson estimatior
#'
#' @param data pass-through population data frame
#' @param prop_prefix character prefix used for RDS sample variable
#'
#' @return Data frame of HT estimates for single study
#' @export
#'
#' @import dplyr
#' @importFrom estimatr lm_robust
get_study_est_ht <- function(data, prop_prefix = "prop") {

  fit_ht <-
    data %>%
    dplyr::filter_at(dplyr::vars(dplyr::all_of(prop_prefix)), ~ . == 1) %>%
    estimatr::lm_robust(hidden ~ 1, weights = 1/prop_share, data = .) %>%
    {c(est = unname(.$coefficients), se = unname(.$std.error))}

  # quiet_HT <- quietly(sampling::HTestimator)
  # quiet_varHT <- quietly(sampling::varHT)
  #
  # fit_ht <-
  #   data %>%
  #   dplyr::filter_at(dplyr::vars(dplyr::all_of(prop_prefix)), ~ . == 1) %>%
  #   {
  #     c(est = quiet_HT(y = .$hidden_visible, pik = unlist(.[,paste0(prop_prefix, "_share")]))$result,
  #       se = sqrt(quiet_varHT(
  #         y = .$hidden_visible,
  #         pikl = sampling::UPmidzunopi2(unlist(.[,paste0(prop_prefix, "_prob")])), method = 1)$result))
  #
  #   }

 return(
   data.frame(estimator_label = "hidden_prev_ht",
              estimate = fit_ht["est"],
              se =  fit_ht["se"],
              estimand_label = "hidden_prev")
 )
}