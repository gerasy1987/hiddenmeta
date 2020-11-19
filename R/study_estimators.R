#' Title
#'
#' @param data pass-through population data frame
#' @param prior_median prior median of hidden population size for SS-PSE estimation
#'
#' @return
#' @export
#'
#' @import dplyr
#' @importFrom sspse posteriorsize
get_study_est_sspse <- function(data, prior_median = 150) {

  quiet_sspse <- quietly(sspse::posteriorsize)

  fit_sspse <-
    data %>%
    dplyr::filter(rds == 1) %>%
    {
      quiet_sspse(s = .$hidden_visible[order(.$rds_t)],
                  interval = 10,
                  median.prior.size = prior_median,
                  verbose = FALSE,
                  max.coupons = 3
      )
    }


  data.frame(estimator_label = c("size_hidden_sspse"),
             estimate = c(unname(fit_sspse$result$N["Median AP"])),
             sd =   c(sd(fit_sspse$result$sample[,"N"])),
             estimand_label = c("size_hidden")
  )
}
