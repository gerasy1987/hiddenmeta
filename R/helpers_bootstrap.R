# function detecting whether RDS respondent received coupon from seed_id (directly)
is_child <- function(ids, seed_id) {
  res <- strsplit(ids, split = "-")
  return(
    plyr::laply(
      .data = res,
      .fun = function(x) {
        length(x) != 1 &&
          paste0(x[1:(length(x) - 1)], collapse = "-") == seed_id
      })
  )
}

# recursive (!!!) function to make chain data suitable for surveybootstrap function
make_chain <- function(seed_id, data, key, is_child_fn = is_child) {

  keys <- data[[key]]
  seed_data <- data[keys == seed_id, ]
  child_ids <- keys[is_child_fn(keys, seed_id)]

  if (length(child_ids) == 0) {
    return(
      list(data = seed_data,
           children = NULL)
    )
  } else {
    return(
      list(data = seed_data,
           children = plyr::llply(child_ids,
                                  make_chain,
                                  data = data,
                                  key = key,
                                  is_child_fn = is_child_fn))
    )
  }

}

#' Prepare RDS data for Bootstrap Re-sampling Following Salganik (2006)
#'
#' @param data pass-through population data frame
#' @param seed_condition character string containing condition to define seeds. Defaults to "rds_from == -999" that applies to simulated RDS samples
#' @param in_coupon character string for variable containing in coupons identifiers. In coupon is coupon used to enroll specific respondent
#' @param out_coupon character string for variable(s) containing out coupons identifiers. Out coupons are coupons given to respondent to enroll other study participants from their network
#' @param is_child_fn function used to identify direct childs of specific seed. Defaults to function that uses "-" in the coupon identifiers
#' @param trait_var variable used to split RDS sample for resampling. Defaults to "hidden_visible_out" that captures number of hidden population members recognized by respondent
#' @param other_vars other variables to preserve in bootstrap samples
#' @param n_boot number of bootstrap resamples
#'
#' @return
#' @export
#'
#' @examples
#'
#' @references Dennis M. Feehan, Matthew J. Salganik. “The surveybootstrap package.” (2016). \url{https://CRAN.R-project.org/package=surveybootstrap}.
#'
#' @import dplyr
#' @import surveybootstrap
#' @importFrom tidyr pivot_longer
get_rds_boot <-
  function(data,
           seed_condition = "rds_from == -999",
           in_coupon = "rds_own_coupon",
           out_coupon = "rds_coupon_",
           is_child_fn = is_child,
           trait_var = "hidden_visible_out",
           other_vars = c("hidden_visible_out",
                          "known_visible_out",
                          "hidden"),
           n_boot = 100) {

    rds_data <-
      data %>%
      # dplyr::filter(get(rds_prefix) == 1) %>%
      dplyr::mutate(
        across(all_of(trait_var), list("b" = ~ as.integer(. > median(., na.rm = TRUE))))) %>%
      `attr<-`(., "key", in_coupon)

    seed_ids <-
      rds_data %>%
      dplyr::filter(eval(parse(text = seed_condition))) %>%
      dplyr::pull(in_coupon)

    rds_chains <-
      plyr::llply(seed_ids,
                  make_chain,
                  key = in_coupon,
                  is_child_fn = is_child_fn,
                  data = rds_data)

    rds_data_nonseed <-
      rds_data %>%
      dplyr::filter(!(get(in_coupon) %in% seed_ids)) %>%
      as.data.frame()

    rds_data_parent <-
      rds_data %>%
      tidyr::pivot_longer(
        cols = names(rds_data)[grep(x = names(rds_data), pattern = out_coupon)],
        names_to = "coupon_id",
        values_to = "out_coupon"
      ) %>%
      `attr<-`(., "key", "out_coupon") %>%
      as.data.frame()

    rds_boot_data <-
      list(chains =
             plyr::llply(
               .data = seed_ids,
               .fun = make_chain,
               key = in_coupon,
               is_child_fn = is_child_fn,
               data = rds_data),
           mm =
             surveybootstrap:::estimate.mixing(
               survey.data = rds_data_nonseed,
               parent.data = rds_data_parent,
               traits = paste0(trait_var, "_b")),
           dd =
             surveybootstrap:::estimate.degree.distns(
               survey.data = rds_data,
               d.hat.vals = trait_var,
               traits = paste0(trait_var, "_b"),
               keep.vars = other_vars))

    return(
      surveybootstrap:::rds.mc.boot.draws(rds_boot_data$chains,
                        rds_boot_data$mm,
                        rds_boot_data$dd,
                        num.reps = n_boot)
    )

  }
