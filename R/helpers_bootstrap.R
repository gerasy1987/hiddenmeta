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
#' @param data sample data frame
#' @param seed_condition character string containing condition to define seeds. Defaults to "rds_from == -999" that applies to simulated RDS samples
#' @param in_coupon character string for variable containing in coupons identifiers. In coupon is coupon used to enroll specific respondent
#' @param out_coupon character string for variable(s) containing out coupons identifiers. Out coupons are coupons given to respondent to enroll other study participants from their network
#' @param is_child_fn function used to identify direct childs of specific seed. Defaults to function that uses "-" in the coupon identifiers
#' @param trait_var variable used to split RDS sample for resampling. Defaults to "hidden_visible_out" that captures number of hidden population members recognized by respondent
#' @param other_vars other variables to preserve in bootstrap samples
#' @param n_boot number of bootstrap resamples
#'
#' @export
#'
#' @references Dennis M. Feehan, Matthew J. Salganik. “The surveybootstrap package.” (2016). \url{https://cran.r-project.org/package=surveybootstrap}.
#'
#' @import tidyselect surveybootstrap
#' @importFrom magrittr `%>%` `%$%`
#' @importFrom dplyr mutate filter select group_by ungroup summarize pull arrange rename_with left_join bind_rows if_all as_tibble
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
      # dplyr::filter(get(prefix) == 1) %>%
      dplyr::mutate(
        across(all_of(trait_var), list("b" = ~ as.integer(. > median(., na.rm = TRUE))))) %>%
      `attr<-`(., "key", in_coupon)

    seed_ids <-
      rds_data %>%
      dplyr::filter(eval(parse(text = seed_condition))) %>%
      dplyr::pull(in_coupon)

    # rds_chains <-
    #   plyr::llply(seed_ids,
    #               make_chain,
    #               key = in_coupon,
    #               is_child_fn = is_child_fn,
    #               data = rds_data)

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


#' Modified Function for Rescaled Bootstrap by Rust and Rao (1996)
#'
#' @param data sample data frame
#' @param survey_design a formula describing the design of the survey
#' @param n_boot number of bootstrap resamples
#'
#' @export
#'
#' @references
#' Dennis M. Feehan, Matthew J. Salganik. “The surveybootstrap package.” (2016). \url{https://cran.r-project.org/package=surveybootstrap}.
#' Rust, Keith F., and J. N. K. Rao. "Variance estimation for complex surveys using replication techniques." Statistical methods in medical research 5, no. 3 (1996): 283-310.
#'
#' @import surveybootstrap
#' @importFrom dplyr group_indices group_by across all_of group_split
#' @importFrom plyr llply alply
get_rescaled_boot <-
  function(data,
           survey_design,
           n_boot = 10) {

    data$.internal_id <- 1:nrow(data)

    design <- extract_terms(survey_design)

    if (length(design$psu) == 1 & design$psu == "1") {
      design$psu <- as.name(".internal_id")
    }

    data$.cluster_id <-
      dplyr::group_indices(dplyr::group_by(data, dplyr::across(dplyr::all_of(design$psu))))

    if (is.null(design$strata)) {
      strata <- list(data)
    } else {
      strata <- dplyr::group_split(dplyr::group_by(data, dplyr::across(dplyr::all_of(design$strata))))
    }

    boot_weights <- plyr::llply(strata, function(stratum_data) {
      res <- surveybootstrap:::resample_stratum(stratum_data$.cluster_id, n_boot)
      colnames(res) <- paste0("rep.", 1:ncol(res))
      res <- cbind(index = stratum_data$.internal_id, res)
      return(res)
    })

    boot_weights <- do.call("rbind", boot_weights)

    res <- plyr::alply(boot_weights[, -1, drop = FALSE], 2, function(this_col) {
      return(data.frame(index = boot_weights[, 1], weight.scale = this_col))
    })

    return(res)
  }

# report.aggregator_ <-
#   function(resp.data, attribute.names, qoi, weights, qoi.name,
#            dropmiss = FALSE) {
#     resp.data <- resp.data
#     wdat <- select_(resp.data, .dots = weights)
#     qdat <- select_(resp.data, .dots = qoi)
#     adat <- select_(resp.data, .dots = attribute.names)
#     df <- bind_cols(wdat, qdat, adat)
#     wgt.col <- as.symbol(names(df)[1])
#     qoi.col <- as.symbol(names(df)[2])
#     grouping.cols <- names(df)[-1]
#     dots <- lapply(grouping.cols, as.symbol)
#     df.summ <-
#       df %>%
#       group_by_(.dots = dots) %>%
#       dplyr::summarise_(mean.qoi = lazyeval::interp(~weighted.mean(a, w = b), a = qoi.col, b = wgt.col),
#                         sum.qoi = lazyeval::interp(~sum(a * b), a = qoi.col, b = wgt.col),
#                         wgt.total = lazyeval::interp(~sum(b), b = wgt.col),
#                         wgt.inv.total = lazyeval::interp(~sum(1/b), b = wgt.col),
#                         num.obs = lazyeval::interp(~length(a), a = qoi.col))
#     toren <- list(~mean.qoi, ~sum.qoi, ~wgt.total, ~wgt.inv.total,
#                   ~num.obs)
#     newnames <- paste0(c("mean.", "sum.", "wgt.total.", "wgt.inv.total.",
#                          "num.obs."), qoi.name)
#     df.summ <- dplyr::rename_(df.summ, .dots = setNames(toren,
#                                                         newnames))
#     return(df.summ)
#   }
#
# kp.estimator <-
#   function(resp.data, known.populations, attribute.names, weights,
#            total.kp.size = NULL, alter.popn.size = NULL) {
#     wdat <- select_(resp.data, .dots = weights)
#     kpdat <- select_(resp.data, .dots = known.populations)
#     adat <- select_(resp.data, .dots = attribute.names)
#     alter.popn.size <- ifelse(is.null(alter.popn.size) || is.null(lazyeval::lazy_eval(alter.popn.size)),
#                               sum(wdat[, 1]), lazyeval::lazy_eval(alter.popn.size))
#     total.kp.size <- ifelse(is.null(total.kp.size), 1, lazyeval::lazy_eval(total.kp.size))
#     kptot <- data_frame(kptot = rowSums(kpdat))
#     df <- bind_cols(kptot, wdat, adat)
#     atnames <- paste(colnames(adat))
#     agg <- report.aggregator_(resp.data = df, attribute.names = atnames,
#                               qoi = "kptot", weights = weights, qoi.name = "y.kp")
#     tograb <- lapply(c(colnames(adat), "sum.y.kp", "wgt.total.y.kp",
#                        "num.obs.y.kp"), as.symbol)
#     sum.y.kp <- NULL
#     sum.y.kp.over.kptot <- NULL
#     wgt.total.y.kp <- NULL
#     res <-
#       select_(agg, .dots = tograb) %>%
#       dplyr::mutate(sum.y.kp.over.kptot = sum.y.kp/total.kp.size,
#                     dbar.Fcell.F = sum.y.kp.over.kptot * (alter.popn.size/wgt.total.y.kp))
#     return(res)
#   }
