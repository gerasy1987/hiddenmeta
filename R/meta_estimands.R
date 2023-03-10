#' Get meta-analysis estimands for specific study-level estimand
#'
#' @param data pass-through meta population
#' @param study_estimand study-level estimand for meta-analysis estimands
#' @param samp_est_benchmark within-study benchmark sampling-estimator strategy
#'
#' @return Data frame of specified meta level estimands
#'
#' @export
#'
#' @import tidyselect
#' @importFrom magrittr `%>%` `%$%`
#' @importFrom dplyr mutate filter select group_by ungroup summarize pull arrange rename_with left_join bind_rows if_all if_any
#' @importFrom tidyr pivot_longer
#' @importFrom purrr pmap_chr
get_meta_estimands <-
  function(data,
           study_estimand = "hidden_size",
           samp_est_benchmark = "pps_ht") {

    est <-
      data %>%
      dplyr::filter(inquiry == study_estimand) %>%
      {
        dplyr::bind_rows(
          # study level estimands
          dplyr::group_by(., sim_ID, study) %>%
            dplyr::summarise(.groups = "drop",
                             estimand = unique(estimand)) %>%
            dplyr::mutate(inquiry = paste0(study, "_", study_estimand)) %>%
            dplyr::select(sim_ID, inquiry = inquiry, estimand),
          # # sample-estimator level absolute estimands
          # dplyr::group_by(., sim_ID, sample, estimator) %>%
          #   dplyr::summarise(.groups = "drop",
          #                    bias = mean(estimate - estimand),
          #                    bias_sd = sd(estimate - estimand)) %>%
          #   tidyr::pivot_longer(cols = c(bias, bias_sd),
          #                       values_to = "estimand", names_to = "inquiry") %>%
          #   dplyr::arrange(inquiry) %>%
          #   dplyr::mutate(
          #     inquiry =
          #       purrr::pmap_chr(
          #         list(x = sample, y = inquiry, z = estimator),
          #         function(x,y,z) paste0(y, "_", paste0(x, collapse = "_"), "_", z)
          #       )) %>%
          #   dplyr::select(sim_ID, inquiry, estimand),
          # sample-estimator level ratio estimands
          dplyr::group_by(., sim_ID, sample, estimator) %>%
            dplyr::summarise(.groups = "drop",
                             rel_bias = mean(estimate/estimand)#,
                             # rel_bias_sd = sd(estimate/estimand)
                             ) %>%
            tidyr::pivot_longer(cols = c(rel_bias #,
                                         # rel_bias_sd
                                         ),
                                values_to = "estimand", names_to = "inquiry") %>%
            dplyr::mutate(
              samp_est =
                purrr::pmap_chr(
                  list(x = sample, y = estimator),
                  function(x,y) paste0(paste0(x, collapse = "_"), "_", y)
                )) %>%
            dplyr::select(sim_ID, samp_est, inquiry, estimand) %>%
            dplyr::group_by(sim_ID, inquiry) %>%
            dplyr::mutate(estimand = estimand/estimand[which(samp_est == samp_est_benchmark)]) %>%
            dplyr::ungroup() %>%
            dplyr::arrange(inquiry) %>%
            dplyr::mutate(inquiry = paste0(inquiry, "_", samp_est)) %>%
            dplyr::select(sim_ID, inquiry, estimand)
        )
      } %>%
      dplyr::group_by(inquiry) %>%
      dplyr::summarise(estimand = mean(estimand, na.rm = TRUE)) %>%
      as.data.frame(stringsAsFactors = FALSE)

    return(est)

  }
