#' Read Study Designs from Google Spreadsheet
#'
#' @param ss Character string giving ID of the Google spreadsheet
#' @param sheet Character string giving name of the sheet in the Google spreadsheet to read
#'
#' @export
#'
#' @import tidyselect
#' @importFrom magrittr `%>%` `%$%`
#' @importFrom dplyr mutate filter select group_by ungroup summarize pull arrange rename_with left_join bind_rows if_all as_tibble
#' @importFrom plyr mapvalues
read_study_params <- function(
  ss = "1HwMM6JwoGALLMTpC8pQRVzaQdD2X71jpZNhfpKTnF8Y",
  sheet = "study_params_readable"
) {

  if (!requireNamespace("googlesheets4", quietly = TRUE)) {
    stop(
      "Package \"googlesheets4\" must be installed to use this function.",
      call. = FALSE
    )
  }

  labels <-
    googlesheets4::read_sheet(
      ss = ss,
      sheet = sheet,
      skip = 0,
      n_max = 4,
      col_names = FALSE,
      .name_repair = "minimal") %>%
    unname %>%
    t %>%
    `colnames<-`(c("family", "type", "name", "description")) %>%
    dplyr::as_tibble() %>%
    dplyr::mutate(
      # prior_type = purrr::map_chr(family, ~ strsplit(.x, "_")[[1]][2]),
      # family = purrr::map_chr(family, ~ strsplit(.x, "_")[[1]][1]),
      colnames = dplyr::if_else(is.na(family), name, paste0(family, "_", name))
    )


  spread_sheet <-
    googlesheets4::read_sheet(
      ss = ss,
      sheet = sheet,
      col_names = FALSE,
      na = "NA",
      col_types = paste0(labels$type, collapse = ""),
      .name_repair = "minimal",
      skip = 4) %>%
    `names<-`(dplyr::pull(labels, colnames))

  .out <-
    sapply(
      spread_sheet$study_id, simplify = FALSE, USE.NAMES = TRUE,
      FUN = function(study) {

        .row <- spread_sheet[which(spread_sheet$study_id == study),]

        if (.row$prior_known_hidden_interact > .row$observed_k_known)
          stop("More known groups that interact with hidden group than known groups considered")

        .group_names <-
          c(paste0("known", 1:(.row$prior_known_hidden_interact)),
            # this essentially precludes from having more than 1 hidden group
            paste0("hidden"))

        .pop <-
          list(
            handler = get_study_population,

            # network structure setup
            network_handler = sim_block_network,
            network_handler_args = list(
              N = .row$observed_n,
              K = .row$prior_known_hidden_interact + .row$observed_k_hidden,
              prev_K =
                `names<-`(
                  c(rep(.row$prior_known_prev, times = .row$prior_known_hidden_interact),
                    rep(.row$prior_hidden_prev, times = .row$observed_k_hidden)),
                  .group_names
                ),
              rho_K =
                if (length(.row$prior_rho) == 1) rep(.row$prior_rho, times = (length(.group_names)^2)/2 - length(.group_names)/2) else .row$prior_rho,
              p_edge_within =
                sapply(.group_names, simplify = FALSE,
                       function(x) {
                         if (grepl(pattern = "^known", x))
                           c(.row$prior_edge_within_known, .row$prior_edge_within_known)
                         else
                           c(.row$prior_edge_within_known, .row$prior_edge_within_hidden)}),
              p_edge_between =
                sapply(.group_names, simplify = FALSE,
                       function(x) {
                         if (grepl(pattern = "^known", x))
                           .row$prior_edge_within_known
                         else
                           .row$prior_edge_between}),
              directed = FALSE),

            # groups
            group_names = .group_names,

            # probability of visibility (show-up) for each group
            p_visible =
              `names<-`(as.list(
                c(rep(1, .row$prior_known_hidden_interact),
                  rep(.row$prior_p_visibility, .row$observed_k_hidden))),
                .group_names),

            # probability of service utilization in hidden population
            # for service multiplier
            add_groups = list()
          )

        # add service use where relevant (e.g. service multiplier is used)
        if (!is.na(.row$prior_p_service_use)) {

          .pop$add_groups <-
            c(
              .pop$add_groups,
              service_use = .row$prior_p_service_use
            )
        }

        # add time-locations presence if TLS sample is present
        if (.row$tls == 1) {
          # this only depends on hidden for now,
          # need to extend this to arbitrary number of hidden groups
          .pop$add_groups <-
            c(
              .pop$add_group,
              paste0(
                "paste0('loc_', 1:", .row$tls_n_time_locations,
                ") := lapply(rep(", .row$prior_p_showup,
                ", times = ", .row$tls_n_time_locations,
                "), function(add) rbinom(.N, 1, 0.05 + hidden * add))")
            )
        }

        if (.row$pps == 1) {
          if (!is.na(.row$pps_cluster)) {
            .pop$add_groups <-
              c(
                .pop$add_groups,
                cluster_for_pps =
                  paste0(
                    "sample(rep(1:ceiling(.N/", .row$pps_n_per_cluster ,
                    "), times = ", .row$pps_n_per_cluster , "))[1:.N]")
              )
          }
        }

        # add known groups that are not correlated with hidden group membership
        if (.row$prior_known_hidden_interact < .row$observed_k_known) {

          .group_names_add <-
            paste0("known", (.row$prior_known_hidden_interact + 1):.row$observed_k_known)

          .pop$add_groups <-
            c(
              .pop$add_group,
              `names<-`(
                as.list(rep(.row$prior_known_prev, times = length(.group_names_add))),
                .group_names_add
              )
            )
        }

        if (!is.na(.row$rds_allow_non_hidden) & .row$rds_allow_non_hidden == 1) {
          .pop$add_groups <-
            c(
              .pop$add_groups,
              `names<-`(as.list(.row$prior_hidden_prev), "hidden2")
            )
        }

        .inquiries <-
          list(handler = get_study_estimands,
               in_groups = if (!is.na(.row$rds_allow_non_hidden) & .row$rds_allow_non_hidden == 1) "^hidden2$" else NULL,
               out_groups = paste0("^known[1-", .row$prior_known_hidden_interact, "]$"),
               hidden_var = "hidden")

        .samples <-
          list(
            rds =
              if (.row$rds == 1)
                list(handler = sample_rds,
                     sampling_variable = "rds",
                     hidden_var = "hidden", # default
                     n_seed = .row$rds_n_seeds,
                     n_coupons = .row$rds_n_coupons,
                     n_waves =
                       if (!is.na(.row$rds_target_n_waves)) .row$rds_target_n_waves,
                     add_seeds = .row$rds_allow_add_seeds,
                     target_type = .row$rds_target_type,
                     target_n_rds = .row$rds_target_n_sample),
            tls =
              if (.row$tls == 1)
                list(handler = sample_tls,
                     sampling_variable = "tls",
                     hidden_var = if (.row$tls_allow_non_hidden == 0) "hidden" else NULL,
                     target_n_clusters = .row$tls_target_n_clusters,
                     target_cluster_type = .row$tls_target_cluster_type,
                     target_per_cluster = .row$tls_target_per_cluster,
                     clusters = paste0("loc_", 1:.row$tls_n_time_locations)),
            pps =
              if (.row$pps == 1)
                list(handler = sample_pps,
                     sampling_variable = "pps",
                     sampling_frame = NULL,
                     n_clusters =
                       if (!is.na(.row$pps_cluster)) .row$pps_cluster,
                     strata =
                       if (!is.na(.row$pps_strata)) "strata_for_pps",
                     cluster =
                       if (!is.na(.row$pps_cluster)) "cluster_for_pps",
                     target_n_pps = .row$pps_target_n_pps)
          ) %>%
          {.[sapply(., function(x) !is.null(x))]}

        .estimators <-
          list(
            rds =
              list(sspse = list(handler = get_study_est_sspse,
                                prior_mean = .row$prior_hidden_prev * .row$observed_n,
                                mcmc_params = list(interval = 5,
                                                   burnin = 2000,
                                                   samplesize = 1000),
                                total = .row$observed_n,
                                n_coupons = .row$rds_n_coupons,
                                prefix = "rds",
                                additional_params = list(alpha = 5),
                                label = "rds_sspse")
                   ),
            tls =
              list(ht = list(handler = get_study_est_ht,
                             weight_var = "tls_weight",
                             survey_design = ~ tls_loc_sampled,
                             prefix = "tls",
                             label = "tls_ht"),
                   nsum = list(handler = get_study_est_nsum,
                               known = paste0("known", 1:.row$observed_k_known),
                               hidden = "hidden_visible_out",
                               survey_design = ~ tls_loc_sampled,
                               n_boot = 100,
                               prefix = "tls",
                               label = "tls_nsum"),
                   recap = list(handler = get_study_est_recapture,
                                capture_parse =
                                  "strsplit(x = unique(na.omit(tls_locs_sampled)), split = ';')[[1]]",
                                sample_condition = "tls == 1",
                                model = "Mt",
                                hidden_variable = "hidden",
                                label = "tls_recap")),
            pps =
              list(ht = list(handler = get_study_est_ht,
                             prefix = "pps",
                             label = "pps_ht"),
                   nsum = list(handler = get_study_est_nsum,
                               known = paste0("known", 1:.row$observed_k_known),
                               hidden = "hidden_visible_out",
                               survey_design = ~ pps_cluster + strata(pps_strata),
                               n_boot = 100,
                               prefix = "pps",
                               label = "pps_nsum")),
            rds_pps =
              list(recap1 = list(handler = get_study_est_recapture,
                                 capture_vars = c("rds", "pps"),
                                 model = "Mt",
                                 hidden_variable = "hidden",
                                 label = "rds_pps_recap")),
            rds_tls =
              list(recap2 = list(handler = get_study_est_recapture,
                                 capture_vars = c("rds"),
                                 capture_parse =
                                   "strsplit(x = unique(na.omit(tls_locs_sampled)), split = ';')[[1]]",
                                 model = "Mt",
                                 hidden_variable = "hidden",
                                 label = "rds_tls_recap"))
          )

        # add service multiplier method where applicable
        if (!is.na(.row$multiplier_n_service_use)) {
          .estimators$rds$multiplier <-
            list(handler = get_study_est_multiplier,
                 service_var = "service_use",
                 total_service =
                   .row$prior_p_service_use * .row$observed_n * .row$prior_hidden_prev,
                 seed_condition = "rds_from == -999",
                 n_boot = 100,
                 prefix = "rds",
                 label = "rds_multi")
        }

        # add chords estimator for regular RDS samples
        if ((.row$rds == 1) & (.row$rds_seed_selection_type == "rds")) {
          .estimators$rds$chords <-
            list(handler = get_study_est_chords,
                 type = "jeffreys",
                 seed_condition = "rds_from == -999",
                 n_boot = 100,
                 prefix = "rds",
                 label = "rds_chords")
        }

        .estimators <-
          .estimators[sapply(strsplit(names(.estimators), "_"),
                             function(x) all(x %in% names(.samples)))]

        return(
          list(
            pop = .pop,
            inquiries = .inquiries,
            samples = .samples,
            estimators = .estimators
          )
        )
      })

  return(.out)
}


#' Generate Data Frame With Required Data by Study
#'
#' @param multi_study_data tibble with two columns, study (character identifier of the study) and population (list of population data frames/tibbles). This can be directly generated by multi-study design declaration
#' @param study_map named character vector of mapping between study labels and names used in the plot
#'
#' @export
#'
#' @import tidyselect
#' @importFrom magrittr `%>%` `%$%`
#' @importFrom dplyr mutate filter select group_by ungroup summarize pull arrange rename_with left_join bind_rows if_all as_tibble full_join rowwise
#' @importFrom plyr mapvalues
get_required_data <- function(
  multi_study_data,
  study_map = c("Brazil (FF)" = "ff_brazil",
                "Tunisia (UMass)" = "umass_tunisia",
                "Brazil (Stanford)" = "stanford_brazil",
                "Pakistan (JHU)" = "jhu_pakistan",
                "Costa Rica (John Jay)" = "johnjay_costarica",
                "Tanzania (John Jay)" = "johnjay_tanzania",
                "Morocco (NORC)" = "norc_marocco",
                "USA (RTI)" = "rti_usa")) {

  .out <-
    lapply(1:nrow(multi_study_data), function(i) {
      multi_study_data$population[[i]] %>%
        dplyr::select(
          name,
          starts_with(c("known", "hidden", "rds", "pps", "tls", "total"))) %>%
        {
          cbind(
            gsub(names(.), pattern = "[0-9]{1,2}", replacement = "XX"),
            sapply(unname(.), class),
            sapply(unname(.), function(x) {
              if (inherits(x, "numeric")) {
                paste0(format(sample(x = na.omit(x), size = 1), scientific = FALSE), collapse = ";")
              } else {
                paste0(na.omit(sample(x = x, size = 1)), collapse = ";")
              }
            }))
        } %>%
        {
          `names<-`(suppressMessages(
            as_tibble(., .name_repair = "universal")),
            c("variable", "type", "example"))
        } %>%
        dplyr::distinct(variable, type, .keep_all = TRUE) %>%
        dplyr::rename(
          "{ multi_study_data$study[i] }" := example
        )
    }) %>%
    purrr::reduce(., dplyr::full_join, by = c("variable", "type")) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      example = paste0(na.omit(unique(c_across(-c(variable,type)))), collapse = ";"),
      example = gsub("\\;{2,}", "", example)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      dplyr::across(all_of(multi_study_data$study), ~ ifelse(is.na(.), "", "X")),
      label =
        plyr::mapvalues(
          warn_missing = FALSE,
          variable,
          from = c(
            "name",
            "knownX", "knownX_visible_out", "knownX_visible_in",
            "hidden", "hidden_visible_out", "hidden_visible_in",

            "rds", "rds_from", "rds_t", "rds_wave", "rds_hidden",
            "rds_own_coupon", "rds_coupon_X",

            "pps", "pps_frame", "pps_weight", "pps_cluster", "pps_strata",

            "tls", "tls_loc_present", "tls_loc_sampled",
            "tls_weight",

            "total", "total_knownX", "total_service_use", "total_loc_X"),
          to = c("Respondent ID",
                 "Known group X member",
                 "Knows how many members of known group X",
                 "Known to be member of known group X to how many people",
                 "Hidden group member",
                 "Knows how many members of hidden group",
                 "Known to be member of hidden group to how many people",
                 "RDS: Sampled",
                 "RDS: ID of recruiting respondent",
                 "RDS: Time of recruitment",
                 "RDS: Wave",
                 "RDS: Member of hidden group",
                 "RDS: Own coupon ID",
                 "RDS: ID of recruitment coupon X",
                 "PPS: Sampled",
                 "PPS: Sampling frame indicator",
                 "PPS: Sampling weight",
                 "PPS: Cluster ID",
                 "PPS: Strata ID",
                 "TLS: Sampled",
                 "TLS: Time-locations that can be visited",
                 "TLS: Time-location at which sampled",
                 "TLS: Sampling weight (including probability of sampling location)",
                 "Total size of population (denominator for prevalence)",
                 "Total size of known groups",
                 "Total service use by hidden group members over time",
                 "Total size of time-locations included in sampling frame"
          ))
    ) %>%
    dplyr::select("Variable" = variable, "Label" = label, "Type" = type, "Example" = example, everything()) %>%
    dplyr::rename_with(~ plyr::mapvalues(., from = unname(study_map), to = names(study_map))) %>%
    dplyr::filter(Variable != Label)

  return(.out)

}


#' Create Diagnosis Plot of Log Normalized RMSE Across Studies
#'
#' @param diagnosands_df diagnosands data frame produced by \code{diagnose_design()}
#' @param study_map named character vector of mapping between study labels and names used in the plot
#' @param sampling_map named character vector of mapping between sampling strategies and names used in the plot
#' @param estimator_map named character vector of mapping between estimators and names used in the plot
#'
#' @export
#'
#' @import tidyselect ggplot2
#' @importFrom magrittr `%>%` `%$%`
#' @importFrom dplyr mutate filter select group_by ungroup summarize pull arrange rename_with left_join bind_rows if_all as_tibble
#' @importFrom purrr map2_chr map_chr
#' @importFrom tidyr nesting expand
#' @importFrom plyr mapvalues
get_rmse_plots <- function(
  diagnosands_df,
  study_map = c("Brazil (FF)" = "ff_brazil",
                "Brazil (Stanford)" = "stanford_brazil",
                "Costa Rica (John Jay)" = "johnjay_costarica",
                "Morocco (NORC)" = "norc_marocco",
                "Pakistan (JHU)" = "jhu_pakistan",
                "Tanzania (John Jay)" = "johnjay_tanzania",
                "Tunisia (UMass)" = "umass_tunisia",
                "USA (RTI)" = "rti_usa"),
  sampling_map = c("PPS" = "pps",
                   "RDS" = "rds",
                   "TLS" = "tls",
                   "RDS/TLS" = "rds_tls",
                   "RDS/PPS" = "rds_pps"),
  estimator_map = c("NSUM" = "nsum",
                    "HT" = "ht",
                    "Recapture" = "recap",
                    "SS-PSE" = "sspse",
                    "Service Multiplier" = "multi",
                    "Chords" = "chords")) {

  plot_df <-
    diagnosands_df %>%
    dplyr::filter(!is.na(estimator)) %>%
    dplyr::mutate(
      study = purrr::map_chr(inquiry, ~ sapply(.x, function(x) {
        unname(study_map)[sapply(unname(study_map), function(stud) grepl(stud, x, fixed = TRUE))]
      })),
      inquiry =
        purrr::map2_chr(study, inquiry, ~ gsub(pattern = paste0(.x, "_"), replacement = "", x = .y, fixed = TRUE)),
      estimator =
        purrr::map2_chr(inquiry, estimator, ~ gsub(pattern = paste0(.x, "_"), replacement = "", x = .y, fixed = TRUE)),
      target = rmse/mean_estimand,
      estimator_full = estimator,
      estimator = purrr::map_chr(estimator_full,
                                 ~ tail(strsplit(.x, "_")[[1]], 1)),
      sampling =
        purrr::map2_chr(estimator, estimator_full,
                        ~ gsub(pattern = paste0("_", .x), replacement = "", x = .y, fixed = TRUE)),
      study = plyr::mapvalues(study, from = unname(study_map), to = names(study_map),
                              warn_missing = FALSE),
      sampling = factor(
        levels = names(sampling_map),
        ordered = TRUE,
        x = plyr::mapvalues(sampling, from = unname(sampling_map), to = names(sampling_map),
                            warn_missing = FALSE)),
      estimator =
        plyr::mapvalues(estimator, from = unname(estimator_map), to = names(estimator_map),
                        warn_missing = FALSE)
    )

  if ("sim_type" %in% names(plot_df)) {

    suppressMessages(suppressWarnings(
      plot_df %>%
        dplyr::select(sim_type, study, inquiry, estimator_full, sampling, estimator, target) %>%
        {
          dplyr::left_join(
            x = tidyr::expand(., sim_type, study, tidyr::nesting(sampling, estimator)),
            y = .)
        } %>%
        dplyr::mutate(empty = ifelse(is.na(target), "X", NA_character_)) %>%
        ggplot2::ggplot(., aes(y = estimator, x = target, fill = sim_type)) +
        ggplot2::geom_col(width = .6, position = "dodge", color = "black", size = .1) +
        ggplot2::scale_x_continuous(trans = "sqrt") +
        ggplot2::facet_grid(sampling~study,
                            labeller = label_context,
                            scales = "free_y",
                            space = "free_y") +
        ggplot2::geom_text(aes(x = 0.1, label=empty), colour="red", size=4) +
        ggplot2::geom_vline(xintercept=0, color = "red", size=.5) +
        ggplot2::labs(x = "Normalized RMSE", y = "") +
        ggplot2::theme_bw() +
        ggplot2::theme(#axis.text.x = element_text(angle = 90),
                       axis.text.y = element_text(angle = 30),
                       legend.position = "bottom") +
        ggplot2::scale_fill_brewer(palette = "Set2") +
        ggplot2::guides(fill = guide_legend(title = "Simulation type",
                                            label.position = "right",
                                            keywidth = unit(.2, "cm")))
    ))

  } else {

    suppressMessages(suppressWarnings(
      plot_df %>%
        dplyr::select(study, inquiry, estimator_full, sampling, estimator, target) %>%
        {
          dplyr::left_join(
            x = tidyr::expand(., study, tidyr::nesting(sampling, estimator)),
            y = .)
        } %>%
        dplyr::mutate(empty = ifelse(is.na(target), "X", NA_character_)) %>%
        ggplot2::ggplot(., aes(y = estimator, x = target)) +
        ggplot2::geom_col(width = .2, position = "dodge", color = "black", size = .1) +
        ggplot2::scale_x_continuous(trans = "sqrt") +
        ggplot2::facet_grid(sampling~study,
                            labeller = label_context,
                            scales = "free_y",
                            space = "free_y") +
        ggplot2::geom_text(aes(x = 0.1, label = empty), colour="red", size=4) +
        ggplot2::geom_vline(xintercept = 0, color = "red", size=.5) +
        ggplot2::labs(x = "Normalized RMSE", y = "") +
        ggplot2::theme_bw() +
        ggplot2::theme(#axis.text.x = element_text(angle = 90),
                       axis.text.y = element_text(angle = 30))
    ))

  }

}
