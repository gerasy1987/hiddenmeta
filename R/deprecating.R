#' #' Simulate single population with given network structure using tidyverse
#' #'
#' #' @param network_handler function that takes several arguments and returns igraph network object with vertex attribute \code{type} in the binary coded format (consists of 0's and 1's only)
#' #' @param network_handler_args list of arguments passed to \code{network_handler}
#' #' @param group_names character vector of groups within population with last group being hidden
#' #' @param p_visible named list of visibility tendencies by group. This is used as mean of Beta distribution (with SD = 0.09) to generate probability of being recognized as member of group, being sampled as seed, etc. The order of objects in list have to follow the order of \code{group_names}
#' #' @param add_groups named list of probabilities of additional group memberships. Examples include probability of service utilization (for service multiplier), being present at particular time-location (for TLS), etc.
#' #'
#' #' @return Population data frame for single study
#' #'
#' #' @keywords internal
#' #'
#' #'
#' #' @examples
#' #' \dontrun{
#' #'   get_study_population_tidy(
#' #'     network_handler = sim_block_network,
#' #'     network_handler_args =
#' #'       list(N = 1000, K = 2, prev_K = c(known = .1, hidden = .2), rho_K = 0,
#' #'            p_edge_within = list(known = c(0.1, 0.1), hidden = c(0.1, 0.3)),
#' #'            p_edge_between = list(known = 0.02, hidden = 0.02),
#' #'            directed = FALSE),
#' #'
#' #'     # groups
#' #'     group_names = c("known", "hidden"),
#' #'
#' #'     # probability of visibility (show-up) for each group
#' #'     p_visible = list(known = 1, hidden = 1),
#' #'
#' #'     # probability of service utilization in hidden population
#' #'     # for service multiplier
#' #'     add_groups =
#' #'       list(
#' #'         service_use = "rbinom(n(), 1, 0.25)",
#' #'         "purrr::map_df(hidden, ~ sapply( `names<-`(rep(0.2, times = 10),
#' #'         paste0('loc_', 1:10)), function(add) rbinom(length(.x), 1, 0.05 + .x * add)))",
#' #'         known_2 = 0.3,
#' #'         "purrr::map_df(known, ~ sapply( `names<-`(rep(0, times = 8),
#' #'         paste0('known_', 3:10)), function(add) rbinom(length(.x), 1, 0.3)))")
#' #'   )
#' #' }
#' #'
#' get_study_population_tidy <-
#'   function(network_handler = sim_block_network,
#'            network_handler_args,
#'            group_names,
#'            p_visible,
#'            add_groups) {
#'
#'     if (length(unique(group_names)) != length(group_names) |
#'         any(names(p_visible) != group_names))
#'       stop("group names have to be unique and consistent")
#'
#'     if (any(grepl("visible|type|_in|_out", c(names(add_groups), group_names))))
#'       stop("some group names include `visible`, `type`, `_in` or `_out`.
#'            This might clash with internal naming conventions and produce inadequate results")
#'
#'     g <- do.call(what = network_handler, args = network_handler_args)
#'
#'     .g_adj <- igraph::as_adj(g)
#'     .g_attr <- igraph::vertex_attr(g)
#'
#'     if (!("igraph" %in% class(g)))
#'       stop("network_handler should produce an igraph object")
#'
#'     .base_df <-
#'       .g_attr$type %>%
#'       stringr::str_split(., pattern = "", simplify = TRUE) %>%
#'       dplyr::as_tibble(.name_repair = function(names) group_names) %>%
#'       dplyr::mutate(
#'         across(all_of(group_names), as.integer),
#'         type = .g_attr$type)
#'
#'     .out_df <-
#'       .base_df %>%
#'       mutate_add_groups_tidy(., add_groups = add_groups)
#'
#'     .add_groups_names <- setdiff(names(.out_df), names(.base_df))
#'
#'     .out_df %>%
#'       mutate_visibility_tidy(., vars = names(p_visible), p_visible = p_visible) %>%
#'       dplyr::mutate(
#'         type_visible =
#'           apply(X = .[,paste0(group_names, "_visible")], MARGIN = 1, FUN = paste0, collapse = "")) %>%
#'       fastDummies::dummy_cols(., select_columns = "type_visible", remove_selected_columns = TRUE) %>%
#'       dplyr::rename_with(
#'         .cols = starts_with("type_visible_"),
#'         ~ sapply(base::strsplit(.x, split = "_"), function(x) paste0(x[c(1,3,2)], collapse = "_"))) %>%
#'       {
#'         dplyr::bind_cols(
#'           .,
#'           # out-reports given visibility of others
#'           `names<-`(
#'             as.data.frame(
#'               as.matrix(
#'                 .g_adj %*% as.matrix(.[,names(.)[grepl("_visible$", names(.))]]))),
#'             paste0(names(.)[grepl("_visible$", names(.))], "_out")),
#'           `names<-`(
#'             as.data.frame(
#'               as.matrix(
#'                 .g_adj %*% as.matrix(.[,.add_groups_names]))),
#'             paste0(.add_groups_names, "_visible_out")),
#'           # in-reports given subject's visibility
#'           `names<-`(
#'             as.data.frame(
#'               as.matrix(
#'                 Matrix::rowSums(.g_adj) * as.matrix(.[,names(.)[grepl("_visible$", names(.))]]))),
#'             paste0(names(.)[grepl("_visible$", names(.))], "_in")),
#'           `names<-`(
#'             as.data.frame(
#'               as.matrix(
#'                 Matrix::rowSums(.g_adj) * as.matrix(.[,.add_groups_names]))),
#'             paste0(.add_groups_names, "_visible_in"))
#'         )
#'       } %>%
#'       dplyr::mutate(
#'         n_visible_out = apply(X = .[,grep("^type_.*_out$", names(.))], MARGIN = 1, FUN = sum),
#'         name = 1:n(),
#'         links = sapply(igraph::as_adj_list(g), paste0, collapse = ";"),
#'         total = n(),
#'         across(all_of(c(group_names, .add_groups_names)),
#'                list(total = sum), .names = "{.fn}_{.col}")) %>%
#'       dplyr::bind_cols(.,
#'                        .g_attr[-which(names(.g_attr) == "type")]) %>%
#'       dplyr::select(name, type, all_of(group_names), links, # main parameters
#'                     .add_groups_names, # additional group memberships
#'                     n_visible_out,
#'                     ends_with("_visible_out"), # out-reports
#'                     ends_with("_visible_in"), # in-reports
#'                     starts_with("p_visible_"), # probabilities of visible
#'                     starts_with("total")) # total numbers by group
#'   }
#'
#'
#' #' Draw probability sample proportional to size (PPS) from single study
#' #'
#' #' Sampling handler for drawing proportional sample with given characteristics from individual study population
#' #'
#' #' @param data pass-through population data frame
#' #' @param sample_label character string that is used as prefix for all variables generated by proportional sampling. Default is 'pps'
#' #' @param drop_nonsampled logical indicating whether to drop units that are not sampled. Default is \code{FALSE}
#' #' @param target_n_pps target size of proportional sample
#' #' @param n_clusters number of clusters
#' #' @param sampling_frame character vector containing all binary vectors identifying sampling frame
#' #' @param cluster character vector containing name(s) of all cluster variables
#' #' @param strata character vector containing name(s) of stratifying variables. Currently not implemented
#' #' @param weights_type character string giving the type of weights to compute. Can be one of "absolute" or "relative". Currently only absolute weights are calculated
#' #'
#' #' @return Population or sample data frame for single study with PPS sample characteristics added
#' #'  \describe{
#' #'   \item{[sample_label]}{Sampling indicator}
#' #'   \item{[sample_label]_frame}{Indicator for sampling frame (respondents with 0 cannot be enrolled)}
#' #'   \item{[sample_label]_strata}{ID of respondent's strata}
#' #'   \item{[sample_label]_strata_prop}{Proportion of sampling frame in the strata to which respondent belongs}
#' #'   \item{[sample_label]_cluster}{ID of respondent's cluster}
#' #'   \item{[sample_label]_cluster_prop}{Proportion of strata in the cluster to which respondent belongs}
#' #'   \item{[sample_label]_weight}{Sampling weights}
#' #'  }
#' #' @keywords internal
#' #'
#' sample_pps_tidy <-
#'   function(data, sample_label = "pps", drop_nonsampled = FALSE,
#'            target_n_pps = 400,
#'            n_clusters = target_n_pps,
#'            sampling_frame = NULL,
#'            strata = NULL,
#'            cluster = NULL,
#'            weights_type = c("absolute", "relative")
#'   ) {
#'
#'     warning("This function was deprecated due to lacking performance. It has been replaced with an optimized version. Consider using sample_pps().")
#'     weights_type <- match.arg(weights_type)
#'
#'     data <- dplyr::mutate(data, temp_id = 1:n())
#'     temp_data <- data %>% dplyr::select(temp_id, all_of(c(sampling_frame, strata, cluster)))
#'     nmax_per_cluster <- ceiling(target_n_pps/n_clusters)
#'
#'
#'     if (!is.null(sampling_frame)) {
#'       temp_data %<>%
#'         dplyr::mutate(frame = apply(.[,..sampling_frame, drop = FALSE], 1, function(x) as.integer(all(x == 1))))
#'     } else {
#'       temp_data %<>% dplyr::mutate(frame = 1)
#'     }
#'
#'     # check that sampling is possible
#'     if (target_n_pps >= sum(temp_data$frame)) {
#'       stop("Requested PPS sample size is larger than population (or sampling frame) size")
#'     }
#'
#'     if (!is.null(strata)) {
#'       temp_data %<>%
#'         tidyr::nest(dat = -dplyr::all_of(strata)) %>%
#'         dplyr::mutate(strata_id = 1:n(),
#'                       strata_prop = purrr::map_int(dat, ~ nrow(.x)),
#'                       strata_prop = strata_prop/sum(strata_prop))%>%
#'         tidyr::unnest(cols = c(dat))
#'     } else {
#'       temp_data %<>%
#'         dplyr::mutate(strata_id = 1,
#'                       strata_prop = 1)
#'     }
#'
#'     nclustmax_per_strat <- ceiling(n_clusters/length(unique(temp_data$strata_id)))
#'
#'     temp_data %<>%
#'       split(., .$strata_id) %>%
#'       plyr::ldply(
#'         .data = .,
#'         .fun = function(strat_df) {
#'           if (!is.null(cluster)) {
#'             strat_df %<>%
#'               tidyr::nest(dat = -dplyr::all_of(cluster)) %>%
#'               dplyr::mutate(cluster_id = 1:n(),
#'                             cluster_prop = purrr::map_int(dat, ~ nrow(.x)),
#'                             cluster_prop = cluster_prop/sum(cluster_prop)) %>%
#'               tidyr::unnest(cols = c(dat))
#'
#'             strat_df %>%
#'               dplyr::filter(frame == 1) %>%
#'               tidyr::nest(dat = -c(cluster_id, cluster_prop)) %>%
#'               dplyr::mutate(
#'                 sampled_cluster =
#'                   ifelse(cluster_id %in% sample(x = cluster_id,
#'                                                 size = min(n(), nclustmax_per_strat),
#'                                                 prob = cluster_prop), 1, 0),
#'                 dat =
#'                   purrr::map2(
#'                     dat, cluster_prop,
#'                     ~ dplyr::mutate(.x,
#'                                     sampled = sample(c(rep(1,min(n(), nmax_per_cluster)),
#'                                                        rep(0,max(0, n() - nmax_per_cluster)))),
#'                                     weight = n()/(sum(sampled) * .y)))
#'               ) %>%
#'               tidyr::unnest(cols = c(dat)) %>%
#'               dplyr::bind_rows(dplyr::filter(strat_df, frame != 1), .) %>%
#'               dplyr::mutate(
#'                 across(c(sampled_cluster, sampled), ~ ifelse(is.na(.), 0, .)),
#'                 sampled = dplyr::if_else(sampled_cluster == 0, 0, sampled),
#'                 weight = weight * cluster_prop)
#'           } else {
#'             strat_df %>%
#'               dplyr::mutate(cluster_id = 1:n(),
#'                             cluster_prop = 1/n()) %>%
#'               dplyr::filter(frame == 1) %>%
#'               dplyr::mutate(sampled =
#'                               sample(c(rep(1,min(n(), target_n_pps)),
#'                                        rep(0,max(0, n() - target_n_pps)))),
#'                             sampled_cluster = sampled,
#'                             weight = sum(frame)/sum(sampled)) %>%
#'               dplyr::bind_rows(dplyr::filter(strat_df, frame != 1), .) %>%
#'               dplyr::mutate(across(c(sampled_cluster, sampled), ~ ifelse(is.na(.), 0, .)))
#'           }
#'         }
#'       )
#'
#'     data <-
#'       suppressMessages(
#'         temp_data %>%
#'           dplyr::rename_with(.cols = c(frame, cluster_prop, strata_prop, sampled_cluster, weight),
#'                              ~ paste0(sample_label, "_", gsub("^sampled$", "", x = .))) %>%
#'           dplyr::mutate("{ sample_label }" := sampled,
#'                         "{ sample_label }_cluster" := cluster_id,
#'                         "{ sample_label }_strata" := strata_id) %>%
#'           dplyr::select(temp_id, starts_with(sample_label)) %>%
#'           dplyr::left_join(data, ., by = "temp_id") %>%
#'           dplyr::select(-temp_id)
#'       )
#'
#'     if (drop_nonsampled) {
#'       data %<>% dplyr::filter(dplyr::if_all(dplyr::all_of(sample_label), ~ . == 1))
#'     }
#'
#'     return(data)
#'
#'   }
#'
#'
#'
#' #' Draw time-location (TLS) sample from single study
#' #'
#' #' Sampling handler for drawing TLS sample with given characteristics from individual study population
#' #'
#' #' @param data pass-through population data frame
#' #' @param sample_label character string that is used as prefix for all variables generated by TLS sampling. Default is 'tls'
#' #' @param drop_nonsampled logical indicating whether to drop units that are not sampled. Default is \code{FALSE}
#' #' @param hidden_var character string specifying hidden variable name (associated probability of visibility should be named \code{p_visible_[hidden_var]}). Defaults to "hidden" for the simulations
#' #' @param target_n_clusters target number of clusters (time-locations). Clusters are always sampled proportionally to their size in terms of number of hidden population members
#' #' @param target_cluster_type character string specifying the type of TLS sampling within each location. Either "prop" in which case \code{target_per_cluster} should give a share or "fixed" in which case \code{target_per_cluster} should be integer. Default is proportional sampling within clusters
#' #' @param target_per_cluster numeric target for within cluster. Either share for proportional sampling or integer for fixed sampling. If in any cluster fixed number of units required for sampling is larger than the number of units in cluster, the whole cluster is sampled and the warning is produced
#' #' @param clusters character string containing names of all locality names in the study population data frame
#' #'
#' #' @return Population or sample data frame for single study with TLS sample characteristics added
#' #'  \describe{
#' #'   \item{[sample_label]}{Sampling indicator}
#' #'   \item{[sample_label]_cluster}{Time-locations at which subject was encountered first}
#' #'   \item{[sample_label]_loc_present}{Sampled time-locations at which subject is present}
#' #'   \item{[sample_label]_loc_sampled}{Sampled time-locations at which subject was sampled}
#' #'   \item{[sample_label]_locs_sampled}{Sampled time-locations}
#' #'   \item{[sample_label]_weight}{Sampling weight without visibility}
#' #'   \item{[sample_label]_weight_visible}{Sampling weight with visibility}
#' #'  }
#' #'
#' #' @keywords internal
#' #'
#' #'
#' sample_tls_tidy <-
#'   function(data,
#'            sample_label = "tls", drop_nonsampled = FALSE,
#'            hidden_var = "hidden",
#'            target_n_clusters,
#'            target_cluster_type = c("prop", "fixed"),
#'            target_per_cluster,
#'            clusters
#'   ) {
#'
#'     warning("This function was deprecated due to lacking performance. It has been replaced with an optimized version. Consider using sample_tls().")
#'
#'
#'     target_cluster_type <- match.arg(target_cluster_type)
#'
#'     if (length(clusters) < target_n_clusters)
#'       stop("Number of requested locations for TLS sample exceeds number of available locations")
#'
#'
#'     if ((target_cluster_type == "prop" & (target_per_cluster > 1 | target_per_cluster < 0)) |
#'         (target_cluster_type == "fixed" & ((target_per_cluster %% 1) != 0))) {
#'       stop("There is a mismatch between type and target specified for within clusters sampling")
#'     }
#'
#'     sampling_probs <-
#'       data %>%
#'       # think about adding visibility here as well
#'       purrr::when(
#'         !is.null(hidden_var) ~ dplyr::filter(., across(all_of(hidden_var), ~ .x == 1)),
#'         ~ .) %>%
#'       dplyr::summarise(dplyr::across(dplyr::all_of(clusters), sum, .names = "{.col}")) %>%
#'       { tibble(loc = colnames(.), n = unlist(.), loc_sampling_prob = unlist(.)/sum(.)) }
#'
#'     sampled_locs <-
#'       sampling_probs %>%
#'       {
#'         sample(x = .$loc, size = target_n_clusters, prob = .$loc_sampling_prob)
#'       }
#'
#'     if (target_cluster_type == "prop") {
#'       sampling_probs %<>%
#'         dplyr::filter(loc %in% sampled_locs) %>%
#'         dplyr::mutate(
#'           target_n = n * target_per_cluster,
#'           target_n =
#'             dplyr::if_else((target_n %% 1) <= median(target_n %% 1),
#'                            floor(target_n),
#'                            ceiling(target_n)),
#'           unit_sampling_prob = target_n/n
#'         )
#'     } else if (target_cluster_type == "fixed") {
#'       sampling_probs %<>%
#'         dplyr::filter(loc %in% sampled_locs) %>%
#'         dplyr::mutate(
#'           target_n = min(target_per_cluster, n),
#'           unit_sampling_prob = target_n/n
#'         )
#'
#'       if (any(target_per_cluster > sampling_probs$n))
#'         warning("For some clusters number of units required by sampling procedure was larger than cluster size. Probability of unit sampling within those clusters is set to 1")
#'     }
#'
#'     data %<>% dplyr::mutate(temp_id = 1:n())
#'
#'     data %<>%
#'       dplyr::filter(if_any(.cols = all_of(sampled_locs), ~ . == 1)) %>%
#'       purrr::when(
#'         !is.null(hidden_var) ~ dplyr::filter(., if_all(all_of(hidden_var), ~ .x == 1)),
#'         ~ .) %>%
#'       {
#'         dplyr::bind_rows(
#'           lapply(sampled_locs,
#'                  function(x) tibble(.[.[[x]] == 1,], loc = x))
#'         )
#'       } %>%
#'       dplyr::left_join(., sampling_probs, by = "loc") %>%
#'       dplyr::group_by(loc) %>%
#'       dplyr::mutate(
#'         sampled = sample(c(rep(1, unique(target_n)),
#'                            rep(0, n() - unique(target_n))),
#'                          prob = if (!is.null(hidden_var)) get(paste0("p_visible_", hidden_var)))
#'       )%>%
#'       dplyr::group_by(temp_id, .add = FALSE) %>%
#'       purrr::when(
#'         !is.null(hidden_var) ~ dplyr::summarise(
#'           .,
#'           sampled = as.integer(any(sampled == 1)),
#'           loc_present = paste0(loc, collapse = ";"),
#'           loc_sampled = dplyr::if_else(any(sampled == 1),
#'                                        paste0(loc[which(sampled == 1)], collapse = ";"),
#'                                        NA_character_),
#'           weight = 1/(1 - prod(c(1 - unit_sampling_prob))),
#'           weight_visible = 1/(1 - prod(1 - get(paste0("p_visible_", hidden_var)) * unit_sampling_prob))
#'         ),
#'         ~ dplyr::summarise(
#'           .,
#'           sampled = as.integer(any(sampled == 1)),
#'           loc_present = paste0(loc, collapse = ";"),
#'           loc_sampled = dplyr::if_else(any(sampled == 1),
#'                                        paste0(loc[which(sampled == 1)], collapse = ";"),
#'                                        NA_character_),
#'           weight = 1/(1 - prod(c(1 - unit_sampling_prob))),
#'           weight_visible = 1/(1 - prod(1 - unit_sampling_prob))
#'         )) %>%
#'       dplyr::ungroup() %>%
#'       dplyr::filter(sampled == 1) %>%
#'       dplyr::select(-sampled) %>%
#'       dplyr::mutate(
#'         locs_sampled = paste0(na.omit(unique(loc_sampled)), collapse = ";")
#'       ) %>%
#'       dplyr::rename_with(~ paste0(sample_label, "_", .), .cols = -temp_id) %>%
#'       {
#'         left_join(
#'           dplyr::mutate(data, "{ sample_label }" := as.integer(temp_id %in% .$temp_id)),
#'           .,
#'           by = "temp_id"
#'         )
#'       } %>%
#'       dplyr::select(-temp_id)
#'
#'     if (drop_nonsampled) data %<>% dplyr::filter(if_all(all_of(sample_label), ~ . == 1))
#'
#'     return(data)
#'
#'   }
#'
#'
#' #' Draw respondent-driven sample (RDS) sample from single study using tidyverse
#' #'
#' #' Sampling handler for drawing RDS sample with given characteristics from individual study population
#' #'
#' #' @param data pass-through population data frame
#' #' @param sample_label character string that is used as prefix for all variables generated by RDS sampling (sample identifier, recruiter ID, wave, time of show-up)
#' #' @param hidden_var character string specifying hidden variable name (associated probability of visibility should be named \code{p_visible_[hidden_var]}). Defaults to "hidden" for the simulations
#' #' @param n_seed number of seeds randomly drawn from members of hidden population (group K)
#' #' @param n_coupons number of unique coupons given to each study participant
#' #' @param n_waves number of waves allowed. Disregarded in \code{target_type = 'sample'}
#' #' @param target_type one of 'sample' or 'waves'
#' #' @param target_n_rds numeric target size of RDS sample. If \code{target_type = "sample"}, this gives maximum number of respondents to be sampled (right now the RDS network can also end before reaching sample size target). If \code{target_type = "waves"}, this gives maximum number of waves of recruitment allowed
#' #' @param add_seeds numeric indicating how many seeds to add at a time if target sample size is not reached with initial seeds. Additional seeds are randomly drawn from non-sampled hidden population members. Defaults to \code{NULL} that does not allow adding seeds
#' #' @param arrival_rate numeric rate of respondent arrival per interval of time (e.g. per hour or day). Defaults to .5
#' #' @param drop_nonsampled logical indicating whether to drop units that are not sampled. Default is \code{FALSE}
#' #'
#' #' @return Population or sample data frame for single study with RDS sample characteristics added
#' #'  \describe{
#' #'   \item{[sample_label]}{Sampling indicator}
#' #'   \item{[sample_label]_from}{ID of respondent who enrolled current respondent}
#' #'   \item{[sample_label]_t}{Time at which respondent was enrolled}
#' #'   \item{[sample_label]_wave}{Number of steps respondent is away from the seed (seeds are wave 0)}
#' #'   \item{[sample_label]_own_coupon}{ID of the coupon with which respondent was enrolled}
#' #'   \item{[sample_label]_coupon_[1-9]}{IDs of the coupons that were given to respondent for enrollment from their network}
#' #'  }
#' #'
#' #' @keywords internal
#' #'
#' #' @examples
#' #' \dontrun{
#' #'    sample_rds_tidy(
#' #'      data = data,
#' #'      sample_label = "rds",
#' #'      hidden_var = "hidden" ,
#' #'      n_seed = 20,
#' #'      n_coupons = 3,
#' #'      add_seeds = 3,
#' #'      target_type = "sample",
#' #'      target_n_rds = 60,
#' #'      arrival_rate = .5
#' #'    )
#' #' }
#' #'
#' #'
#' sample_rds_tidy <-
#'   function(data,
#'            sample_label = "rds",
#'            hidden_var,
#'            target_type = c("sample", "waves"),
#'            n_seed, n_coupons, n_waves = NULL,
#'            target_n_rds,
#'            add_seeds = NULL,
#'            arrival_rate = .5,
#'            drop_nonsampled = FALSE) {
#'
#'     warning("This function was deprecated due to lacking performance. It has been replaced with an optimized version. Consider using sample_rds().")
#'     target_type <- match.arg(target_type)
#'
#'     data %<>%
#'       dplyr::as_tibble() %>%
#'       dplyr::mutate(links_list = retrieve_graph(links) %>% igraph::as_adj_list())
#'
#'     # generate arrival times using rate of 2 per unit of time
#'     # consider making rate a lambda parameter later
#'     .arrival_time <-
#'       rexp(n = length(data$name[dplyr::pull(data, hidden_var) == 1]),
#'                   rate = arrival_rate) %>%
#'       base::cumsum()
#'
#'     if (length(data$name[dplyr::pull(data, hidden_var) == 1]) <= n_seed) {
#'       .seeds <- data$name[dplyr::pull(data, hidden_var) == 1]
#'     } else {
#'       .seeds <-
#'         sample(
#'           # sample out of all people in hidden population
#'           x = data$name[dplyr::pull(data, hidden_var) == 1],
#'           # select only prescribed number of subjects
#'           size = n_seed,
#'           prob =
#'             dplyr::pull(data, paste0("p_visible_", hidden_var))[dplyr::pull(data, hidden_var) == 1],
#'           replace = FALSE)
#'     }
#'
#'     # at t=1 only seeds are sampled
#'     .sampled <-
#'       dplyr::tibble(
#'         name = .seeds,
#'         from = -999,
#'         t = .arrival_time[1:length(.seeds)],
#'         wave = 1,
#'         !!hidden_var := dplyr::pull(data, hidden_var)[data$name %in% .seeds],
#'         own_coupon = as.character(1:length(.seeds)))
#'
#'     for (j in 1:n_coupons) {
#'       .sampled[paste0("coupon_", j)] <- paste0(.sampled$own_coupon, "-", j)
#'     }
#'
#'     # n_coupons of their links (not just in hidden pop) are eligible
#'     .eligible <-
#'       .seeds %>%
#'       # sample from each seed links
#'       lapply(.,
#'              function(x){
#'
#'                # presume that only hidden population links can be sampled
#'                .available_links <-
#'                  data %>%
#'                  dplyr::filter(
#'                    name %in% data$links_list[[which(data$name == x)]],
#'                    eval(expr = parse(text = paste0(hidden_var, " == 1")))
#'                  ) %>%
#'                  .$name
#'
#'
#'                if (length(.available_links) > 0) {
#'                  if (length(.available_links) == 1) {
#'                    dplyr::tibble(
#'                      from = x,
#'                      to = as.integer(.available_links),
#'                      own_coupon = sample(
#'                        x = unname(unlist(.sampled[.sampled$name == x,
#'                                                   grep("^coupon\\_", names(.sampled))])),
#'                        size = 1))
#'                  } else {
#'                    dplyr::tibble(
#'                      from = x,
#'                      to = sample(x = as.integer(.available_links),
#'                                  size = min(length(.available_links), n_coupons),
#'                                  replace = FALSE),
#'                      own_coupon = sample(
#'                        x =
#'                          unname(
#'                            unlist(.sampled[.sampled$name == x,
#'                                            grep("^coupon\\_", names(.sampled))])),
#'                        size = min(length(.available_links), n_coupons)))
#'                  }
#'                } else {
#'                  tibble(from = integer(), to = integer(), own_coupon = character())
#'                }
#'              }) %>%
#'       dplyr::bind_rows() %>%
#'       # join in eligible showup rates
#'       dplyr::left_join(., data[, c("name", paste0(c("p_visible_", ""), hidden_var))],
#'                        by = c(to = "name")) %>%
#'       # drop those who won't show up and those who were already sampled
#'       # also record wave - number of links from seed
#'       dplyr::mutate(
#'         showup = rbinom(n = dplyr::n(), size = 1,
#'                         prob = get(paste0("p_visible_", hidden_var))),
#'         wave = 2) %>%
#'       dplyr::filter(showup == 1,
#'                     !(to %in% .sampled$name),
#'                     !duplicated(to)) %>%
#'       # do we exclude non-members of hidden population?
#'       dplyr::select(from, to, wave, dplyr::all_of(hidden_var), own_coupon)
#'
#'     # Next, run the loop with the same procedure
#'     .t <- (length(.seeds) + 1)
#'
#'     repeat {
#'
#'       # if ran out of links add seeds at random from those not sampled
#'       # this also allows for adding links if initial seeds have no connections
#'       if (target_type == "sample" & (nrow(.eligible) == 0)) {
#'         if (is.numeric(add_seeds)) {
#'
#'           # get nodes that were not sampled yet
#'           .nonsampled <-
#'             setdiff(data$name[dplyr::pull(data, hidden_var) == 1], .sampled$name)
#'
#'           if (length(.nonsampled) == 0) {
#'             break
#'           } else {
#'
#'             # check whether we have enough remaining seeds to sample according to add_seeds
#'             .new_seeds <- min(length(.nonsampled), add_seeds)
#'
#'             # sample new seeds according to number of additional seeds specified
#'             .new_seeds <-
#'               sample(
#'                 x = .nonsampled,
#'                 size = .new_seeds,
#'                 prob =
#'                   dplyr::pull(data, paste0("p_visible_", hidden_var))[data$name %in% .nonsampled],
#'                 replace = FALSE)
#'
#'             .eligible <-
#'               dplyr::tibble(
#'                 from = -999, to = .new_seeds, wave = 1,
#'                 !!hidden_var := dplyr::pull(data, hidden_var)[data$name %in% .new_seeds],
#'                 own_coupon = as.character((n_seed+1):(n_seed + length(.new_seeds))))
#'
#'             # update number of seeds
#'             n_seed <- n_seed + length(.new_seeds)
#'           }
#'         } else {
#'           break
#'         }
#'       }
#'
#'       # if ran out of links and waves are target - break
#'       if (target_type == "waves" & (nrow(.eligible) == 0))
#'         break
#'
#'       # sample 1 individual weighting by wave
#'       # (so that earlier waves are more likely to be sampled)
#'       .new <-
#'         dplyr::slice_sample(.eligible, n = 1, replace = FALSE, weight_by = 1/wave)
#'
#'       # move new sampled from eligible to sampled
#'       .sampled <-
#'         dplyr::tibble(name = .new$to,
#'                       from = .new$from,
#'                       t = .arrival_time[.t],
#'                       wave = .new$wave,
#'                       !!hidden_var := dplyr::pull(.new, hidden_var),
#'                       own_coupon = unname(.new$own_coupon)) %>%
#'         dplyr::bind_cols(.,
#'                          `names<-`(as.list(paste0(.new$own_coupon, "-", 1:n_coupons)),
#'                                    paste0("coupon_", 1:n_coupons))) %>%
#'         dplyr::bind_rows(.sampled, .)
#'
#'       .eligible %<>% dplyr::filter(to != .new$to)
#'
#'       # presume that only hidden population links can be sampled
#'       if ((target_type == "waves") & ifelse(is.null(n_waves), FALSE, (.new$wave == n_waves))) {
#'         .new_available_links <- c()
#'       } else {
#'         .new_available_links <-
#'           data %>%
#'           dplyr::filter(
#'             name %in% data$links_list[[which(data$name == .new$to)]],
#'             eval(expr = parse(text = paste0(hidden_var, " == 1")))
#'           ) %>%
#'           .$name
#'       }
#'
#'       # add new eligible links using the same procedure as above
#'       if (length(.new_available_links) > 0) {
#'
#'         .eligible <-
#'           .new_available_links %>%
#'           purrr::when(
#'             length(.new_available_links) == 1 ~
#'               dplyr::tibble(from = .new$to,
#'                             to = .,
#'                             own_coupon = sample(
#'                               x = unname(unlist(.sampled[.sampled$name == .new$to,
#'                                                          grep("^coupon\\_", x = names(.sampled))])),
#'                               size = 1)),
#'             dplyr::tibble(from = .new$to,
#'                           to = sample(x = as.integer(.), size = min(length(.), n_coupons),
#'                                       replace = FALSE),
#'                           own_coupon = sample(
#'                             x =
#'                               unname(
#'                                 unlist(.sampled[.sampled$name == .new$to,
#'                                                 grep("^coupon\\_", x = names(.sampled))])),
#'                             size = min(length(.), n_coupons)))) %>%
#'           # join in eligible showup rates
#'           dplyr::left_join(., data[, c("name", paste0(c("p_visible_", ""), hidden_var))],
#'                            by = c(to = "name")) %>%
#'           # drop those who won't show up and those who were already samples
#'           dplyr::mutate(showup = rbinom(n = dplyr::n(), size = 1,
#'                                         prob = get(paste0("p_visible_", hidden_var))),
#'                         wave = .new$wave + 1) %>%
#'           dplyr::filter(showup == 1,
#'                         !(to %in% .sampled$name),
#'                         !duplicated(to)) %>%
#'           # do we exclude non-members of hidden population?
#'           # dplyr::filter(hidden == 1) %>%
#'           dplyr::select(from, to, wave, dplyr::all_of(hidden_var), own_coupon) %>%
#'           dplyr::bind_rows(.eligible, .)
#'       }
#'
#'       # break if reached desired sample size (for waves or sample target)
#'       if (nrow(.sampled) >= target_n_rds)
#'         break
#'
#'       .t <- .t + 1
#'
#'     }
#'
#'     data[,sample_label] <- as.integer(data$name %in% .sampled$name)
#'     names(.sampled) %<>% {c(.[1], paste0(sample_label, "_", .[-1]))}
#'
#'     data %<>%
#'       dplyr::left_join(., .sampled, by = "name") %>%
#'       dplyr::select(-links_list)
#'
#'     if (drop_nonsampled) data %<>%
#'       dplyr::filter(dplyr::if_all(dplyr::all_of(sample_label), ~ . == 1))
#'
#'     return(data)
#'
#'   }
#'
#'
#'
#'
#' #' Get individual study estimands
#' #'
#' #' @param data pass-through population data frame
#' #' @param known_pattern character string containing regular expression to match known group names in the study population dataset
#' #' @param hidden_var character string containing hidden group name in the study population data frame
#' #'
#' #' @keywords internal
#' #'
#' #' @return Estimands data frame for single study
#' #'
#' get_study_estimands_tidy <- function(data,
#'                                      known_pattern = "^known(\\_\\d|\\d)?$",
#'                                      hidden_var = "hidden") {
#'
#'
#'   known_groups <- names(data)[grepl(x = names(data), pattern = known_pattern)]
#'
#'   .out <-
#'     data %>%
#'     # dplyr::mutate(
#'     #   degree_all = purrr::map_int(links, ~ length(strsplit(.x, split = ";")[[1]])),
#'     #   degree_hidden =
#'     #     purrr::map_int(links,
#'     #                    ~ sum(unlist(data[data$name %in% strsplit(.x, split = ";")[[1]],
#'     #                                      grep(pattern = paste0("^", hidden_var, "$"),
#'     #                                           names(data))])))
#'     # ) %>%
#'     dplyr::summarise(
#'       dplyr::across(c(dplyr::matches(known_pattern),
#'                       dplyr::matches(paste0("^", hidden_var, "$"))),
#'                     list(size = sum), na.rm  = TRUE),
#'       dplyr::across(c(dplyr::matches(known_pattern),
#'                       dplyr::matches(paste0("^", hidden_var, "$"))),
#'                     list(prev = mean), na.rm  = TRUE),
#'       # dplyr::across(c(degree_all, degree_hidden), mean, na.rm  = TRUE)
#'     ) %>%
#'     {
#'       data.frame(
#'         inquiry = names(.),
#'         estimand = unname(t(.)),
#'         stringsAsFactors = FALSE)
#'     }# %>%
#'   # {.[!grepl(pattern = known_pattern, .$inquiry),]}
#'
#'   if (length(known_groups) != 0) {
#'     for (i in seq_along(known_groups)) {
#'       .out <-
#'         data %>%
#'         dplyr::filter(eval(parse(text = paste0(known_groups[i], " == 1")))) %>%
#'         dplyr::summarise(
#'           dplyr::across(c(dplyr::matches(paste0("^", hidden_var, "$"))),
#'                         list(size = sum), na.rm  = TRUE),
#'           dplyr::across(c(dplyr::matches(paste0("^", hidden_var, "$"))),
#'                         list(prev = mean), na.rm  = TRUE),
#'           # dplyr::across(c(degree_all, degree_hidden), mean, na.rm  = TRUE)
#'         ) %>%
#'         {
#'           data.frame(
#'             inquiry = paste0(names(.), "_in_", known_groups[i]),
#'             estimand = unname(t(.)),
#'             stringsAsFactors = FALSE)
#'         } %>%
#'         bind_rows(.out, .)
#'     }
#'   }
#'
#'   return(.out)
#' }
#'
#'
#'
#' #' NSUM estimatior
#' #'
#' #' @param data pass-through population data frame
#' #' @param known character vector containing names of known groups
#' #' @param hidden character vector containing names of hidden groups
#' #' @param survey_design a formula describing the design of the survey
#' #' @param degree_ratio numeric value between 0 and 1 representing degree ratio
#' #' @param transmission_rate numeric value between 0 and 1 representing information transmission rate
#' #' @param n_boot number of bootstrap resamples
#' #' @param parallel_boot logical, whether to compute bootstrap samples in parallel using \code{foreach} package
#' #' @param prefix character prefix used for PPS sample variables
#' #' @param label character string describing the estimator
#' #'
#' #' @return Data frame of NSUM estimates for a single study with PPS sample
#' #'
#' #' @keywords internal
#' #'
#' #' @reference {Dennis M. Feehan, Matthew J. Salganik. “The networkreporting package.” (2014). \url{https://cran.r-project.org/package=networkreporting}.
#' #' @reference {Dennis M. Feehan, Matthew J. Salganik. “The surveybootstrap package.” (2016). \url{https://cran.r-project.org/package=surveybootstrap.}}
#' #' @reference {Salganik, Matthew J. "Variance estimation, design effects, and sample size calculations for respondent-driven sampling." Journal of Urban Health 83, no. 1 (2006): 98. \url{https://doi.org/10.1007/s11524-006-9106-x}}
#' #'
#' get_study_est_nsum_old <- function(data,
#'                                known,
#'                                hidden,
#'                                survey_design = ~ pps_cluster + strata(pps_strata),
#'                                degree_ratio = 1,
#'                                transmission_rate = 1,
#'                                n_boot = 1000,
#'                                parallel_boot = FALSE,
#'                                prefix = "pps",
#'                                label = "nsum") {
#'
#'   if (parallel_boot) {
#'     requireNamespace(c("doParallel", "parallel"))
#'     doParallel::registerDoParallel(cores = parallel::detectCores() - 1)
#'   }
#'
#'   .data_mod <- data[get(prefix) == 1,]
#'
#'   .known_pops <-
#'     .data_mod %>%
#'     dplyr::select(total, all_of(paste0("total_", known))) %>%
#'     dplyr::rename_with(
#'       .cols = starts_with("total_"), ~ paste0(gsub("^total\\_", "", .), "_visible_out")) %>%
#'     apply(., MARGIN = 2, unique)
#'
#'   .data_mod$d_est <-
#'     purrr::quietly(networkreporting::kp.individual.estimator)(
#'       resp.data = .data_mod,
#'       alter.popn.size = .known_pops[1],
#'       known.populations = names(.known_pops[2:length(.known_pops)]),
#'       total.kp.size = sum(.known_pops[2:length(.known_pops)]))$result$dbar.Fcell.F
#'
#'   .fit_nsum <-
#'     purrr::quietly(networkreporting::nsum.estimator)(
#'       survey.data = .data_mod,
#'       d.hat.vals = "d_est",
#'       total.popn.size = .known_pops[1],
#'       killworth.se = FALSE,
#'       y.vals = hidden,
#'       weights = paste0(prefix, "_weight"),
#'       missing = "ignore",
#'       deg.ratio = degree_ratio,
#'       tx.rate = transmission_rate)$result
#'
#'   .fit_nsum_boot <-
#'     get_rescaled_boot(data = .data_mod,
#'                       survey_design = survey_design,
#'                       n_boot = n_boot) %>%
#'     plyr::llply(
#'       .data = .,
#'       .fun = function(wgt) {
#'         .data_mod %>%
#'           dplyr::mutate(index = 1:dplyr::n()) %>%
#'           dplyr::left_join(., wgt, by = "index") %>%
#'           dplyr::mutate(rescaled_weight = weight.scale * get(paste0(prefix, "_weight"))) %>%
#'           networkreporting::nsum.estimator(
#'             survey.data = .,
#'             d.hat.vals = "d_est",
#'             total.popn.size = .known_pops[1],
#'             killworth.se = FALSE,
#'             y.vals = hidden,
#'             weights = "rescaled_weight",
#'             missing = "ignore",
#'             deg.ratio = degree_ratio,
#'             tx.rate = transmission_rate)
#'       },
#'       .parallel = parallel_boot) %>%
#'     dplyr::bind_rows()
#'
#'   return(
#'     # data.frame(
#'     #   estimator = paste0(c("hidden_size_", "degree_"), label),
#'     #   estimate = c(unname(.fit_nsum$estimate),
#'     #                unname(.fit_nsum$sum.d.hat/.fit_nsum$estimate)),
#'     #   se = c(sd(.fit_nsum_boot$estimate),
#'     #          sd(.fit_nsum_boot$sum.d.hat/.fit_nsum_boot$estimate)),
#'     #   inquiry = c("hidden_size", "degree_all"))
#'     data.frame(
#'       estimator = paste0(c("hidden_size_"), label),
#'       estimate = c(unname(.fit_nsum$estimate)),
#'       se = c(sd(.fit_nsum_boot$estimate)),
#'       inquiry = c("hidden_size"))
#'   )
#' }
#'
#'
#' #' Get Parameters and Design Features of Specific Study
#' #'
#' #' @param study Character string giving name of the study to pull. Have to match the names provided in Google Spreadsheet
#' #' @param ss Character string giving ID of the Google spreadsheet
#' #' @param sheet Character string giving name of the sheet in the Google spreadsheet to read
#' #'
#' #' @keywords internal
#' #'
#' read_single_study_params <- function(
#'     study,
#'     ss = "1HwMM6JwoGALLMTpC8pQRVzaQdD2X71jpZNhfpKTnF8Y",
#'     sheet = "study_params_readable"
#' ) {
#'
#'   if (!requireNamespace("googlesheets4", quietly = TRUE)) {
#'     stop(
#'       "Package \"googlesheets4\" must be installed to use this function.",
#'       call. = FALSE
#'     )
#'   }
#'
#'   labels <-
#'     suppressMessages(
#'       googlesheets4::read_sheet(
#'         ss = ss,
#'         sheet = sheet,
#'         skip = 0,
#'         n_max = 4,
#'         col_names = FALSE,
#'         .name_repair = "minimal")
#'     ) %>%
#'     unname %>%
#'     t %>%
#'     `colnames<-`(c("family", "type", "name", "description")) %>%
#'     dplyr::as_tibble() %>%
#'     dplyr::mutate(
#'       # prior_type = purrr::map_chr(family, ~ strsplit(.x, "_")[[1]][2]),
#'       # family = purrr::map_chr(family, ~ strsplit(.x, "_")[[1]][1]),
#'       colnames = dplyr::if_else(is.na(family), name, paste0(family, "_", name))
#'     )
#'
#'
#'   study_data <-
#'     suppressMessages(
#'       googlesheets4::read_sheet(
#'         ss = ss,
#'         sheet = sheet,
#'         col_names = FALSE,
#'         na = "NA",
#'         col_types = paste0(labels$type, collapse = ""),
#'         .name_repair = "minimal",
#'         skip = 4) %>%
#'         `names<-`(dplyr::pull(labels, colnames))
#'     ) %>%
#'     dplyr::filter(study_id == study) %>%
#'     {dplyr::tibble(colnames = colnames(.), params := t(.)[,1])}
#'
#'   .out <- dplyr::left_join(labels, study_data, by = "colnames")
#'
#'   for (i in c("rds", "pps", "tls")) {
#'     if (.out$params[which(.out$colnames == i)] == 0) {
#'       .out <- dplyr::filter(.out, is.na(family) | (family != i))
#'     }
#'   }
#'
#'   return(
#'     dplyr::select(.out,
#'                   "Name" = name,
#'                   "Description" = description,
#'                   "Value" = params)
#'   )
#'
#' }
#'
#'
#' #' Get Parameters and Design Features of Specific Study
#' #'
#' #' @param study Character string giving name of the study to pull. Have to match the names provided in Google Spreadsheet
#' #' @param ss Character string giving ID of the Google spreadsheet
#' #' @param sheet Character string giving name of the sheet in the Google spreadsheet to read
#' #'
#' #' @keywords internal
#' #'
#' #' @importFrom magrittr `%>%` `%$%`
#' #' @importFrom DeclareDesign declare_population declare_sampling declare_inquiry declare_estimator
#' get_single_study_design <- function(
#'     study = "johnjay_tanzania",
#'     ss = "1HwMM6JwoGALLMTpC8pQRVzaQdD2X71jpZNhfpKTnF8Y",
#'     sheet = "study_params_readable") {
#'
#'   .type_map <-
#'     list(population = DeclareDesign::declare_population,
#'          sample     = DeclareDesign::declare_sampling,
#'          inquiry    = DeclareDesign::declare_inquiry,
#'          est        = DeclareDesign::declare_estimator)
#'
#'   .design_list <- read_study_params(ss, sheet)[[study]]
#'
#'   .declare_list <-
#'     list(study_population   = .design_list$pop,
#'          study_sample_rds   = .design_list$samples$rds,
#'          study_sample_pps   = .design_list$samples$pps,
#'          study_sample_tls   = .design_list$samples$tls,
#'          study_inquiry      = .design_list$inquiries,
#'          est_sspse          = .design_list$estimators$rds$sspse,
#'          est_chords         = .design_list$estimators$rds$chords,
#'          est_multi          = .design_list$estimators$rds$multiplier,
#'          est_ht_pps         = .design_list$estimators$pps$ht,
#'          est_nsum_pps       = .design_list$estimators$pps$nsum,
#'          est_recap_rds_pps  = .design_list$estimators$rds_pps$recap1,
#'          est_ht_tls         = .design_list$estimators$tls$ht,
#'          est_nsum_tls       = .design_list$estimators$tls$nsum,
#'          est_recap_tls      = .design_list$estimators$tls$recap
#'     )
#'
#'   for (i in seq_along(.declare_list)) {
#'     if (!is.null(.declare_list[[i]])) {
#'       .declare_fun <- .type_map[[stringr::str_extract(pattern = paste0(names(.type_map), collapse = "|"),
#'                                                       string = names(.declare_list)[i])]]
#'       assign(names(.declare_list)[i], eval(as.call(c(list(.declare_fun), .declare_list[[i]]))))
#'     }
#'   }
#'
#'   return(
#'     eval(parse(text = paste0(names(.declare_list)[sapply(.declare_list, function(x) !is.null(x))], collapse = " + ")))
#'   )
#'
#' }
#' #' Handler for Simulation of Network Using Block Model
#' #'
#' #' Draw a Simulation of Network Using Block Matrix of Probabilities of Edges Based on Group Memberships
#' #'
#' #' @param N number of units in population
#' #' @param K number of groups
#' #' @param prev_K named numeric vector of prevalence for each group with last group being hidden
#' #' @param rho_K numeric vector of correlations in group memberships
#' #' @param p_edge_within named list of numeric vectors giving probability of link between in-group members and out-group members for each of groups. The order of objects in list have to follow the order of \code{prev_K}
#' #' @param p_edge_between named list of numeric values giving probability of link between in- and out-group member for each of groups. The order of objects in list have to follow the order of \code{prev_K}
#' #' @param directed logical, whether links are directed or undirected
#' #'
#' #' @return igraph network object with vertex attribute \code{type} in the binary coded format (consists of 0's and 1's only)
#' #'
#' #' @examples
#' #' \dontrun{
#' #'  sim_block_network(N = 2000, K = 3,
#' #'                    prev_K = c(frame = .7,
#' #'                               known = .2,
#' #'                               hidden = .05),
#' #'                    rho_K = c(.05, .01, .01),
#' #'                    p_edge_within =
#' #'                      list(frame = c(.3,.3),
#' #'                           known = c(.3,.3),
#' #'                           hidden = c(.3,.9)),
#' #'                    p_edge_between =
#' #'                      list(frame = 0.3,
#' #'                           known = 0.3,
#' #'                           hidden = 0.05),
#' #'                    directed = FALSE)
#' #' }
#' #'
#' #' @importFrom igraph sample_pref vertex_attr
#' #' @importFrom magrittr `%>%` `%$%`
#' #' @importFrom plyr mapvalues
#' #'
#' sim_block_network_old <-
#'   function(N, K,
#'            prev_K,
#'            rho_K,
#'            p_edge_within,
#'            p_edge_between,
#'            directed = FALSE) {
#'
#'     if (is.null(names(prev_K)) | length(unique(names(prev_K))) != K)
#'       stop("prev_K have to be named vector of prevalences with unique names identifying each group. The last element have to correspond to hidden group prevalence")
#'
#'     type_names <- gen_group_types(K = K)
#'
#'     g <-
#'       igraph::sample_pref(nodes = N, types = 2^K,
#'                           type.dist = gen_group_sizes(N = N,
#'                                                       prev_K = prev_K, rho_K = rho_K,
#'                                                       .ord = type_names),
#'                           pref.matrix = gen_block_matrix(p_edge_within = p_edge_within,
#'                                                          p_edge_between = p_edge_between,
#'                                                          .ord = type_names),
#'                           fixed.sizes = TRUE, directed = directed, loops = FALSE)
#'
#'     igraph::vertex_attr(g)$type <-
#'       igraph::vertex_attr(g)$type %>%
#'       {plyr::mapvalues(x = ., from = unique(.), to = type_names)}
#'
#'     return(g)
#'
#'   }
#'
#'
#' get_multi_estimands <- function(data, study_designs) {
#'
#'   # if (is.null(data$study) | is.null(data$sample)) stop("Data does not contain study or population variables, probably because it was not generated by get_multi_samples()")
#'
#'   return(
#'     data %>%
#'       dplyr::mutate(
#'         estimands =
#'           purrr::map2(study, population,
#'                       ~ apply_args(study = .x,
#'                                    args = study_designs[[.x]][["inquiries"]],
#'                                    df = .y))
#'       ) %>%
#'       tidyr::unnest(estimands) %>%
#'       # dplyr::mutate(inquiry = paste0(study, "_", sample, "_", inquiry)) %>%
#'       dplyr::mutate(inquiry = paste0(study, "_", inquiry)) %>%
#'       dplyr::select(inquiry, estimand) %>%
#'       dplyr::distinct(.keep_all = TRUE) %>%
#'       as.data.frame(x = ., stringsAsFactors = FALSE)
#'   )
#'
#' }
#'
#'
#'
#' get_multi_estimates <- function(data, study_designs) {
#'
#'   # if (is.null(data$study) | is.null(data$sample)) stop("Data does not contain study or population variables, probably because it was not generated by get_multi_samples()")
#'
#'   .est_out <-
#'     data.frame(estimator = character(0),
#'                estimate = numeric(0),
#'                se = numeric(0),
#'                inquiry = character(0))
#'
#'   for (st in names(study_designs)) {
#'     for (samp in names(study_designs[[st]][["samples"]])) {
#'       for (est in names(study_designs[[st]][["estimators"]][[samp]])) {
#'
#'         # cat(st, ",", samp, ",", est, "\n")
#'
#'         .est_out <-
#'           apply_args(study = st,
#'                      df = data[study == st][["population"]][[1]],
#'                      args = study_designs[[st]][["estimators"]][[samp]][[est]]) %>%
#'           dplyr::mutate(inquiry = paste0(st, "_", inquiry),
#'                         estimator = as.character(estimator)) %>%
#'           dplyr::bind_rows(.est_out, .)
#'       }
#'     }
#'   }
#'
#'   return(as.data.frame(.est_out, stringsAsFactors = FALSE))
#'
#' }
#'
#'
#'
#' #' Use Stan model to estimate bias and estimands
#' #'
#' #' @param data pass-through meta population or meta sample data frame
#' #' @param sample_label name of variable storing meta analysis sampling information
#' #' @param which_estimand name of study level estimand for meta analysis
#' #' @param benchmark named list of length 2 giving benchmark sampling-estimator pair (only accepts one value across studies for now)
#' #' @param stan_handler function that takes stan_data as input and produces compilable stan model object
#' #' @param hidden_prior list of two hyperpriors, on means and standard errors of each included sampling-estimator pairs. Names of list objects should be "mean" and "se". If one number provided for a hyperprior it gets expanded to all sampling-estimator pairs
#' #' @param rel_bias_prior list of two hyperpriors, on means and standard errors of relative bias. Names of list objects should be "mean" and "se".
#' #' @param control_params list of additional parameters to pass to Stan fit function. These can include number of iterations, chains, thinning, seed and number of cores to use
#' #'
#' #' @return Data frame of meta level estimates and pertaining estimand names
#' #'
#' #'
#' #' @import tidyselect
#' #' @import rstan
#' #' @importFrom magrittr `%>%` `%$%`
#' #' @importFrom dplyr mutate filter select group_by ungroup summarize pull arrange rename_with left_join bind_rows if_all if_any
#' #' @importFrom parallel detectCores
#' #' @importFrom plyr alply
#' get_meta_estimates_old <- function(
#'     data,
#'     sample_label = "meta",
#'     which_estimand = "hidden_size",
#'     benchmark = list(sample = "pps", estimator = "ht"),
#'     stan_handler = get_meta_stan3,
#'     hidden_prior = NULL,
#'     rel_bias_prior = list(mean = 1, se = 10),
#'     control_params = list(
#'       iter = 8000, chains = 8, thin = 10,
#'       seed = 872312,
#'       cores = 1)
#' ) {
#'
#'   # allow for NULL sample label to keep all data
#'   if (is.null(sample_label)) {
#'     .stan_data <- data
#'   } else if (is.character(sample_label)) {
#'     .stan_data <-
#'       data %>%
#'       dplyr::filter(dplyr::if_any(sample_label, ~ . == 1),
#'                     inquiry %in% which_estimand)
#'   }
#'
#'   .samp_ests <- unique(.stan_data[,c("sample", "estimator")])
#'
#'
#'   .benchmark_sample_est <-
#'     which(lapply(.samp_ests$sample, paste0, collapse = "_") ==
#'             paste0(benchmark$sample, collapse = "_") &
#'             .samp_ests$estimator == benchmark$estimator)
#'   .samp_ests <-
#'     rbind(.samp_ests[.benchmark_sample_est,],
#'           .samp_ests[-.benchmark_sample_est,])
#'   .samp_est_names <-
#'     paste0(sapply(.samp_ests$sample, paste0, collapse = "_"), "_", .samp_ests$estimator)
#'   .studies <- unique(.stan_data$study)
#'
#'   .N <- length(.studies)
#'   .K <- nrow(.samp_ests)
#'
#'   if (!is.null(hidden_prior) & is.list(hidden_prior)) {
#'     .alpha_prior <- do.call(cbind, hidden_prior)
#'
#'     if (nrow(.alpha_prior) == 1) {
#'       .alpha_prior <- .alpha_prior[rep(1, times = .N),]
#'     } else if (nrow(.alpha_prior) == .N) {
#'       .alpha_prior <- .alpha_prior[.studies,]
#'     } else {
#'       stop("There is mismatch in length of hyperpriors on mean and std. error of target parameters. Length of each element should be equal either to length of number of studies included or to 1")
#'     }
#'   } else if (is.null(hidden_prior)) {
#'     .alpha_prior <-
#'       .stan_data %>%
#'       dplyr::group_by(study) %>%
#'       dplyr::summarize(mean = sum(estimate / (se * sum(1/se))), se = 2 * mean(se / (se * sum(1/se)))) %>%
#'       { .[order(match(.$study,.studies)),] } %>%
#'       dplyr::select(-study) %>%
#'       as.matrix()
#'   } else {
#'     stop("Hyperpriors on target parameters are provided in wrong format (should be NULL or list of two numeric objects)")
#'   }
#'
#'   # get ids of unique samp-est pairs in observed data
#'   .samp_est_ids <-
#'     mapply(
#'       x = .samp_ests$sample, y = .samp_ests$estimator,
#'       function(x, y) {
#'         lapply(.stan_data$sample, paste0, collapse = "_") ==
#'           paste0(x, collapse = "_") & .stan_data$estimator == y
#'       },
#'       SIMPLIFY = FALSE
#'     )
#'
#'   .stan_data <-
#'     .samp_est_ids %>%
#'     {
#'       c(lapply(X = ., FUN = sum),
#'         lapply(X = ., FUN = function(x) which(.studies %in% .stan_data$study[x])),
#'         lapply(X = ., FUN = function(x) .stan_data$estimate[x]),
#'         lapply(X = ., FUN = function(x) .stan_data$se[x]))
#'     } %>%
#'     `names<-`(
#'       .,
#'       apply(X = expand.grid(1:.K, c("n", "observed", "est", "est_se")),
#'             MARGIN = 1,
#'             FUN = function(x) paste0(x[2], as.integer(x[1])))
#'     ) %>%
#'     c(list(N = .N, K = .K,
#'            alpha_mean = unname(.alpha_prior[,"mean"]), alpha_se = unname(.alpha_prior[,"se"]),
#'            rel_bias_mean = rel_bias_prior$mean, rel_bias_se = rel_bias_prior$se),
#'       .)
#'
#'   .stan_model <-
#'     rstan::stan_model(model_code = stan_handler(.stan_data))
#'
#'   .fit <-
#'     do.call(rstan::sampling, c(list(object = .stan_model,
#'                                     data = .stan_data,
#'                                     verbose = FALSE,
#'                                     refresh = 0),
#'                                control_params)) %>%
#'     rstan::extract(.)
#'
#'   .biases <-
#'     .fit$rel_bias %>%
#'     cbind(1, .) %>%
#'     plyr::alply(
#'       .margins = 1,
#'       .fun = function(x) {
#'         c(matrix(x, nrow = .K, ncol = .K, byrow = TRUE) /
#'             matrix(x, nrow = .K, ncol = .K, byrow = FALSE))
#'       }) %>%
#'     do.call(rbind, .) %>%
#'     plyr::aaply(.margins = 2,
#'                 .fun = function(x) c(est = mean(x),
#'                                      se = sd(x))) %>%
#'     plyr::alply(.margins = 2,
#'                 function(x) {
#'                   matrix(data = x,
#'                          nrow = .K, ncol = .K,
#'                          byrow = FALSE,
#'                          dimnames = list(.samp_est_names,
#'                                          .samp_est_names))
#'                 })
#'
#'
#'   .study_ests <-
#'     .fit$alpha %>%
#'     plyr::aaply(.margins = 2,
#'                 .fun = function(x) c(est = mean(x),
#'                                      se = sd(x)))
#'
#'   return(
#'     data.frame(estimator = c(paste0("rel_bias_", .samp_est_names, "_meta"),
#'                              paste0(.studies, "_meta")),
#'                estimate = c(unname(.biases[[1]][1,]), unname(.study_ests[,1])),
#'                se =   c(unname(.biases[[2]][1,]), unname(.study_ests[,2])),
#'                inquiry = c(paste0("rel_bias_", .samp_est_names),
#'                            paste0(.studies, "_", which_estimand)),
#'                # note assumes first item is benchmark
#'                prior = c(1, rep(rel_bias_prior[[1]], .K-1), .alpha_prior[,1]) %>% unlist,
#'                prior_sd = c(0, rep(rel_bias_prior[[2]], .K-1), .alpha_prior[,2]) %>% unlist,
#'                stringsAsFactors = FALSE
#'     )
#'   )
#'
#' }
#'
#' #' @title Helper for meta-analysis
#' #'
#' #' @description Helper for meta-analysis
#' #'
#' #' @param stan_data A list of data for Stan
#' #'
#' #' @return A string of Stan code
#' #'
#' #' @export
#' #'
#' get_meta_stan_size <- function(
#'     stan_data
#' ) {
#'
#'   return(
#'     "data {
#'          int<lower=0> N;
#'          int<lower=0> K;
#'          int<lower=0> n_ests;
#'          real<lower=0> alpha_prior_mean[N];
#'          real<lower=0> alpha_prior_sd[N];
#'          real<lower=0> bias_prior_mean[K];
#'          real<lower=0> bias_prior_sd[K];
#'          int<lower=1> study[n_ests];
#'          int<lower=1> design[n_ests];
#'          vector<lower=0>[n_ests] ests;
#'          vector<lower=0>[n_ests] ses;
#'        }
#'        parameters {
#'          vector<lower=0.01>[K] bias;
#'          vector<lower=0>[N] alpha;
#'        }
#'
#'        model {
#'
#'          target += normal_lpdf(ests | bias[design].*alpha[study], ses);
#'          alpha ~ normal(alpha_prior_mean, alpha_prior_sd);
#'          bias ~ lognormal(log(bias_prior_mean), bias_prior_sd);
#'
#'        }"
#'   )
#' }
#'
#' #' @title Helper for meta-analysis
#' #'
#' #' @description Helper for meta-analysis
#' #'
#' #' @param stan_data A list of data for Stan
#' #'
#' #' @return A string of Stan code
#' #'
#' #' @export
#' #'
#' get_meta_stan_prev <- function(
#'     stan_data
#' ) {
#'
#'   return(
#'     "data {
#'          int<lower=0> N;
#'          int<lower=0> K;
#'          int<lower=0> n_ests;
#'          real<lower=0> alpha_prior_mean[N];
#'          real<lower=0> alpha_prior_sd[N];
#'          real<lower=0> bias_prior_mean[K];
#'          real<lower=0> bias_prior_sd[K];
#'          int<lower=1> study[n_ests];
#'          int<lower=1> design[n_ests];
#'          vector<lower=0>[n_ests] ests;
#'          vector<lower=0>[n_ests] ses;
#'        }
#'        parameters {
#'          vector<lower=0.01>[K] bias;
#'          vector<lower=0,upper=1>[N] alpha;
#'        }
#'
#'        model {
#'
#'          target += normal_lpdf(ests | bias[design].*alpha[study], ses);
#'          alpha ~ normal(alpha_prior_mean, alpha_prior_sd);
#'          bias ~ lognormal(log(bias_prior_mean), bias_prior_sd);
#'
#'        }"
#'   )
#' }
