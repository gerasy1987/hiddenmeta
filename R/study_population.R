#' Simulate single population with given network structure
#'
#' @param network_handler function that takes several arguments and returns igraph network object with vertex attribute \code{type} in the binary coded format (consists of 0's and 1's only)
#' @param network_handler_args list of arguments passed to \code{network_handler}
#' @param group_names character vector of groups within population with last group being hidden
#' @param p_visible named list of visibility tendencies by group. This is used as mean of Beta distribution (with SD = 0.09) to generate probability of being recognized as member of group, being sampled as seed, etc. The order of objects in list have to follow the order of \code{group_names}
#' @param add_groups named list of probabilities of additional group memberships. Examples include probability of service utilization (for service multiplier), being present at particular time-location (for TLS), etc.
#'
#' @return Population data frame for single study
#' @export
#'
#' @examples
#' \dontrun{
#'  get_study_population(
#'    network_handler = sim_block_network,
#'    network_handler_args =
#'      list(N = 2000, K = 3,
#'           prev_K = c(frame = .7,
#'                      known = .2,
#'                      hidden = .05),
#'           rho_K = c(.05, .01, .01),
#'           p_edge_within =
#'             list(frame = c(.3,.3),
#'                  known = c(.3,.3),
#'                  hidden = c(.3,.9)),
#'           p_edge_between =
#'             list(frame = 0.3,
#'                  known = 0.3,
#'                  hidden = 0.05),
#'           directed = FALSE),
#'    group_names = c("frame", "known", "hidden"),
#'    p_visible = list(frame = 1, known = 1, hidden = .7),
#'    add_groups = list(service_use = 0.3,
#'                      loc_1 = 0.3, loc_2 = 0.1, loc_3 = 0.2,
#'                      known_2 = 0.1, known_3 = 0.2))
#' }
#'
#' @import dplyr
#' @importFrom igraph vertex_attr as_adj as_adj_list
#' @importFrom Matrix rowSums
#' @importFrom fastDummies dummy_cols
#' @importFrom stringr str_split
get_study_population <-
  function(network_handler = sim_block_network,
           network_handler_args,
           group_names,
           p_visible,
           add_groups) {

    if (length(unique(group_names)) != length(group_names) |
        any(names(p_visible) != group_names))
      stop("group names have to be unique and consistent")

    if (any(grepl("visible|type|_in|_out", c(names(add_groups), group_names))))
      stop("some group names include `visible`, `type`, `_in` or `_out`.
           This might clash with internal naming conventions and produce inadequate results")

    g <- do.call(what = network_handler, args = network_handler_args)

    .g_adj <- igraph::as_adj(g)
    .g_attr <- igraph::vertex_attr(g)

    if (!("igraph" %in% class(g)))
      stop("network_handler should produce an igraph object")

    .base_df <-
      .g_attr$type %>%
      stringr::str_split(., pattern = "", simplify = TRUE) %>%
      dplyr::as_tibble(.name_repair = function(names) group_names) %>%
      dplyr::mutate(
        across(all_of(group_names), as.integer),
        type = .g_attr$type)

    .out_df <-
      .base_df %>%
      mutate_add_groups(., add_groups = add_groups)

    .add_groups_names <- setdiff(names(.out_df), names(.base_df))

    .out_df %>%
      mutate_visibility(., vars = names(p_visible), p_visible = p_visible) %>%
      dplyr::mutate(
        type_visible =
          apply(X = .[,paste0(group_names, "_visible")], MARGIN = 1, FUN = paste0, collapse = "")) %>%
      fastDummies::dummy_cols(., select_columns = "type_visible", remove_selected_columns = TRUE) %>%
      dplyr::rename_with(
        .cols = starts_with("type_visible_"),
        ~ sapply(base::strsplit(.x, split = "_"), function(x) paste0(x[c(1,3,2)], collapse = "_"))) %>%
      {
        dplyr::bind_cols(
          .,
          # out-reports given visibility of others
          `names<-`(
            as.data.frame(
              as.matrix(
                .g_adj %*% as.matrix(.[,names(.)[grepl("_visible$", names(.))]]))),
            paste0(names(.)[grepl("_visible$", names(.))], "_out")),
          `names<-`(
            as.data.frame(
              as.matrix(
                .g_adj %*% as.matrix(.[,.add_groups_names]))),
            paste0(.add_groups_names, "_visible_out")),
          # in-reports given subject's visibility
          `names<-`(
            as.data.frame(
              as.matrix(
                Matrix::rowSums(.g_adj) * as.matrix(.[,names(.)[grepl("_visible$", names(.))]]))),
            paste0(names(.)[grepl("_visible$", names(.))], "_in")),
          `names<-`(
            as.data.frame(
              as.matrix(
                Matrix::rowSums(.g_adj) * as.matrix(.[,.add_groups_names]))),
            paste0(.add_groups_names, "_visible_in"))
        )
      } %>%
      dplyr::mutate(
        n_visible_out = apply(X = .[,grep("^type_.*_out$", names(.))], MARGIN = 1, FUN = sum),
        name = 1:n(),
        links = sapply(igraph::as_adj_list(g), paste0, collapse = ";"),
        total = n(),
        across(all_of(c(group_names, .add_groups_names)),
               list(total = sum), .names = "{.fn}_{.col}")) %>%
      dplyr::bind_cols(.,
                       .g_attr[-which(names(.g_attr) == "type")]) %>%
      dplyr::select(name, type, all_of(group_names), links, # main parameters
                    .add_groups_names, # additional group memberships
                    n_visible_out,
                    ends_with("_visible_out"), # out-reports
                    ends_with("_visible_in"), # in-reports
                    starts_with("p_visible_"), # probabilities of visible
                    starts_with("total")) # total numbers by group
  }
