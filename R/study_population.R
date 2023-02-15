#' Simulate single population with given network structure using tidyverse
#'
#' @param network_handler function that takes several arguments and returns igraph network object with vertex attribute \code{type} in the binary coded format (consists of 0's and 1's only)
#' @param network_handler_args list of arguments passed to \code{network_handler}
#' @param group_names character vector of groups within population with last group being hidden
#' @param p_visible named list of visibility tendencies by group. This is used as mean of Beta distribution (with SD = 0.09) to generate probability of being recognized as member of group, being sampled as seed, etc. The order of objects in list have to follow the order of \code{group_names}
#' @param add_groups named list of probabilities of additional group memberships. Examples include probability of service utilization (for service multiplier), being present at particular time-location (for TLS), etc.
#'
#' @return Population data frame for single study
#'
#' @examples
#' \dontrun{
#'   get_study_population_tidy(
#'     network_handler = sim_block_network,
#'     network_handler_args =
#'       list(N = 1000, K = 2, prev_K = c(known = .1, hidden = .2), rho_K = 0,
#'            p_edge_within = list(known = c(0.1, 0.1), hidden = c(0.1, 0.3)),
#'            p_edge_between = list(known = 0.02, hidden = 0.02),
#'            directed = FALSE),
#'
#'     # groups
#'     group_names = c("known", "hidden"),
#'
#'     # probability of visibility (show-up) for each group
#'     p_visible = list(known = 1, hidden = 1),
#'
#'     # probability of service utilization in hidden population
#'     # for service multiplier
#'     add_groups =
#'       list(
#'         service_use = "rbinom(n(), 1, 0.25)",
#'         "purrr::map_df(hidden, ~ sapply( `names<-`(rep(0.2, times = 10),
#'         paste0('loc_', 1:10)), function(add) rbinom(length(.x), 1, 0.05 + .x * add)))",
#'         known_2 = 0.3,
#'         "purrr::map_df(known, ~ sapply( `names<-`(rep(0, times = 8),
#'         paste0('known_', 3:10)), function(add) rbinom(length(.x), 1, 0.3)))")
#'   )
#' }
#'
#' @import dplyr
#' @importFrom igraph vertex_attr as_adj as_adj_list
#' @importFrom Matrix rowSums
#' @importFrom fastDummies dummy_cols
#' @importFrom stringr str_split
get_study_population_tidy <-
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
      mutate_add_groups_tidy(., add_groups = add_groups)

    .add_groups_names <- setdiff(names(.out_df), names(.base_df))

    .out_df %>%
      mutate_visibility_tidy(., vars = names(p_visible), p_visible = p_visible) %>%
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

#' Simulate single population with given network structure using data.table
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
#'   get_study_population_dt(
#'     network_handler = sim_block_network,
#'     network_handler_args =
#'       list(N = 1000, K = 2, prev_K = c(known = .1, hidden = .2), rho_K = 0,
#'            p_edge_within = list(known = c(0.1, 0.1), hidden = c(0.1, 0.3)),
#'            p_edge_between = list(known = 0.02, hidden = 0.02),
#'            directed = FALSE),
#'
#'     # groups
#'     group_names = c("known", "hidden"),
#'
#'     # probability of visibility (show-up) for each group
#'     p_visible = list(known = 1, hidden = .5),
#'
#'     # probability of service utilization in hidden population
#'     # for service multiplier
#'     add_groups =
#'       list(
#'         service_use = "rbinom(.N, 1, 0.25)",
#'         "paste0('loc_', 1:10) := lapply(rep(.2, times = 10), function(add) rbinom(.N, 1, 0.05 + hidden * add))",
#'         known_2 = 0.3,
#'         "paste0('known_', 3:10) := lapply(3:10, function(x) rbinom(.N, 1, 0.3))")
#'   )
#' }
#'
#' @import data.table
#' @importFrom igraph vertex_attr as_adj as_adj_list
#' @importFrom Matrix rowSums as.matrix
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


    .out_df <-
      data.table::data.table(
        type = .g_attr$type,
        stringsAsFactors = FALSE
      )[, (group_names) := data.table::tstrsplit(type, "", type.convert=TRUE)]

    .type_names <- paste0("type_", unique(.out_df$type))

    .g_attr <- .g_attr[-which(names(.g_attr) == "type")]

    mutate_add_groups(.out_df, add_groups = add_groups)

    .add_groups_names <- setdiff(names(.out_df), c("type", group_names))

    .out_cols <-
      list(
        main = c("name", "type", group_names, "links", # main parameters
                 .add_groups_names),
        visible_tot = "n_visible_out",
        visible_type = c("type_visible", paste0(c(group_names, .type_names), "_visible")),
        visible_out = paste0(c(group_names, .type_names, .add_groups_names), "_visible_out"),
        visible_in = paste0(c(group_names, .type_names, .add_groups_names), "_visible_in"),
        p_visible = paste0("p_visible_", group_names),
        totals = c("total", paste0("total_", c(group_names, .add_groups_names)))
      )

    mutate_visibility(.out_df, vars = names(p_visible), p_visible = p_visible)[
      , type_visible :=
        apply(X = .SD[,paste0(group_names, "_visible"), with=FALSE],
              MARGIN = 1, FUN = paste0, collapse = "")
    ]

    .out_df[
      , (.out_cols$visible_type[(length(group_names)+2):length(.out_cols$visible_type)]) :=
        lapply(unique(.SD[["type_visible"]]), function(x) as.integer(.SD[["type_visible"]] == x))
    ][
      , (.out_cols$visible_out) :=
        apply(.g_adj %*% Matrix::as.matrix(.SD), MARGIN = 2, FUN = function(x) x,
              simplify = FALSE)
      , .SDcols = c(.out_cols$visible_type[-1], .add_groups_names)
    ][
      , (.out_cols$visible_in) :=
        apply(Matrix::rowSums(.g_adj) * Matrix::as.matrix(.SD), MARGIN = 2, FUN = function(x) x,
              simplify = FALSE)
      , .SDcols = c(.out_cols$visible_type[-1], .add_groups_names)
    ][
      , (paste0("total_", c(group_names, .add_groups_names))) :=
        lapply(.SD[, c(group_names, .add_groups_names), with=FALSE], sum)
    ][
      , `:=`(
        n_visible_out = apply(.SD, MARGIN = 1, FUN = sum),
        name = .I,
        links = sapply(igraph::as_adj_list(g), paste0, collapse = ";"),
        total = .N)
      , .SDcols = paste0(.type_names, "_visible_out")
    ]

    suppressWarnings(.out_df[, (names(.g_attr)) := .g_attr])

    data.table::setcolorder(.out_df, do.call(c, args = unname(.out_cols[-3])))

    return(.out_df)
  }
