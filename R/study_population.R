#' Simulate single population with given network structure
#'
#' @param network_handler function that takes several arguments and returns igraph network object with vertex attribute \code{type} in the binary coded format (consists of 0's and 1's only)
#' @param network_handler_args list of arguments passed to \code{network_handler}
#' @param group_names named numeric vector of prevalence for each group with last group being hidden
#' @param p_visible named list of visibility tendencies by group. This is used as mean of Beta distribution (with SD = 0.09) to generate probability of being recognized as member of group, being sampled as seed, etc. The order of objects in list have to follow the order of \code{group_names}
#' @param add_groups named list of probabilities of additional group memberships. Examples include probability of service utilization (for service multiplier), being present at particular time-location (for TLS), etc.
#'
#' @return Population data frame for single study
#' @export
#'
#' @examples
#' \dontrun{
#' get_study_population(
#'   network_handler = sim_block_network,
#'   network_handler_args =
#'     list(N = 2000, K = 2, prev_K = c(known = .3, hidden = .1), rho_K = .05,
#'          p_edge_within = list(known = c(0.05, 0.05), hidden = c(0.05, 0.9)),
#'          p_edge_between = list(known = 0.05, hidden = 0.01),
#'          directed = FALSE),
#'   group_names = c("known", "hidden"),
#'   p_visible = list(known = .99, hidden = .7),
#'   add_groups = list(service_use = 0.3,
#'                     loc_1 = 0.3, loc_2 = 0.1, loc_3 = 0.2,
#'                     known_2 = 0.1, known_3 = 0.2))
#' }
#'
#' @import dplyr
#' @importFrom igraph vertex_attr as_adj as_adj_list gorder
#' @importFrom magrittr `%>%`
#' @importFrom fastDummies dummy_cols
#' @importFrom stringr str_split
get_study_population <-
  function(network_handler = sim_block_network,
           network_handler_args =
             list(N = 2000, K = 3,
                  prev_K = c(frame = .7,
                             known = .2,
                             hidden = .05),
                  rho_K = c(.05, .01, .01),
                  p_edge_within =
                    list(frame = c(.3,.3),
                         known = c(.3,.3),
                         hidden = c(.3,.9)),
                  p_edge_between =
                    list(frame = 0.3,
                         known = 0.3,
                         hidden = 0.05),
                  directed = FALSE),
           group_names = c("frame", "known", "hidden"),
           p_visible = list(frame = 1, known = 1, hidden = .7),
           add_groups = list(service_use = 0.3,
                             loc_1 = 0.3, loc_2 = 0.1, loc_3 = 0.2,
                             known_2 = 0.1, known_3 = 0.2)) {

    g <- do.call(what = network_handler,
                 args = network_handler_args)

    if (class(g) != "igraph") stop("network_handler should produce an igraph object")

    igraph::vertex_attr(g)$type %>%
      stringr::str_split(., pattern = "", simplify = TRUE) %>%
      dplyr::as_tibble(.name_repair = function(names) group_names) %>%
      dplyr::mutate(
        across(all_of(group_names), as.integer),
        type = igraph::vertex_attr(g)$type
      ) %>%
      mutate_visibility(., vars = group_names, p_visible = p_visible) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        type_visible = paste0(dplyr::c_across(dplyr::ends_with("_visible")), collapse = "")) %>%
      dplyr::ungroup() %>%
      fastDummies::dummy_cols(., select_columns = "type_visible",
                              remove_selected_columns = TRUE) %>%
      dplyr::mutate(
        dplyr::across(c(dplyr::starts_with("type_visible"), dplyr::ends_with("_visible")),
                      ~ colSums(as.matrix((igraph::as_adj(g) * .) == 1)))) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(n_visible = sum(dplyr::c_across(dplyr::starts_with("type_visible_")))) %>%
      dplyr::ungroup() %>%
      dplyr::bind_cols(., igraph::vertex_attr(g)[-which(names(igraph::vertex_attr(g)) == "type")]) %>%
      mutate_add_groups(., add_groups = add_groups) %>%
      dplyr::mutate(
        name = 1:n(),
        links = igraph::as_adj_list(g),
        total = n(),
        across(all_of(c(group_names, names(add_groups))),
               list(total = sum), .names = "{.fn}_{.col}"),
        across(all_of(names(add_groups)),
               list(visible = ~ colSums(as.matrix((igraph::as_adj(g) * .x) == 1))),
               .names = "{.col}_{.fn}")
      ) %>%
      dplyr::select(name, type, links,
                    all_of(group_names), names(add_groups),
                    ends_with("_visible"), starts_with("type_visible"),
                    starts_with("p_visible_"), starts_with("total"))

  }
