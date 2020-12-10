#' Simulate single population with given network structure
#'
#' @param N number of units in population
#' @param K number of groups
#' @param prev_K named numeric vector of prevalence for each group with last group being hidden
#' @param rho_K numeric vector of correlations in group memberships
#' @param p_edge_within named list of numeric vectors giving probability of link between in-group members and out-group members for each of groups. The order of objects in list have to follow the order of \code{prev_K}
#' @param p_edge_between named list of numeric values giving probability of link between in- and out-group member for each of groups. The order of objects in list have to follow the order of \code{prev_K}
#' @param p_visibility named list of visibility tendencies by group. This is used as mean of Beta distribution (with SD = 0.09) to generate probability of being recognized as member of group, being sampled as seed, etc. The order of objects in list have to follow the order of \code{prev_K}
#' @param add_groups named list of probabilities of additional group memberships. Examples include probability of service utilization (for service multiplier), going to particular location (for TLS), etc.
#' @param directed logical, whether links are directed or undirected
#'
#' @return Population data frame for single study
#' @export
#'
#' @examples
#' \dontrun{
#' get_study_population(
#'   # total population size for one study
#'   N = 1000,
#'   # number of groups
#'   # (K-th group is hidden population we are interested in)
#'   K = 2,
#'   # probability of membership in each of the groups (prev_K[K] is the true prevalence)
#'   prev_K = c(known = .3, hidden = .1),
#'   # correlation matrix of group memberships
#'   rho_K = .05,
#'   # block edge probabilities depending on group memberships
#'   # 1 - list of in- and out-group probability of links for each group
#'   # 2 - probability of link between in- and out-group members
#'   p_edge_within = list(known = c(0.05, 0.05), hidden = c(0.05, 0.9)),
#'   p_edge_between = list(known = 0.05, hidden = 0.01),
#'   # probability of visibility (show-up) for each group
#'   p_visibility = list(hidden = 1, known = 1),
#'   # probability of service utilization in hidden population
#'   # for service multiplier
#'   p_service = 0.3)
#' }
#'
#' @import dplyr
#' @importFrom igraph sample_pref vertex_attr as_adj as_adj_list
#' @importFrom magrittr `%>%`
#' @importFrom fastDummies dummy_cols
#' @importFrom plyr mapvalues
#' @importFrom stringr str_split
#' @importFrom stats rbinom
get_study_population <-
  function(N = 1000, K = 2, prev_K = c(known = .3, hidden = .1), rho_K = .05,
           p_edge_within = list(known = c(0.05, 0.05), hidden = c(0.05, 0.9)),
           p_edge_between = list(known = 0.05, hidden = 0.01),
           p_visibility = list(hidden = .7, known = .99),
           add_groups = list(service_use = 0.3, loc_1 = .3, loc_2 = .3,
                             loc_3 = .2), directed = FALSE) {

    type_names <- gen_group_types(K = K)

    if (is.null(names(prev_K)) | length(unique(names(prev_K))) != K)
      stop("prev_K have to be named vector of prevalences with unique names identifying each group. The last element have to correspond to hidden group prevalence")

    group_names <- names(prev_K)

    g <-
      igraph::sample_pref(nodes = N, types = 2^K,
                          type.dist = gen_group_sizes(N = N,
                                                      prev_K = prev_K, rho_K = rho_K,
                                                      .ord = type_names),
                          fixed.sizes = TRUE, directed = FALSE, loops = FALSE,
                          pref.matrix = gen_block_matrix(p_edge_within = p_edge_within,
                                                         p_edge_between = p_edge_between,
                                                         .ord = type_names))


    igraph::vertex_attr(g)$type <-
      igraph::vertex_attr(g)$type %>%
      {plyr::mapvalues(x = ., from = unique(.), to = type_names)}

    igraph::vertex_attr(g) <-
      igraph::vertex_attr(g)$type %>%
      stringr::str_split(., pattern = "", simplify = TRUE) %>%
      base::apply(X = ., MARGIN = 2, as.integer) %>%
      {`colnames<-`(., group_names)} %>%
      dplyr::as_tibble(.) %>%
      dplyr::mutate(type = igraph::vertex_attr(g)$type) %>%
      mutate_visibility(., p_visibility = p_visibility) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        type_visible = paste0(dplyr::c_across(dplyr::ends_with("visible")), collapse = "")) %>%
      dplyr::ungroup() %>%
      fastDummies::dummy_cols(., select_columns = "type_visible",
                              remove_selected_columns = TRUE) %>%
      dplyr::mutate_at(
        dplyr::vars(dplyr::contains("visible")),
        ~ colSums(as.matrix((igraph::as_adj(g) * .) == 1))) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(total_visible = sum(dplyr::c_across(dplyr::contains("visible_")))) %>%
      dplyr::ungroup()

    data <-
      g %>%
      {
        dplyr::tibble(name = 1:N,
                      dplyr::as_tibble(igraph::vertex_attr(g)),
                      links = igraph::as_adj_list(.))
      } %>%
      mutate_add_groups(., add_groups = add_groups) %>%
      dplyr::mutate_at(
        names(add_groups),
        list(visible = ~ colSums(as.matrix((igraph::as_adj(g) * .) == 1))))

    return(data)

  }
