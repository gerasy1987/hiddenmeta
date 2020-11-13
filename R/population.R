#' Simulate signle population with given network structure
#'
#' @param N number of units in population
#' @param K number of groups
#' @param prev_K numeric vector of prevalence for each group with last group being hidden population
#' @param rho_K numeric vector of correlations in group memberships
#' @param p_edge_within named list of numeric vectors giving probability of link between in-group members and out-group members for each of groups
#' @param p_edge_between named list of numeric values giving probability of link between in- and out-group member for each of groups
#' @param p_visibility list of probabilities of visibility (show-up)
#' @param p_service probability of service utilization in hidden population
#' @param directed logical, whether links are directed or undirected
#'
#' @return
#' @export
#'
#' @examples
#' population_network(
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
#'
#'
#' @import dplyr
#' @importFrom igraph sample_pref vertex_attr as_adj as_adj_list
#' @importFrom plyr mapvalues
#' @importFrom stringr str_split
#' @importFrom magrittr `%>%`
#' @importFrom fastDummies dummy_cols
population_network <-
  function(N, K, prev_K, rho_K, p_edge_within, p_edge_between,
           p_visibility, p_service, directed = FALSE) {

    type_names <-
      lapply(1:K, function(x) 0:1) %>%
      expand.grid() %>%
      apply(X = ., MARGIN = 1, FUN = paste0, collapse = "")

    # this can only handle 2 groups at the moment !!!
    g <-
      igraph::sample_pref(nodes = N, types = 2^K,
                          type.dist = gen_group_sizes(N = N, prev_K = prev_K, rho_K = rho_K),
                          fixed.sizes = TRUE, directed = FALSE, loops = FALSE,
                          pref.matrix= gen_block_matrix(p_edge_within = p_edge_within,
                                                        p_edge_between = p_edge_between))

    igraph::vertex_attr(g)$type <-
      igraph::vertex_attr(g)$type %>%
      { plyr::mapvalues(x = ., from = unique(.), to = type_names) }


    igraph::vertex_attr(g) <-
      igraph::vertex_attr(g)$type %>%
      stringr::str_split(., pattern = "", simplify = TRUE) %>%
      apply(X = ., MARGIN = 2, as.integer) %>%
      { `colnames<-`(., c(paste0("known_", 1:(ncol(.)-1)), "hidden"))} %>%
      dplyr::as_tibble(.) %>%
      dplyr::mutate_at(dplyr::vars(dplyr::starts_with("known_")),
                       list(known_visible = ~ rbinom(n(), 1, p_visibility$known) * .)) %>%
      dplyr::mutate(
        hidden_visible = rbinom(dplyr::n(), 1, p_visibility$hidden) * hidden,
        type = igraph::vertex_attr(g)$type,
        p_visibility_hidden = p_visibility$hidden,
        p_visibility_known = p_visibility$known) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        type_visible = paste0(dplyr::c_across(dplyr::ends_with("visible")), collapse = "")) %>%
      dplyr::ungroup() %>%
      fastDummies::dummy_cols(.data = .,
                              select_columns = "type_visible",
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
      dplyr::group_by(hidden) %>%
      dplyr::mutate(service = sample(c(rep(0, ceiling(n() * (1 - p_service))),
                                       rep(1, floor(n() * p_service))))) %>%
      dplyr::ungroup()

    return(data)

  }

#' Simulate population with given network structure for multiple studies with varying designs
#'
#' @param pop_args named list of named lists of arguments to population_network()
#'
#' @return
#' @export
#'
#' @examples
#'
#' @import tidyverse
#' @importFrom magrittr `%>%`
#' @importFrom pbapply pblapply
get_populations <- function(pop_args) {

  pbapply::pblapply(
    X = names(pop_args),
    cl = 2,
    FUN = function(x) {

      do.call(what = population_network,
              args = pop_args[[x]][names(pop_args[[x]]) %in% formalArgs(population_network)]) %>%
        dplyr::mutate(study = x)

    }) %>%
    dplyr::bind_rows(.) %>%
    dplyr::select(study, everything())

}
