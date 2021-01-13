## Network simulation helpers

#' Handler for Simulation of Network Using Block Model
#'
#' Simulate network using block model matrix based on the group memberships
#'
#' @param N number of units in population
#' @param K number of groups
#' @param prev_K named numeric vector of prevalence for each group with last group being hidden
#' @param rho_K numeric vector of correlations in group memberships
#' @param p_edge_within named list of numeric vectors giving probability of link between in-group members and out-group members for each of groups. The order of objects in list have to follow the order of \code{prev_K}
#' @param p_edge_between named list of numeric values giving probability of link between in- and out-group member for each of groups. The order of objects in list have to follow the order of \code{prev_K}
#' @param directed logical, whether links are directed or undirected
#'
#' @return igraph network object with vertex attribute \code{type} in the binary coded format (consists of 0's and 1's only)
#' @export
#'
#' @examples
#' \dontrun{
#' sim_block_network(
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
#'   directed = FALSE)
#' }
#'
#' @importFrom igraph sample_pref vertex_attr
#' @importFrom magrittr `%>%`
#' @importFrom plyr mapvalues
sim_block_network <-
  function(N = 2000, K = 2, prev_K = c(known = .3, hidden = .1), rho_K = .05,
           p_edge_within = list(known = c(0.05, 0.05), hidden = c(0.05, 0.9)),
           p_edge_between = list(known = 0.05, hidden = 0.01),
           directed = FALSE) {

    if (is.null(names(prev_K)) | length(unique(names(prev_K))) != K)
      stop("prev_K have to be named vector of prevalences with unique names identifying each group. The last element have to correspond to hidden group prevalence")

    type_names <- gen_group_types(K = K)

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

    return(g)

  }
