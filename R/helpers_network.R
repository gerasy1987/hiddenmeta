## Network simulation helpers

#' Handler for Simulation of Network Using Block Model
#'
#' Draw a Simulation of Network Using Block Matrix of Probabilities of Edges Based on Group Memberships
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



#' Handler for Simulation of Network Using ERGM
#'
#' Draw a Simulation of Network Using Distribution Of An Exponential Family Random Graph Model Based on Observed Data
#'
#' @param fit object of ergm class produced by fitting model using \code{ergm} function
#' @param type_function function of igraph object \code{g} that takes existing vertex attributes and transforms them into \code{type} attribute in the binary coded format (consists of 0's and 1's only)
#'
#' @return igraph network object with vertex attribute \code{type} in the binary coded format (consists of 0's and 1's only)
#' @export
#'
#' @examples
#' \dontrun{
#' library(network)
#' library(ergm)
#'
#' data('faux.magnolia.high')
#'
#' fit <- ergm(faux.magnolia.high ~
#'               edges +
#'               gwesp(0.25,fixed=T) +
#'               nodematch("Grade") + nodematch("Sex"),
#'             control = control.ergm(MCMC.interval = 10000),
#'             verbose = F)
#'
#' sim_ergm_network(fit = fit,
#'                  type_function = function(g) {
#'                    paste0( as.integer(igraph::vertex_attr(g)$Sex == "M"),
#'                            as.integer(igraph::vertex_attr(g)$Grade == 9))})
#'
#' }
#'
#' @importFrom stats simulate
#' @importFrom intergraph asIgraph
#' @importFrom igraph vertex_attr
#' @importFrom magrittr `%>%` `%<>%`
#' @importFrom plyr mapvalues
sim_ergm_network <-
  function(fit,
           type_function) {

    if (class(fit) != "ergm") stop("fit should be ergm object")

    g <-
      stats::simulate(fit, nsim = 1) %>%
      intergraph::asIgraph(.)

    igraph::vertex_attr(g) %<>%
      { c(list(type = type_function(g)), .) }

    return(g)
  }
