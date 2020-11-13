## Helper functions

#' Generate sizes of groups in population
#'
#' @param N sample size
#' @param prev_K numeric vector of prevalence for each group with last group being hidden population
#' @param rho_K numeric vector of correlations in group memberships
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' gen_group_sizes(N = 1000,
#'                 prev_K = c(known = .3, hidden = .1),
#'                 rho_K = .05)
#' }
#'
gen_group_sizes <- function (N, prev_K, rho_K) {

  # can only handle 2 groups at the moment !!!
  {rho_K * sqrt(prod(prev_K * (1 - prev_K))) + prod(prev_K)} %>%
    round(.,floor(log10(N))) %>%
    {N * c(1 - sum(prev_K) + ., prev_K - ., .)} %>%
    `names<-`(
      apply(
        expand.grid(
          lapply(
            seq_along(prev_K), function(x) 0:1)), 1, paste0, collapse = ""))
}


#' Generate block matrix of link probabilities
#'
#' Returns full block matrix of link probabilities given the group memberships based on probabilities of links within each of the groups and between in- and out-group members
#'
#' @param p_edge_within named list of numeric vectors giving probability of link between in-group members and out-group members for each of groups
#' @param p_edge_between named list of numeric values giving probability of link between in- and out-group member for each of groups
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' gen_block_matrix(p_edge_within = list(known = c(0.05, 0.05), hidden = c(0.05, 0.9)),
#'                  p_edge_between = list(known = 0.05, hidden = 0.01))
#' }
#'
#' @importFrom magrittr `%>%`
#' @importFrom RcppAlgos permuteGeneral
gen_block_matrix <- function(p_edge_within, p_edge_between) {

  if (length(p_edge_within) != length(p_edge_between))
    stop("p_edge_within needs to be a list of length K")

  within_mat <-
    lapply(rev(p_edge_within), function(x) 1 - diag(1 - x)) %>%
    Reduce(kronecker, x = .)

  between_mat <-
    lapply(rev(p_edge_between), function(x) RcppAlgos::permuteGeneral(c(1,x))) %>%
    Reduce(kronecker, x = .)

  out_mat <- within_mat * between_mat

  colnames(out_mat) <- rownames(out_mat) <-
    lapply(1:length(p_edge_within), function(x) 0:1) %>%
    do.call(expand.grid, .) %>%
    apply(1, paste0, collapse = "")

  return(out_mat)
}

#' Random generation from Beta distribution
#'
#' Reparametrization of base R rbeta to take location and scale rather than 2 shape parameters
#'
#' @param n number of observations
#' @param mu location parameter
#' @param sd scale parameter
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' rbeta_mod(4, .5, .1)
#' }
#'
#' @importFrom stats rbeta
rbeta_mod <- function(n, mu, sd) {
  alpha <- ((1 - mu) / sd^2 - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(rbeta(n = n, shape1 = alpha, shape2 = beta))
}

merge_vertices <- function(graph_list) {
  if (!all(c("edges", "vertices") %in% names(graph_list)))
    stop("Graph list has to have edges and vertices elements.")

  graph_list %>%
    {
      left_join(
        left_join(.$edges, cbind(.$vertices, from = 1:N),  by = "from"),
        cbind(.$vertices, to = 1:N), by = "to", suffix = c("_from", "_to")
      )
    }

}

#' Generate vector of type labels
#'
#' @param k number of groups
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' gen_marginal_types(3)
#' }
#'
gen_marginal_types <- function(k) {
  sapply(0:k,
         function(x) paste0(paste0(rep("*", times = x), collapse = ""),
                            "1",
                            paste0(rep("*", times = k-x), collapse = "")))
}
