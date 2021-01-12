#' Simulate single population with given network structure
#'
#' @param g igraph network object with vertex attribute \code{type} in the binary coded format (consists of 0's and 1's only)
#' @param group_names named numeric vector of prevalence for each group with last group being hidden
#' @param p_visibility named list of visibility tendencies by group. This is used as mean of Beta distribution (with SD = 0.09) to generate probability of being recognized as member of group, being sampled as seed, etc. The order of objects in list have to follow the order of \code{group_names}
#' @param add_groups named list of probabilities of additional group memberships. Examples include probability of service utilization (for service multiplier), being present at particular time-location (for TLS), etc.
#'
#' @return Population data frame for single study
#' @export
#'
#' @examples
#' \dontrun{
#' get_study_population(
#'   g = sim_block_network(),
#'   group_names = c("known", "hidden"),
#'   p_visibility = list(known = .99, hidden = .7),
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
  function(g,
           group_names = c("known", "hidden"),
           p_visibility = list(known = .99, hidden = .7),
           add_groups = list(service_use = 0.3,
                             loc_1 = 0.3, loc_2 = 0.1, loc_3 = 0.2,
                             known_2 = 0.1, known_3 = 0.2)) {

    if (class(g) != "igraph") stop("g in supplied to study population should be an igraph object")

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
      dplyr::mutate(n_visible = sum(dplyr::c_across(dplyr::contains("visible_")))) %>%
      dplyr::ungroup() %>%
      dplyr::mutate()

    return(
      g %>%
        {
          dplyr::tibble(name = 1:igraph::gorder(.),
                        dplyr::as_tibble(igraph::vertex_attr(g)),
                        links = igraph::as_adj_list(.),
                        total = igraph::gorder(.))
        } %>%
        mutate_add_groups(., add_groups = add_groups) %>%
        dplyr::mutate_at(c(group_names, names(add_groups)), list(total = sum)) %>%
        dplyr::rename_at(vars(ends_with("_total")),
                         ~ paste0("total_", gsub("_total$", "", .))) %>%
        dplyr::mutate_at(
          names(add_groups),
          list(visible = ~ colSums(as.matrix((igraph::as_adj(g) * .) == 1))))
    )

  }
