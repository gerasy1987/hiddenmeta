#' Simulate populations with given network structure for multiple studies with (possibly) varying designs
#'
#' @param pop_args named list of named lists of arguments to \code{pop_network()} (the set of arguments can be excessive)
#'
#' @return
#' @export
#'
#' @import dplyr
#' @importFrom magrittr `%>%`
#' @importFrom pbapply pblapply pboptions
#' @importFrom parallel detectCores
get_pop_network <- function(pop_args) {

  pbapply::pboptions(type = "none")

  pbapply::pblapply(
    X = names(pop_args),
    cl = parallel::detectCores()/2,
    FUN = function(x) {

      do.call(what = pop_network,
              args = pop_args[[x]][names(pop_args[[x]]) %in% formalArgs(pop_network)]) %>%
        dplyr::mutate(study = x)

    }) %>%
    dplyr::bind_rows(.) %>%
    dplyr::select(study, everything())

}
