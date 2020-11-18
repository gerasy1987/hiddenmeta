#' Simulate populations with given network structure for multiple studies with (possibly) varying designs
#'
#' @param pop_args named list of named lists of arguments to \code{get_study_population()} (the set of arguments can be excessive)
#'
#' @return
#' @export
#'
#' @import dplyr
#' @importFrom magrittr `%>%`
#' @importFrom methods formalArgs
#' @importFrom pbapply pblapply pboptions
#' @importFrom parallel detectCores
get_meta_population <- function(pop_args) {

  pbapply::pboptions(type = "none")

  pbapply::pblapply(
    X = names(pop_args),
    cl = parallel::detectCores()/2,
    FUN = function(x) {

      do.call(what = get_study_population,
              args = pop_args[[x]][names(pop_args[[x]]) %in% methods::formalArgs(get_study_population)]) %>%
        dplyr::mutate(study = x)

    }) %>%
    dplyr::bind_rows(.) %>%
    dplyr::select(study, everything())

}
