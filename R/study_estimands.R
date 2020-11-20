#' Get individual study estimands
#'
#' @param data pass-through population data frame
#'
#' @return Estimands data frame for single study
#' @export
#'
#' @import dplyr
#' @importFrom magrittr `%>%`
#' @importFrom purrr map_int
get_study_estimands <- function(data) {

  data %>%
    dplyr::mutate(degree = purrr::map_int(links, ~ length(.x)),
                  degree_hidden = purrr::map_int(links, ~ sum(data$hidden[data$name %in% .x]))) %>%
    {
      dplyr::bind_cols(
        dplyr::summarise_at(., dplyr::vars(dplyr::starts_with("known"), hidden, -dplyr::contains("visible")),
                            list(size = sum), na.rm  = TRUE),
        dplyr::summarise_at(., dplyr::vars(dplyr::starts_with("known"), hidden, -dplyr::contains("visible")),
                            list(prev = mean), na.rm  = TRUE),
        dplyr::summarise_at(., dplyr::vars(degree, degree_hidden, -dplyr::contains("visible")),
                            list(average = mean), na.rm  = TRUE)
      )
    } %>%
    {
      data.frame(
        estimand_label = names(.),
        estimand = unname(t(.)),
        stringsAsFactors = FALSE)
    }

}
