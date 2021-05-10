#' Get individual study estimands
#'
#' @param data pass-through population data frame
#' @param known_pattern character string containing regular expression to match known group names in the study population dataset
#' @param hidden_pattern character string containing regular expression to match hidden group name in the study population dataset
#'
#' @return Estimands data frame for single study
#' @export
#'
#' @import dplyr
#' @importFrom magrittr `%>%`
#' @importFrom purrr map_int
get_study_estimands <- function(data, known_pattern = "^known", hidden_pattern = "^hidden$") {

  data %>%
    dplyr::mutate(
      degree = purrr::map_int(links, ~ length(strsplit(.x, split = ";")[[1]])),
      degree_hidden =
        purrr::map_int(links,
                       ~ sum(unlist(data[data$name %in% strsplit(.x, split = ";")[[1]],
                                         grep(pattern = hidden_pattern, names(data))])))
    ) %>%
    {
      dplyr::bind_cols(
        dplyr::summarise_at(., dplyr::vars(dplyr::matches(known_pattern),
                                           dplyr::matches(hidden_pattern),
                                           -dplyr::contains("visible")),
                            list(size = sum), na.rm  = TRUE),
        dplyr::summarise_at(., dplyr::vars(dplyr::matches(known_pattern),
                                           dplyr::matches(hidden_pattern),
                                           -dplyr::contains("visible")),
                            list(prev = mean), na.rm  = TRUE),
        dplyr::summarise_at(., dplyr::vars(degree, degree_hidden,
                                           -dplyr::contains("visible")),
                            list(average = mean), na.rm  = TRUE)
      )
    } %>%
    {
      data.frame(
        inquiry_label = names(.),
        estimand = unname(t(.)),
        stringsAsFactors = FALSE)
    } %>%
    {.[!grepl(pattern = known_pattern, .$inquiry_label),]}

}
