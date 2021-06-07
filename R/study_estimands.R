#' Get individual study estimands
#'
#' @param data pass-through population data frame
#' @param known_pattern character string containing regular expression to match known group names in the study population dataset
#' @param hidden_var character string containing hidden group name in the study population data frame
#'
#' @return Estimands data frame for single study
#' @export
#'
#' @import dplyr
#' @importFrom magrittr `%>%`
#' @importFrom purrr map_int
get_study_estimands <- function(data,
                                known_pattern = "^known(\\_\\d|\\d)?$",
                                hidden_var = "hidden") {

  data %>%
    # dplyr::mutate(
    #   degree_all = purrr::map_int(links, ~ length(strsplit(.x, split = ";")[[1]])),
    #   degree_hidden =
    #     purrr::map_int(links,
    #                    ~ sum(unlist(data[data$name %in% strsplit(.x, split = ";")[[1]],
    #                                      grep(pattern = paste0("^", hidden_var, "$"),
    #                                           names(data))])))
    # ) %>%
    dplyr::summarise(
      dplyr::across(c(dplyr::matches(known_pattern),
                      dplyr::matches(paste0("^", hidden_var, "$"))),
                    list(size = sum), na.rm  = TRUE),
      dplyr::across(c(dplyr::matches(known_pattern),
                      dplyr::matches(paste0("^", hidden_var, "$"))),
                    list(prev = mean), na.rm  = TRUE),
      # dplyr::across(c(degree_all, degree_hidden), mean, na.rm  = TRUE)
      ) %>%

    {
      data.frame(
        inquiry_label = names(.),
        estimand = unname(t(.)),
        stringsAsFactors = FALSE)
    } %>%
    {.[!grepl(pattern = known_pattern, .$inquiry_label),]}

}
