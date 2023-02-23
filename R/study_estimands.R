#' Get individual study estimands
#'
#' @param data pass-through population data frame
#' @param known_pattern character string containing regular expression to match known group names in the study population dataset
#' @param hidden_var character string containing hidden group name in the study population data frame
#'
#' @return Estimands data frame for single study
#'
#' @import dplyr
#' @importFrom magrittr `%>%`
#' @importFrom purrr map_int
get_study_estimands_tidy <- function(data,
                                known_pattern = "^known(\\_\\d|\\d)?$",
                                hidden_var = "hidden") {


  known_groups <- names(data)[grepl(x = names(data), pattern = known_pattern)]

  .out <-
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
        inquiry = names(.),
        estimand = unname(t(.)),
        stringsAsFactors = FALSE)
    }# %>%
    # {.[!grepl(pattern = known_pattern, .$inquiry),]}

  if (length(known_groups) != 0) {
    for (i in seq_along(known_groups)) {
      .out <-
        data %>%
        dplyr::filter(eval(parse(text = paste0(known_groups[i], " == 1")))) %>%
        dplyr::summarise(
          dplyr::across(c(dplyr::matches(paste0("^", hidden_var, "$"))),
                        list(size = sum), na.rm  = TRUE),
          dplyr::across(c(dplyr::matches(paste0("^", hidden_var, "$"))),
                        list(prev = mean), na.rm  = TRUE),
          # dplyr::across(c(degree_all, degree_hidden), mean, na.rm  = TRUE)
        ) %>%
        {
          data.frame(
            inquiry = paste0(names(.), "_in_", known_groups[i]),
            estimand = unname(t(.)),
            stringsAsFactors = FALSE)
        } %>%
        bind_rows(.out, .)
    }
  }

  return(.out)
}



#' Get individual study estimands using data.table
#'
#' @param data pass-through population data frame
#' @param known_pattern character string containing regular expression to match known group names in the study population dataset
#' @param hidden_var character string containing hidden group name in the study population data frame
#'
#' @return Estimands data frame for single study
#' @export
#'
#' @import data.table
get_study_estimands <- function(data,
                                   known_pattern = "^known(\\_\\d|\\d)?$",
                                   hidden_var = "hidden") {

  if (!data.table::is.data.table(data)) data <- data.table::setDT(data)

  known_groups <-
    names(data)[grepl(x = names(data), pattern = known_pattern)]

  .out <-
    data[, grep(pattern = paste0(known_pattern, "|^", hidden_var, "$"),
                       x = names(data), value = TRUE)
       , with = FALSE]

  .out <-
    do.call(
      rbind,
      lapply(
        X = names(.out),
        FUN = function(x)
          data.frame(inquiry = paste0(x, c("_size", "_prev")),
                     estimand = c(sum(!is.na(.out[[x]])),
                                  mean(.out[[x]], na.rm = TRUE)))))

  if (length(known_groups) != 0) {
    for (i in seq_along(known_groups)) {
      .temp <-
        data[
          get(known_groups[i]) == 1
          , grep(pattern = paste0("^", hidden_var, "$"), names(data), value = TRUE)
          , with = FALSE
        ]

      .out <-
        rbind(.out,
              data.frame(
                inquiry = paste0(names(.temp), c("_size", "_prev"), "_in_", known_groups[i]),
                estimand = c(sum(!is.na(unlist(.temp))),
                             mean(unlist(.temp), na.rm = TRUE))))
    }
  }

  return(.out)
}
