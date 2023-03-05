#' Get individual study estimands using data.table
#'
#' @param data pass-through population data frame
#' @param in_groups character string containing regular expression to match group names for which prevalence is estimated within hidden group
#' @param out_groups character string containing regular expression to match group names within which hidden group prevalence is estimated
#' @param hidden_var character string containing hidden group name in the study population data frame
#'
#' @return Estimands data frame for single study
#' @export
#'
#' @import data.table
get_study_estimands <- function(data,
                                in_groups = NULL,
                                out_groups = "^known(\\_\\d|\\d)?$",
                                hidden_var = "hidden") {

  if (!data.table::is.data.table(data)) data <- data.table::setDT(data)

  .out_groups_vec <-
    names(data)[grepl(x = names(data), pattern = out_groups)]

  if (!is.null(in_groups)) {
    .in_groups_vec <-
      names(data)[grepl(x = names(data), pattern = in_groups)]
  }


  .out <-
    data[
      , grep(
        pattern = paste0(c(in_groups, out_groups, paste0("^", hidden_var, "$")),
                         collapse = "|"),
        x = names(data), value = TRUE)
      , with = FALSE]

  .out <-
    do.call(
      rbind,
      lapply(
        X = names(.out),
        FUN = function(x)
          data.frame(inquiry = paste0(x, c("_size", "_prev")),
                     estimand = c(sum(.out[[x]], na.rm = TRUE),
                                  mean(.out[[x]], na.rm = TRUE)))))

  # get prevalences and sizes of other groups within hidden group
  if (!is.null(in_groups)) {
    for (i in seq_along(.in_groups_vec)) {
      .temp <-
        data[get(hidden_var) == 1, .in_groups_vec[i], with = FALSE]

      .out <-
        rbind(.out,
              data.frame(
                inquiry = paste0(names(.temp), c("_size", "_prev"), "_in_", hidden_var),
                estimand = c(sum(unlist(.temp), na.rm = TRUE),
                             mean(unlist(.temp), na.rm = TRUE))))
    }
  }

  # get prevalences and sizes of hidden group within other groups
  if (!is.null(out_groups)) {
    for (i in seq_along(.out_groups_vec)) {
      .temp <-
        data[
          get(.out_groups_vec[i]) == 1
          , grep(pattern = paste0("^", hidden_var, "$"), names(data), value = TRUE)
          , with = FALSE
        ]

      .out <-
        rbind(.out,
              data.frame(
                inquiry = paste0(names(.temp), c("_size", "_prev"), "_in_", .out_groups_vec[i]),
                estimand = c(sum(unlist(.temp), na.rm = TRUE),
                             mean(unlist(.temp), na.rm = TRUE))))
    }
  }

  return(.out)
}
