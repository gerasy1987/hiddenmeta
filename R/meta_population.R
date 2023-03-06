#' Simulate meta study population using simulation of several studies
#'
#' @param multi_design design class object declaring full design of the set of studies with all possible sampling-estimator pairs
#' @param n_sim number of simulations to use
#' @param parallel logical. Whether to parallelize design simulations
#'
#' @return
#' @export
#'
#' @import tidyselect
#' @importFrom magrittr `%>%` `%$%`
#' @importFrom dplyr mutate filter select group_by ungroup summarize pull arrange bind_rows
#' @importFrom DeclareDesign simulate_design
#' @importFrom stringr str_extract str_remove
#' @importFrom plyr llply
#' @importFrom purrr map map_chr
get_meta_population <-
  function(multi_design,
           n_sim = 2,
           parallel = FALSE) {

    if (parallel) {

      requireNamespace(c("doParallel", "parallel"))
      doParallel::registerDoParallel(cores = parallel::detectCores() - 1)

      data <-
        plyr::llply(
          as.list(1:n_sim),
          .fun = function(x) {
            base::suppressWarnings(
              DeclareDesign::simulate_design(multi_design, sims = 1) %>%
                dplyr::mutate(
                  sim_ID = x
                )
            )
          },
          .parallel = TRUE
        ) %>%
        dplyr::bind_rows()

    } else {
      data <-
        base::suppressWarnings(
          DeclareDesign::simulate_design(multi_design, sims = n_sim)
        )
    }

    data <-
      data %>%
      dplyr::filter(!is.na(estimate)) %>%
      dplyr::mutate(
        study = stringr::str_extract(inquiry, "^study\\_(\\d+)"),
        inquiry = stringr::str_remove(pattern = paste0(study, "_"), string = inquiry),
        estimator =
          strsplit(stringr::str_remove(pattern = paste0(inquiry, "_"), string = estimator),
                   split = "_"),
        sample = purrr::map(estimator, ~ .x[1:(length(.x)-1)]),
        estimator = purrr::map_chr(estimator, ~ .x[length(.x)]),
      ) %>%
      dplyr::select(
        sim_ID, study, inquiry, sample, estimator, estimand, estimate, se
      )

    return(data)

  }
