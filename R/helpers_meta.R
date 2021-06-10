#' Parse stan model for meta analysis
#'
#' @param stan_data list of Stan model data inputs
#'
#' @return
#' @export
get_meta_stan <- function(stan_data) {

  .sizes <-
    stan_data %>%
    names %>%
    grep(pattern = "^n[0-9]", x = .) %>%
    stan_data[.] %>%
    sapply(., function(x) x == 1)

  .K <- length(sizes)

  .dims <- ifelse(.sizes, ";", paste0("[n", 1:.K, "];"))

  return(
    paste0(
      "data {\nint<lower=0> N;\nint<lower=0> K;\n",
      paste0("int<lower=0,upper=N> n", 1:.K, ";", collapse = "\n"), "\n",
      paste0("int<lower=0,upper=N> observed", 1:.K, .dims, collapse = "\n"), "\n",
      paste0("real<lower=0> est", 1:.K, .dims, collapse = "\n"), "\n",
      paste0("real<lower=0> est_se", 1:.K, .dims, collapse = "\n"), "\n",
      "}\nparameters {\nreal error[K];\nvector<lower=0>[N] alpha;\n}\nmodel {\n",
      paste0(
        "target += normal_lpdf(est", 1:.K, "| error[", 1:.K, "] + alpha[observed", 1:.K, "], est_se", 1:.K, ");",
        collapse = "\n"),
      "\nerror ~ normal(0, 100);\nalpha ~ normal(0, 100);\n}"
    )
  )
}
