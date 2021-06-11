#' Parse stan model for meta analysis
#'
#' @param stan_data list of Stan model data inputs
#'
#' @return
#' @export
get_meta_stan <- function(
  stan_data
) {

  .sizes <-
    stan_data %>%
    names %>%
    grep(pattern = "^n[0-9]", x = .) %>%
    stan_data[.] %>%
    sapply(., function(x) x == 1)

  .K <- length(.sizes)

  .dims <- ifelse(.sizes, ";", paste0("[n", 1:.K, "];"))

  return(
    paste0(
      "data {\nint<lower=0> N;\nint<lower=0> K;\n",
      paste0("int<lower=0,upper=N> n", 1:.K, ";", collapse = "\n"), "\n",
      paste0("int<lower=0,upper=N> observed", 1:.K, .dims, collapse = "\n"), "\n",
      paste0("real<lower=0> est", 1:.K, .dims, collapse = "\n"), "\n",
      paste0("real<lower=0> est_se", 1:.K, .dims, collapse = "\n"), "\n",
      "}\nparameters {\nreal<lower=0.01> rel_bias[K-1];\nvector<lower=1>[N] alpha;\n}\nmodel {\n",
      paste0(
        "target += normal_lpdf(est", 1:.K, " | ",
        c("", paste0("rel_bias[", 1:(.K - 1), "] * ")),
        "alpha[observed", 1:.K, "], est_se", 1:.K, ");",
        collapse = "\n"),
      "\nrel_bias ~ normal(1, 10);\nalpha ~ normal(300, 100);\n}" #,
      # "\ngenerated quantities {\nmatrix<lower=0>[K,K] rel_biases;\nvector<lower=0>[K] rel_bias_full;\nrel_bias_full[1] = 1;\nfor (i in 2:(K - 1)) {\nrel_bias_full[i] = rel_bias[i - 1];\n}\nfor (i in 1:K) {\nfor (j in 1:K) {\nrel_biases[i,j] = rel_bias_full[i]/rel_bias_full[j];\n}\n}\n}"
    )
  )
}


#' @export
sim_meta <- function(
  possible_samp_est =
    tibble(
      sample = list("pps", "pps",
                    "rds", "rds", "rds",
                    "tls", "tls", "tls",
                    c("rds", "pps"), c("rds", "tls")),
      estimator = c("ht", "nsum",
                    "sspse", "chords", "multiplier",
                    "ht", "nsum", "recapture",
                    "recapture", "recapture"),
      rel_bias = c(1, runif(9, 0.5, 3)),
      absolute_bias = .8 * rel_bias,
      se = (1 + (1 - rel_bias)^2) * 20),
  estimands =
    tibble(estimand = rnbinom(n = 7, size = 10, mu = 300),
           study = paste0("study_", 1:7))
) {

  tidyr::expand_grid(possible_samp_est, estimands) %>%
    dplyr::mutate(estimate = floor(estimand * absolute_bias),
                  sim_ID = 1,
                  inquiry = "hidden_size") %>%
    dplyr::arrange(inquiry, study, sample, estimator) %>%
    get_meta_sample(data = .) %>%
    dplyr::filter(meta == 1)


}
