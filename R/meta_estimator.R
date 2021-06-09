get_meta_estimators = function(data, which = "hidden_size") {

  stan_data <- list(N = nrow(data),
                    K = (ncol(data)-1)/2)

  for (k in 1:stan_data$K) {
    stan_data[[ paste0("observed",k) ]] <- which(!is.na(data[,(2 * k)]))
    stan_data[[ paste0("N",k) ]] <- length(stan_data[[paste0("observed",k)]])
    stan_data[[ paste0("est",k) ]] <-
      data[stan_data[[paste0("observed",k)]],(2 * k)]
    stan_data[[ paste0("est",k,"_sd") ]] <-
      data[stan_data[[paste0("observed",k)]],(1 + 2 * k)]
  }

  fit <-
    rstan::stan(fit = stan_model_meta,
                data = stan_data,
                iter = 4000) %>%
    extract

  data.frame(estimator_label = c(paste0("prev_", 1:N)),
             estimate = c(apply(fit$alpha, 2, mean)),
             sd =   c(apply(fit$alpha, 2, sd)),
             inquiry_label = c(paste0("hidden_prev", 1:N)),
             big_Rhat = big_Rhat
  )

}
