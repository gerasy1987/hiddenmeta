#' Link-tracing sample permutation
#'
#' @param data pass through sample data frame
#' @return list of IDs for observations in each sampling wave
#'
#' @export


lt_permute <- function(data){

  n_initial <- length(which(data$rds_wave == 1))
  n_waves <- length(unique(data$rds_wave))

  wave_samples <- vector(mode = "list", length = n_waves)
  wave_samples[1] <- sample(data$name, n_initial)

  for(i in 2:length(wave_samples)){
    wave_samples[i] <- unique(unlist(data[wave_samples[[i - 1]],"links_list"][[1]]))%>%
      setdiff(., unlist(wave_samples[1:i-1]))
  }

  return(wave_samples)
}

#' Link-tracing MCMC
#'
#'
#'

lt_mcmc <- function(data, chain_samples, chain_burnin, prior_n){

  data_p <- data
  data_p_waves <- hiddenmeta::lt_permute(data = data_p)

  l <- matrix(0,chain_samples,n_strata)
  b <- array(0, c(n_strata, n_strata, chain_samples))
  n <- matrix(0, chain_samples, 1)

  l[1,] <- l_0
  b[,,1] <- b_0
  n[1,1] <- n_0

  for(t in 2:chain_samples){

    ## generate new N
    strata_count <- data%>%
      dplyr::filter(., name %in% unlist(data_p_waves[1:length(data_p_waves)-1]))%>%
      dplyr::count(strata)

    no_link = rep(1, n_strata)
    for(i in 1:n_strata){
      for(j in 1:n_strata){
        no_link[i] <- no_link[i] * (1-b[,,t-1][j,i])^strata_count[j]
      }

    }

    no_link_1 <- l[t-1,] * no_link
    no_link <- sum(no_link_1)

    n_post_range <- c(length(unlist(data_p_waves)):total * 5)

    nn_0 <- length(data_p_waves[[1]])
    nn <- length(unlist(data_p_waves))

    n[t,1] <- sapply(n_post_range, function(k) sum(log((n_post_range[k] + 1 - nn):(n_post_range[k] - nn_0))) +
                       (n_post_range[k] - nn) * log(no_link) - prior_n * log(n_post_range[k]))%>%
      - max(.)%>%
      exp(.)%>%
      sample(n_post_range, 1, prob = .)


    ## assign strata to non sampled units
    not_sampled <- !data_p$name %in% unlist(data_p_waves)

    data_p[not_sampled,strata] <- sample(1:n_strata, length(not_sampled), replace = FALSE, prob = no_link_1/no_link)


    ##

  }

}


#' Gibbs sampler
#'
#' @param
#' @param
#'
#' @export
#'

lt_gibbs <- function(){

}



