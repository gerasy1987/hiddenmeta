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
  wave_samples[[1]] <- sample(data$name, n_initial)

  for(i in 2:length(wave_samples)){
    wave_samples[[i]] <- sort(unique(unlist(data[data$name %in% wave_samples[[(i - 1)]],"links_list"]))%>%
                                setdiff(., unlist(wave_samples[1:(i-1)])))
  }
  return(wave_samples)
}


#' Link-tracing gibbs sampler
#'
#' @param data pass through of sammple data created in get_study_est_linktrace
#' @param y_samp pass through of adjecency matrix generated in get_study_est_linktrace
#' @param strata pass through of the strata column indicator from get_study_est_linktrace
#' @param n_strata pass through of n_strata variable created in get_study_est_linktrace indicating number of unique strata
#' @param n_waves pass through of n_waves variable created in get_study_est_linktrace indicating number of rds waves
#' @param total pass through of total variable from get_study_est_linktrace
#' @param chain_samples pass through of chain_samples variable from get_study_est_linktrace indicating number of samples for gibbs sampler
#' @param chain_burnin pass through of chain_burnin variable from get_study_est_linktrace indicating number of burnin samples for gibbs sampler
#' @param priors pass through of priors from get_study_est_linktrace
#' @param param_init pass through of initial paramters list for population size, lambda and beta from get_study_est_linktrace
#' @return list of gibbs sampler estimates for population size, lambda and beta
#'
#' @export

lt_gibbs <- function(data, y_samp ,strata, n_strata, n_waves, total, chain_samples, chain_burnin, priors, param_init){

  ## permute data
  data_p_waves <- lt_permute(data = data)
  data_p <- data

  n_p <- lapply(data_p_waves,length)
  data_p_reorder <- vector(mode = "list", length = length(data_p_waves))
  data_p_reorder[[1]] <- c(1:n_p[[1]])

  for(i in 2:length(data_p_reorder)){
    data_p_reorder[[i]] <- (sum(unlist(n_p[1:(i-1)])) + 1):(sum(unlist(n_p[1:i])))
  }

  ## initialize parameters and priors
  l <- matrix(0,chain_samples,n_strata)
  b <- array(0, c(n_strata, n_strata, chain_samples))
  n <- matrix(0, chain_samples, 1)

  l[1,] <- param_init[["l_0"]]
  b[,,1] <- param_init[["b_0"]]
  n[1,1] <- param_init[["n_0"]]

  prior_n <- priors[["p_n"]]
  prior_l <- priors[["p_l"]]
  prior_b <- priors[["p_b"]]

  ##start chain
  for(t in 2:chain_samples){

    ## generate new N
    strata_count <- data_p[unlist(data_p_waves[1:(length(data_p_waves)-1)]),strata]%>%
      table()%>%
      as.data.frame()%>%
      .[["Freq"]]

    no_link = rep(1, n_strata)
    for(i in 1:n_strata){
      for(j in 1:n_strata){
        no_link[i] <- no_link[i] * (1-b[,,t-1][j,i])^strata_count[j]
      }
    }

    no_link_1 <- l[t-1,] * no_link
    no_link <- sum(no_link_1)

    nn_0 <- length(data_p_waves[[1]])
    nn <- length(unlist(data_p_waves))

    n_post_range <- c(nn:(total * 5))

    n[t,1] <- sapply(n_post_range, function(k) sum(log((k + 1 - nn):(k - nn_0))) +
                       (k - nn) * log(no_link) - prior_n * log(k))%>%
      - max(.)%>%
      exp(.)%>%
      sample(n_post_range, 1, prob = .)


    ## assign strata to non sampled units
    not_sampled <- setdiff(c(1:n[t,1]), unlist(data_p_reorder))
    stratum <- rep(NA, n[t,1])
    stratum[unlist(data_p_reorder)] <- data_p[unlist(data_p_waves),strata]
    stratum[not_sampled] <- sample(1:n_strata, length(not_sampled), replace = TRUE, prob = no_link_1/no_link)

    link_comb <- purrr::cross2(c(1:n_waves),c(1:n_waves))%>%
      .[-length(.)]

    y <- matrix(0,n[t,1],n[t,1])

    for(i in link_comb){
      y[data_p_reorder[[i[[1]]]],data_p_reorder[[i[[2]]]]] <- y_samp[data_p_waves[[i[[1]]]],data_p_waves[[i[[2]]]]]
    }

    link_pairs <- t(combn(setdiff(1:n[t,1], unlist(data_p_reorder[1:(length(data_p_reorder)-1)])),2))
    n_pairs <- dim(link_pairs)[1]
    link_prob <- runif(n_pairs)

    assigned <- which(link_prob < b[,,t-1][cbind(stratum[link_pairs[1:n_pairs,1]], stratum[link_pairs[1:n_pairs,2]])])
    y[cbind(link_pairs[assigned,1],link_pairs[assigned,2])] <-
      y[cbind(link_pairs[assigned,2],link_pairs[assigned,1])] <- 1

    ## generate new lambda value
    strata_count <- as.data.frame(table(stratum))[["Freq"]]
    l[t,] <- as.vector(dirmult::rdirichlet(1,c(strata_count + prior_l)))


    ## generate new beta values
    strata_link_count <- matrix(0, n_strata, n_strata)
    node_pairs <- t(combn(1:n[t,1],2))
    n_pairs <- dim(node_pairs)[1]

    for(k in unique(stratum)){
      for(j in unique(stratum)){
        strata_pairs <- which(stratum[node_pairs[1:n_pairs,1]] == k & stratum[node_pairs[1:n_pairs,2]] == j)
        strata_link_count[k,j] <- sum(y[node_pairs[strata_pairs,]])
      }
    }

    b_i <- matrix(0, n_strata, n_strata)
    for(i in 1:n_strata){
      b_i[i,i] <- stats::rbeta(1,
                               strata_link_count[i,i] + prior_b,
                               choose(strata_count[i],2) - strata_link_count[i,i] + prior_b)
    }

    for(i in  1:(n_strata - 1)){
      for(j in (i+1):n_strata){
        b_i[i,j] <- b_i[j,i] <-
          stats::rbeta(1,
                       strata_link_count[i,j] + strata_link_count[j,i] + prior_b,
                       strata_count[i] * strata_count[j] - strata_link_count[i,j] - strata_link_count[j,i] + prior_b)
      }
    }

    b[,,t] <- b_i

  }

  lt_gibbs_out <- list(n = mean(n[(chain_burnin + 1):chain_samples]),
                       l = apply(l[(chain_burnin + 1):chain_samples,], 2, mean),
                       b = apply(b[,,(chain_burnin + 1):chain_samples], 1:2, mean))

  return(lt_gibbs_out)
}






