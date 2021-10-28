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
    wave_samples[[i]] <- unique(unlist(data[wave_samples[[i - 1]],"links_list"][[1]]))%>%
      setdiff(., unlist(wave_samples[1:i-1]))
  }

  return(wave_samples)
}


#' Link-tracing gibbs sampler
#'
#' @param data pass through of sammple data created in get_study_est_linktrace
#' @param strata pass through of the strata column indicator from get_study_est_linktrace
#' @param n_strata pass through of n_strata variable created in get_study_est_linktrace indicating number of unique strata
#' @param total pass through of total variable from get_study_est_linktrace
#' @param chain_samples pass through of chain_samples variable from get_study_est_linktrace indicating number of samples for gibbs sampler
#' @param chain_burnin pass through of chain_burnin variable from get_study_est_linktrace indicating number of burnin samples for gibbs sampler
#' @param priors pass through of priors from get_study_est_linktrace
#' @param param_init pass through of initial paramters list for population size, lambda and beta from get_study_est_linktrace
#' @return list of gibbs sampler estimates for population size, lambda and beta
#'
#' @export

lt_gibbs <- function(data, strata, n_strata, total, chain_samples, chain_burnin, priors, param_init){

  data_p <- data
  data_p_waves <- hiddenmeta::lt_permute(data = data_p)

  l <- matrix(0,chain_samples,n_strata)
  b <- array(0, c(n_strata, n_strata, chain_samples))
  n <- matrix(0, chain_samples, 1)

  l[1,] <- param_init[["l_0"]]
  b[,,1] <- param_init[["b_0"]]
  n[1,1] <- param_init[["n_0"]]

  prior_n <- priors[["p_n"]]
  prior_l <- priors[["p_l"]]
  prior_b <- priors[["p_b"]]

  for(t in 2:chain_samples){

    ## generate new N
    strata_count <- data%>%
      dplyr::filter(., name %in% unlist(data_p_waves[1:length(data_p_waves)-1]))%>%
      group_by_at(strata)%>%
      summarise(n = n())%>%
      .[["n"]]

    no_link = rep(1, n_strata)
    for(i in 1:n_strata){
      for(j in 1:n_strata){
        no_link[i] <- no_link[i] * (1-b[,,t-1][j,i])^strata_count[j]
      }

    }

    no_link_1 <- l[t-1,] * no_link
    no_link <- sum(no_link_1)

    n_post_range <- c(length(unlist(data_p_waves)):(total * 5))

    nn_0 <- length(data_p_waves[[1]])
    nn <- length(unlist(data_p_waves))

    n[t,1] <- sapply(n_post_range, function(k) sum(log((n_post_range[k] + 1 - nn):(n_post_range[k] - nn_0))) +
                       (n_post_range[k] - nn) * log(no_link) - prior_n * log(n_post_range[k]))%>%
      - max(.)%>%
      exp(.)%>%
      sample(n_post_range, 1, prob = .)


    ## assign strata to non sampled units
    not_sampled <- setdiff(c(1:n[t,1]), unlist(data_p_waves))
    stratum <- rep(NA, n[t,1])
    stratum[unlist(data_p_waves)] <- data_p[unlist(data_p_waves),strata]
    stratum[not_sampled] <- sample(1:n_strata, length(not_sampled), replace = FALSE, prob = no_link_1/no_link)


    ## assign links between units not in first w-1 waves
    y <- data_p %>%
      {
        nodes <- data.frame(c(1:n[t,1]))

        edges <- tidyr::drop_na(dplyr::select(., c(rds_from, name)))%>%
          rbind(.,
                t(combn(setdiff(1:n[t,1], unlist(data_p_waves[1:length(data_p_waves)-1])), 2)) %>%
                  {
                    n_pairs <- dim(.)[1]
                    link_prob <- runif(n_pairs)

                    as.data.frame(.)%>%
                      dplyr::slice(., which(link_prob < b[,,t-1][cbind(stratum[.[1:n_pairs,1]],
                                                                       stratum[.[1:n_pairs,2]])]))%>%
                      dplyr::rename("rds_from" = V1, "name" = V2)

                  }
                )

        igraph::graph_from_data_frame(edges, nodes, directed = FALSE) %>%
          igraph::as_adjacency_matrix()%>%
          as.matrix()
      }


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
                               choose(strata_count,2) - strata_link_count[i,i] + prior_b)
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



#' Link tracing gibbs sampler wrapper
#'
#' @param n_samples pass through of n_samples variable from from get_study_est_linktrace
#' @param data pass through of sammple data created in get_study_est_linktrace
#' @param n_strata pass through of n_strata variable created in get_study_est_linktrace indicating number of unique strata
#' @param chain_samples pass through of chain_samples variable from get_study_est_linktrace indicating number of samples for gibbs sampler
#' @param chain_burnin pass through of chain_burnin variable from get_study_est_linktrace indicating number of burnin samples for gibbs sampler
#' @param priors pass through of priors from get_study_est_linktrace
#' @param param_init pass through of initial paramters list for population size, lambda and beta from get_study_est_linktrace
#' @return
#'
#' @export

lt_gibbs_sample <- function(n_samples, data, n_strata, chain_samples, chain_burnin, priors, param_init){

}




