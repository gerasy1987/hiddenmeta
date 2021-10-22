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
