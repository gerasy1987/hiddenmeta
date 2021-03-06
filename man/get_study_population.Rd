% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/study_population.R
\name{get_study_population}
\alias{get_study_population}
\title{Simulate single population with given network structure}
\usage{
get_study_population(
  network_handler = sim_block_network,
  network_handler_args = list(N = 2000, K = 2, prev_K = c(known = 0.3, hidden = 0.1),
    rho_K = 0.05, p_edge_within = list(known = c(0.05, 0.05), hidden = c(0.05, 0.9)),
    p_edge_between = list(known = 0.05, hidden = 0.01), directed = FALSE),
  group_names = c("known", "hidden"),
  p_visibility = list(known = 0.99, hidden = 0.7),
  add_groups = list(service_use = 0.3, loc_1 = 0.3, loc_2 = 0.1, loc_3 = 0.2, known_2 =
    0.1, known_3 = 0.2)
)
}
\arguments{
\item{network_handler}{function that takes several arguments and returns igraph network object with vertex attribute \code{type} in the binary coded format (consists of 0's and 1's only)}

\item{network_handler_args}{list of arguments passed to \code{network_handler}}

\item{group_names}{named numeric vector of prevalence for each group with last group being hidden}

\item{p_visibility}{named list of visibility tendencies by group. This is used as mean of Beta distribution (with SD = 0.09) to generate probability of being recognized as member of group, being sampled as seed, etc. The order of objects in list have to follow the order of \code{group_names}}

\item{add_groups}{named list of probabilities of additional group memberships. Examples include probability of service utilization (for service multiplier), being present at particular time-location (for TLS), etc.}
}
\value{
Population data frame for single study
}
\description{
Simulate single population with given network structure
}
\examples{
\dontrun{
get_study_population(
  network_handler = sim_block_network,
  network_handler_args =
    list(N = 2000, K = 2, prev_K = c(known = .3, hidden = .1), rho_K = .05,
         p_edge_within = list(known = c(0.05, 0.05), hidden = c(0.05, 0.9)),
         p_edge_between = list(known = 0.05, hidden = 0.01),
         directed = FALSE),
  group_names = c("known", "hidden"),
  p_visibility = list(known = .99, hidden = .7),
  add_groups = list(service_use = 0.3,
                    loc_1 = 0.3, loc_2 = 0.1, loc_3 = 0.2,
                    known_2 = 0.1, known_3 = 0.2))
}

}
