% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/study_estimands.R
\name{get_study_estimands}
\alias{get_study_estimands}
\title{Get individual study estimands}
\usage{
get_study_estimands(
  data,
  known_pattern = "^known",
  hidden_pattern = "^hidden$"
)
}
\arguments{
\item{data}{pass-through population data frame}

\item{known_pattern}{character string containing regular expression to match known group names in the study population dataset}

\item{hidden_pattern}{character string containing regular expression to match hidden group name in the study population dataset}
}
\value{
Estimands data frame for single study
}
\description{
Get individual study estimands
}
