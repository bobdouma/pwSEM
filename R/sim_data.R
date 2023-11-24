#' sim_data:
#' Simulated data with correlated errors involving endogenous
#' variables and without any grouping structure.  Data generated
#' using the following mixed acyclic graph:
#' X1->X2->X3->X4 and X2<->X4
#' @format ## 'sim_data'
#' A data set with 100 rows and 4 columns
#' \describe{
#' \item{X1}{A Guassian variable}
#' \item{X2}{A Poisson variable}
#' \item{X3}{A Poisson variable}
#' \item{X4}{A Guassian variable}
#' }
"sim_data"
