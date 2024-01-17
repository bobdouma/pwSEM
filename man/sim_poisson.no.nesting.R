#' sim_poisson.no.nesting:
#' Simulated data with correlated errors involving endogenous
#' variables and without any grouping structure.  Includes Poisson distributed variables.
#' Data generated using this mixed acyclic graph:
#' X1->X2->X3->X4 and X2<->X4
#' @format ## 'sim_poisson.no.nesting'
#' A data set with 100 rows and 4 columns
#' \describe{
#' \item{X1}{A Gaussian variable}
#' \item{X2}{A Poisson variable}
#' \item{X3}{A Poisson variable}
#' \item{X4}{A Poisson variable}
#' }
"sim_poisson.no.nesting"
