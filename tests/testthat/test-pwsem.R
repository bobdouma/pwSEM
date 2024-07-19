test()# THIS IS  CODE to create the documentation for the two data
# sets to be distributed with the package
#
#' @format ## 'sim_data'
#' A simulated data set with 100 rows and 4 columns, including both Gaussian and Poisson variables,
#' but without a grouping structure requiring a mixed model context.
#' The data were generated following the following mixed acyclic graph:
#' X1->X2->X3->X4 and X2<->X4
#' \describe{
#' \item{X1}{A Guassian variable}
#' \item{X2}{A Poisson variable}
#' \item{X3}{A Poisson variable}
#' \item{X4}{A Guassian variable}
#' }
#' @source sim_data: A simulated data set
#' @format ## 'nested_data'
#' Example with nesting structure.
#' A data set with 1309 rows and 8 columns, including both Gaussian and Poisson variables,
#' and including three columns (year, nest, ind) indicating
#' the grouping structure of the data, therefore required a mixed model context.
#' \describe{
#' \item{year}{Year observation is collected}
#' \item{nest}{Nest in which the observation is collected}
#' \item{ind}{the individual within a given nest}
#' \item{XR}{A binary variable}
#' \item{XM}{A Gaussian variable}
#' \item{XH}{A Gaussian variable}
#' \item{HP}{A Guassian variable}
#' \item{XF}{A Gaussian variable}
#' }
#' @source nested_data: private
