
#' Example with correlated errors involving endogenous variables
#' X1->X2->X3->X4 and X2<->X4
#' @format ## 'sim_data'
#' A data set with 100 rows and 4 columns
#' \describe{
#' \item{X1}{A Guassian variable}
#' \item{X2}{A Poisson variable}
#' \item{X3}{A Poisson variable}
#' \item{X4}{A Poisson variable}
#' }
set.seed(100)
N<-100
L<-stats::rnorm(N,0,2)
X1<-stats::rnorm(N,0,0.4)
lambda.X2<- 0.25*X1+0.25*L #This is the link function
X2<-rpois(n=100,lambda=exp(lambda.X2))
lambda.X3<-0.25*X2 #This is the link function
X3<-rpois(n=100,lambda=exp(lambda.X3))
lambda.X4<-0.25*X3+0.25*L #This is the link function
X4<-rpois(n=100,lambda=exp(lambda.X4))
sim_data<-data.frame(X1,X2,X3,X4)
pairs(sim_data)
usethis::use_data(sim_data,overwrite=TRUE)

#' Example with nesting structure
#' @format ## 'nested_data'
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
blue.tits<-read.table("C://Users//shij1401//Documents//Book A user's guide to SEM//bluetits.txt",header=T)
names(blue.tits)
nested_data<-data.frame(year=blue.tits$year,nest=blue.tits$nest,
           ind=blue.tits$ind,XR=blue.tits$recruted,
           XM=blue.tits$mass,XH=blue.tits$hemato,
           XP=blue.tits$protos,XF=blue.tits$frass)
pairs(nested_data)
usethis::use_data(nested_data,overwrite = TRUE)
