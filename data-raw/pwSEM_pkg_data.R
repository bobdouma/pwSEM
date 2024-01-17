#MAG: X1->X2->X3->X4 and X2<->X4
#CASE 1: normal variables and no nesting
set.seed(100)
N<-100
L<-stats::rnorm(N,0,1)
X1<-stats::rnorm(N,0,1)
X2.link<-0.5*X1+0.5*L
X2<-stats::rnorm(N,X2.link,sqrt(1-2*0.5^2))
X3.link<-0.5*X2
X3<-stats::rnorm(N,X3.link,sqrt(1-0.5^2))
X4.link<-0.5*X3+0.5*L
X4<-stats::rnorm(N,X3.link,sqrt(1-2*0.5^2))
sim_normal.no.nesting<-data.frame(X1=X1,X2=X2,X3=X3,X4=X4)
pairs(sim_normal.no.nesting)
usethis::use_data(sim_normal.no.nesting,overwrite=TRUE)

#CASE 2: Poisson variables and no nesting

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
sim_poisson.no.nesting<-data.frame(X1,X2,X3,X4)
pairs(sim_poisson.no.nesting)
usethis::use_data(sim_poisson.no.nesting,overwrite=TRUE)

#CASE 3: normal variables with 2-level nesting
set.seed(100)
N<-100
#10 groups, labeled A to J
g<-rep(c("A","B","C","D","E","F","G","H","I","J"),10)
g<-sort(g)
group.names<-unique(g)
group.means<-rnorm(10,4,2)
sim_normal.with.nesting<-data.frame(X1=rep(NA,100),X2=rep(NA,100),
      X3=rep(NA,100),X4=rep(NA,100),group=g)
for(i in 1:10){
  L<-stats::rnorm(N/10,0,1)
  sim_normal.with.nesting$X1[g==group.names[i]]<-stats::rnorm(N/10,group.means[i],1)
  X2.link<-0.5*sim_normal.with.nesting$X1[g==group.names[i]]+0.5*L
  sim_normal.with.nesting$X2[g==group.names[i]]<-stats::rnorm(N/10,X2.link,sqrt(1-2*0.5^2))
  X3.link<-0.5*X2
  sim_normal.with.nesting$X3[g==group.names[i]]<-stats::rnorm(N/10,X3.link,sqrt(1-0.5^2))
  X4.link<-0.5*X3+0.5*L
  sim_normal.with.nesting$X4[g==group.names[i]]<-stats::rnorm(N/10,X3.link,sqrt(1-2*0.5^2))
}
sim_normal.with.nesting
pairs(sim_normal.with.nesting[,1:4])
usethis::use_data(sim_normal.with.nesting,overwrite=TRUE)


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
