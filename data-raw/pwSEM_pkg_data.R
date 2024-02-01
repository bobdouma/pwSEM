#MAG: X1->X2->X3->X4 and X2<->X4
#CASE 1: normal variables and no nesting
set.seed(101)
N<-100
L<-stats::rnorm(N,0,1)
X1<-stats::rnorm(N,0,1)
X2<-0.5*X1+0.5*L + stats::rnorm(N,0,sqrt(1-2*0.5^2))
X3<-0.5*X2 + stats::rnorm(N,0,sqrt(1-0.5^2))
X4<-0.5*X3+0.5*L + stats::rnorm(N,0,sqrt(1-2*0.5^2))
sim_normal.no.nesting<-data.frame(X1=X1,X2=X2,X3=X3,X4=X4)
pairs(sim_normal.no.nesting)
my.list<-list(mgcv::gam(X1~1,data=sim_normal.no.nesting,family=gaussian),
              mgcv::gam(X2~X1,data=sim_normal.no.nesting,family=gaussian),
              mgcv::gam(X3~X2,data=sim_normal.no.nesting,family=gaussian),
              mgcv::gam(X4~X3,data=sim_normal.no.nesting,family=gaussian))
out<-pwSEM(sem.functions=my.list,dependent.errors=list(X4~~X2),
           data=sim_normal.no.nesting,use.permutations = TRUE)
summary(out,structural.equations=TRUE)
usethis::use_data(sim_normal.no.nesting,overwrite=TRUE)

#CASE 2: Poisson variables and no nesting

set.seed(100)
N<-100
L<-stats::rnorm(N,0,2)
X1<-stats::rnorm(N,0,2)
lambda.X2<- 0.25*X1+0.25*L #This is the link function
X2<-rpois(n=100,lambda=exp(lambda.X2))
lambda.X3<-0.25*X2 #This is the link function
X3<-rpois(n=100,lambda=exp(lambda.X3))
lambda.X4<-0.25*X3+0.25*L #This is the link function
X4<-rpois(n=100,lambda=exp(lambda.X4))
sim_poisson.no.nesting<-data.frame(X1,X2,X3,X4)
pairs(sim_poisson.no.nesting)
my.list<-list(mgcv::gam(X1~1,data=sim_poisson.no.nesting,family=gaussian),
              mgcv::gam(X2~X1,data=sim_poisson.no.nesting,family=poisson),
              mgcv::gam(X3~X2,data=sim_poisson.no.nesting,family=poisson),
              mgcv::gam(X4~X3,data=sim_poisson.no.nesting,family=poisson))
out<-pwSEM(sem.functions=my.list,dependent.errors=list(X4~~X2),
           data=sim_poisson.no.nesting,use.permutations = TRUE)
summary(out,structural.equations=TRUE)
usethis::use_data(sim_poisson.no.nesting,overwrite=TRUE)

#CASE 3: normal variables with 2-level nesting
set.seed(101)
N<-1000
#10 groups, labeled A to J
g<-rep(c("A","B","C","D","E","F","G","H","I","J"),N/10)
g<-sort(g)
group.names<-unique(g)
group.means<-rnorm(10,4,2)
sim_normal.with.nesting<-data.frame(X1=rep(NA,N),X2=rep(NA,N),
      X3=rep(NA,N),X4=rep(NA,N),group=g)
for(i in 1:10){
  #group means follow MAG with slopes of 0.5
  L<-stats::rnorm(N/10,0,1)
  X1.g<-stats::rnorm(N/10,0,1)
  X2.g<- 0.5*L + 0.5*X1.g + stats::rnorm(N/10,0,sqrt(1-2*0.5^2))
  X3.g<- 0.5*X2.g + stats::rnorm(N/10,0,sqrt(1-0.5^2))
  X4.g<- 0.5*L + 0.5*X3.g + stats::rnorm(N/10,0,sqrt(1-2*0.5^2))
  for(j in 1:10){
    #individual values follow DAG with slopes of -0.5
    L<-stats::rnorm(N/10,0,1)
    sim_normal.with.nesting$X1[g==group.names[i]]<-stats::rnorm(N/10,X1.g[i],0.5)
    sim_normal.with.nesting$X2[g==group.names[i]]<- -1*L +
      -0.5*sim_normal.with.nesting$X1[g==group.names[i]]+stats::rnorm(N/10,X2.g[i],0.5)
    sim_normal.with.nesting$X3[g==group.names[i]]<-
      -0.5*sim_normal.with.nesting$X2[g==group.names[i]]+stats::rnorm(N/10,X3.g[i],0.5)
    sim_normal.with.nesting$X4[g==group.names[i]]<- 1*L +
      -0.5*sim_normal.with.nesting$X3[g==group.names[i]]+stats::rnorm(N/10,X4.g[i],0.5)
  }
}
pairs(sim_normal.with.nesting[,1:4])
my.list<-list(gamm4::gamm4(X1~1,random=~(1|group),data=sim_normal.with.nesting,family=gaussian),
              gamm4::gamm4(X2~X1,random=~(1|group),data=sim_normal.with.nesting,family=gaussian),
              gamm4::gamm4(X3~X2,random=~(1|group),data=sim_normal.with.nesting,family=gaussian),
              gamm4::gamm4(X4~X3,random=~(1|group),data=sim_normal.with.nesting,family=gaussian))
out<-pwSEM(sem.functions=my.list,dependent.errors=list(X4~~X2),
           data=sim_normal.with.nesting,use.permutations = FALSE,
           do.smooth=FALSE,all.grouping.vars=c("group"))
summary(out,structural.equations=TRUE)
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

#Example 1 of MEE paper
set.seed(101)
N<-500
X1<-rnorm(N,0,1)
X2<-0.5*X1+rnorm(N,0,sqrt(1-0.5^2))
X3<-0.5*X2+rnorm(N,0,sqrt(1-0.5^2))
X4<-0.5*X3+rnorm(N,0,sqrt(1-0.5^2))
dat1<-data.frame(X1,X2,X3,X4)
library(ggm)
dag1<-ggm::DAG(X2~X1,X3~X2,X4~X3)
bs1<-ggm::basiSet(dag1)
bs1
R1<-residuals(mgcv::gam(X1~X2,data=dat1))
R2<-residuals(mgcv::gam(X3~X2,data=dat1))
generalized.covariance(R1,R2)$prob
C1<- -2*sum(log(c(0.2916,0.5722,0.3302)))
1-pchisq(C1,6)
model.list<-list(mgcv::gam(X1~1,data=dat1,family=gaussian),
                 mgcv::gam(X2~X1,data=dat1,family=gaussian),
                 mgcv::gam(X3~X2,data=dat1,family=gaussian),
                 mgcv::gam(X4~X3,data=dat1,family=gaussian))
out<-pwSEM(sem.functions=model.list,data=dat1,use.permutations=FALSE)
summary(out,structural.equations=TRUE)
#Example 2 of MEE paper
mag2.1<-ggm::DAG(X2~X1,X3~X2,X4~X3,order=TRUE)
mag2.1[2,4]<-mag2.1[4,2]<-100
mag2.1
dag2<-CauseAndCorrelation::MAG.to.DAG(mag2.1)
mag2.2<-CauseAndCorrelation::DAG.to.MAG(dag2,latents="L1")
mag2.2
pwSEM::basiSet.MAG(mag2.2)
R1<-residuals(mgcv::gam(X1~X2,data=sim_normal.no.nesting))
R2<-residuals(mgcv::gam(X3~X2,data=sim_normal.no.nesting))
perm.generalized.covariance(R1,R2,nperm=10000)
-2*log(0.0813)
1-pchisq(5.0192,2)

library(lavaan)
lav.list<-'X2~X1
               X3~X2
               X4~X3
               X2~~X4'
summary(sem(lav.list,data=sim_normal.no.nesting))

#Example 3 of MEE paper
my.list<-list(mgcv::gam(X1~1,data=sim_poisson.no.nesting,family=gaussian),
              mgcv::gam(X2~X1,data=sim_poisson.no.nesting,family=poisson),
              mgcv::gam(X3~X2,data=sim_poisson.no.nesting,family=poisson),
              mgcv::gam(X4~X3,data=sim_poisson.no.nesting,family=poisson))
out<-pwSEM(sem.functions=my.list,dependent.errors=list(X4~~X2),
           data=sim_poisson.no.nesting,use.permutations = TRUE,n.perms=10000)
summary(out,structural.equations=TRUE)
R1<-residuals(mgcv::gam(X1~X2,data=sim_poisson.no.nesting,family=gaussian))
R2<-residuals(mgcv::gam(X3~X2,data=sim_poisson.no.nesting,family=poisson,type="working"))
plot(R2~R1)
perm.generalized.covariance(R1,R2,nperm=10000)
-2*log(0.165)
1-pchisq(3.6036,2)
summary(mgcv::gam(X3~X2+X1,data=sim_poisson.no.nesting,family=poisson))
