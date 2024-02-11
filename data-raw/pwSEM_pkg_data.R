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
-2*log(0.5545)
1-pchisq(1.1794,2)
summary(mgcv::gam(X3~X2+X1,data=sim_poisson.no.nesting,family=poisson))

# Case 1 of Appendix 1 of MEE paper
set.seed(99)
N<-c(10,50,100,200,500,1000)
nsim<-5000
p<-matrix(NA,nrow=nsim,ncol=6)
for(j in 1:nsim){
  for(i in 1:6){
    X2<-rnorm(N[i])
    X1<-0.5*X2+rnorm(N[i],0,sqrt(1-0.5^2))
    X3<-0.5*X2+rnorm(N[i],0,sqrt(1-0.5^2))
    r1<-residuals(mgcv::gam(X1~X2))
    r2<-residuals(mgcv::gam(X3~X2))
    p[j,i]<-generalized.covariance(r1,r2)$prob
  }
}
for(i in 1:6){
    print(round(quantile(p[,i],probs=c(0.025,0.05,0.1,0.5)),4))
}

# Case 2 of Appendix 1 of MEE paper
#X1<--X2-->X3
set.seed(99)
X2<-rchisq(10000,2)
X1<-rpois(10000,0.5*X2)
X3<-rpois(10000,0.5*X2)
fit1<-mgcv::gam(X1~s(X2),family="poisson")
fit2<-mgcv::gam(X3~s(X2),family="poisson")
summary(mgcv::gam(X1~X2+X3,family="poisson"))
summary(mgcv::gam(X3~X2+X1,family="poisson"))
#Without smooths...
fit1<-mgcv::gam(X1~X2,family="poisson")
fit2<-mgcv::gam(X3~X2,family="poisson")
plot(residuals(fit1,type="response")~residuals(fit2,type="response"))
generalized.covariance(residuals(fit1,type="response"),residuals(fit2,type="response"))
#p=0.057
#with smooths using defaults...
fit1<-mgcv::gam(X1~s(X2),family="poisson")
fit2<-mgcv::gam(X3~s(X2),family="poisson")
plot(residuals(fit1,type="response")~residuals(fit2,type="response"))
generalized.covariance(residuals(fit1,type="response"),residuals(fit2,type="response"))
#p=0.42
#with smooths using more knots...
#k=3, p=0.55
#k=5, p=0.10
#k=10 (default), p=0.42
#k=15, p=0.48
#k=20, p=0.47
#k=30, p=0.53
#k=50, p=0.53
#k=70, p=0.52
fit1<-mgcv::gam(X1~s(X2,k=10),family="poisson")
fit2<-mgcv::gam(X3~s(X2,k=10),family="poisson")
plot(residuals(fit1,type="response")~residuals(fit2,type="response"))
generalized.covariance(residuals(fit1,type="response"),residuals(fit2,type="response"))
#p=0.42
#Using their regression as gam...
gcm.test(X=X1,Y=X3,Z=X2,regr.method="gam",plot.residuals = TRUE)
#p=0.49

set.seed(99)
N<-c(10,50,100,200)
nsim<-1000
p<-p1<-p2<-p3<-p4<-matrix(NA,nrow=nsim,ncol=length(N))
for(j in 1:nsim){
  for(i in 1:length(N)){
#testing increasing k with N
    k.max<-max(4,round(sqrt(N[i]),0))
    k.max<-min(k.max,20)
    X2<-rchisq(N[i],2)
    X1<-rpois(N[i],0.5*X2)
    X3<-rpois(N[i],0.5*X2)
#testing Lefcheck method
    fit1<-mgcv::gam(X1~X2+X3,family="poisson",control=mgcv::gam.control(maxit=1000))
    fit2<-mgcv::gam(X3~X2+X1,family="poisson",control=mgcv::gam.control(maxit=1000))
    p2[j,i]<-max(summary(fit1)$p.table[3,4],summary(fit2)$p.table[3,4],na.rm=TRUE)
    p3[j,i]<-min(summary(fit1)$p.table[3,4],summary(fit2)$p.table[3,4],na.rm=TRUE)
#testing with no smoothing
    r1<-residuals(mgcv::gam(X1~X2,family="poisson"),type="response")
    r2<-residuals(mgcv::gam(X3~X2,family="poisson"),type="response")
    p4[j,i]<-generalized.covariance(r1,r2)$prob
#testing with smoothing and varying smoother edf
    r1<-residuals(mgcv::gam(X1~s(X2,k=k.max),family="poisson"),type="response")
    r2<-residuals(mgcv::gam(X3~s(X2,k=k.max),family="poisson"),type="response")
    p[j,i]<-generalized.covariance(r1,r2)$prob
#testing with smoothing and default smoother edf
    r1<-residuals(mgcv::gam(X1~s(X2),family="poisson"),type="response")
    r2<-residuals(mgcv::gam(X3~s(X2),family="poisson"),type="response")
    p1[j,i]<-generalized.covariance(r1,r2)$prob
  }
}
#using generalized.covariance with k varying
for(i in 1:length(N)){
  print(round(quantile(p[,i],probs=c(0.025,0.05,0.1,0.5,0.975)),4))
}
#using generalized.covariance with default k
for(i in 1:length(N)){
  print(round(quantile(p1[,i],probs=c(0.025,0.05,0.1,0.5,0.975)),4))
}
#using generalized.covariance with no smoothing
for(i in 1:length(N)){
  print(round(quantile(p4[,i],probs=c(0.025,0.05,0.1,0.5,0.975)),4))
}
#using maximum prob from two regressions
for(i in 1:length(N)){
  print(round(quantile(p2[,i],probs=c(0.025,0.05,0.1,0.5,0.975),na.rm=TRUE),4))
}
#using minimum prob from two regressions
for(i in 1:length(N)){
  print(round(quantile(p3[,i],probs=c(0.025,0.05,0.1,0.5,0.975),na.rm=TRUE),4))
}
plot(p[,4]~p2[,4])



