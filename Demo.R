## This R file will reproduce the examples and figures in the paper
## "Computing confidence intervals from Massive Functional Data via 
## Penalized Quantile Regression Splines" by L. Zhang, E. del Castillo, 
## A. Berglund, M. Tingley, and N. Govind
## It also demonstrates how to use the BLB_CVB method, described in the paper
## and implemented in the package ConfidenceQuant
## See also the Vignette and users's manual included with the package in folder
## \inst\doc


###---------------0. First load the package---------------
library(ConfidenceQuant)


###---------------1. Generate data set of interest-------------------
## This data set is shown in Figure 1 and is used throughout our paper
set.seed(100)
f1<-function(x) sin(2*pi*x)
tau=0.5
n<-5000
x<-runif(n,-1,1)
e0<-rnorm(n,sd=0.5)
y<-x^2+e0
sim.data<-data.frame(x=x,y=y)



###------------2. Look at the effective degrees of freedom (EDF)------------
## This section will showcase the erraticism of EDF and SIC (Figure 2).
## We will look at the median level.
tau<-0.5
Lambda<-seq(0.0005,0.1,by=0.0005)
Edf<-rep(0,length(Lambda))
Fid<-rep(0,length(Lambda))
for(i in 1:length(Lambda)){
  mod <- rqss(y ~ qss(x, lambda = Lambda[i]), tau = tau, data = sim.data)
  Fid[i]<-mod$fidelity
  Edf[i]<-mod$edf
}

plot(Lambda[-(1)],Edf[-(1)],type="l",xlab=expression(lambda),ylab=expression(p(lambda)),main="Effective Degrees of Freedom")
plot(Lambda[-(1)],(log(Fid/n)+(1/2)*(1/n)*Edf*log(n))[-(1)],type="l",xlab=expression(lambda),ylab="SIC",main=expression(paste("SIC(",lambda,") criterion")))
plot(Lambda,Fid/(n-Edf),type="l",xlab=expression(lambda),ylab="GACV",main="GACV criterion")


## Chosse two lambdas that has very different SICs.
l1<-Lambda[16]
l2<-Lambda[17]

mod1 <- rqss(y ~ qss(x, lambda = l1), tau = tau, data = sim.data)
mod2 <- rqss(y ~ qss(x, lambda = l2), tau = tau, data = sim.data)

## Duplicate Figure 1.
#----- The first lambda
plot(sim.data,col="grey")
# these are the 'zero' residual points
subset1<-which(abs(resid.rqss(mod1))<1e-5) 
points(sim.data[subset1,],pch=20,col="red")

# predict the model on a grid of x values
new.data<-data.frame(x=sort(c(seq(-0.999,0.999,by=0.0001),unique(sim.data$x))))
pred1<-predict(mod1,newdata=new.data)
lines(as.vector(t(new.data)),pred1)


#----- The second lambda
plot(sim.data,col="grey")
# these are the 'zero' residual points
subset2<-which(abs(resid.rqss(mod2))<1e-5)
points(sim.data[subset2,],pch=20,col="blue")

# predict the model on a grid of x values
new.data<-data.frame(x=sort(c(seq(-0.999,0.999,by=0.0001),unique(sim.data$x))))
pred2<-predict(mod2,newdata=new.data)
lines(as.vector(t(new.data)),pred2)


###---------------------3. MCV---------------------
## This section will demonstrate the good properties of MCV for λ selection.

## Divide the data set into 10 folds.
nfolds<-10
## Return a vector that indicates the fold number of each row.
cv.idx<-foldsGenerator(sim.data,n=nrow(sim.data),nfolds=nfolds)  

## Calculate MCV by looping through each fold and calculating prediction score
## Parallelization is applied for the iteration (takes around 4 mins)
require(doParallel)
require(foreach)
registerDoParallel(cores=3)


Lambda<-seq(0.001,4,by=0.02)
MCV<-rep(NA,length(Lambda))

for(m in 1:length(Lambda)){
  ## loop through, holding out one set of data at each point
  Fid<-foreach(i = 1:nfolds,.combine="c")%dopar%{
    hold.out.idx=which(cv.idx==i)
    keep.in.data<-sim.data[-hold.out.idx,]
    mod <- rqss(y ~ qss(x, lambda =Lambda[m]), tau = tau, data = keep.in.data)
    new.data<-sim.data[hold.out.idx,]
    extrapolated <- new.data$x > max(keep.in.data$x) | new.data$x < min(keep.in.data$x)
    
    y.pred<-as.vector(predict(mod,newdata=new.data[!extrapolated,]))
    mean(rho(tau=tau,x=new.data$y[!extrapolated]-y.pred))
  }
  
  MCV[m]<-mean(Fid)
}

## plot the MCV results and show that MCV is an easier criterion to locate minima
plot(Lambda[-1], MCV[-1],type='l',xlab=expression(lambda),ylab='MCV')
abline(v=Lambda[which.min(MCV)],lty=2,col="blue")
text(0.9,0.1982,expression(paste(lambda,"*=0.541")))


l<-Lambda[which.min(MCV)] #optimum parameter
mod <- rqss(y ~ qss(x, lambda = l), tau = tau, data = sim.data)

## Duplicate Figure 3.
plot(sim.data,col="grey")
# predict the model on a grid of x values
new.data<-data.frame(x=sort(c(seq(-0.999,0.999,by=0.0001),unique(sim.data$x))))
pred1<-predict(mod,newdata=new.data)
lines(as.vector(t(new.data)),pred1,lwd=2)
# also plot the truth
f<-function(x) x^2
curve(f, add=TRUE,col="red",lty=2,lwd=2)
legend("bottomleft",col=c("black","red"),lwd=2,lty=c(1,2),legend=c("Fitted median","True median"))



###---------------------4. MCV fail for BLB---------------------
## This section will demonstrate that MCV will fail for BLB repeated data sets.
## See algorithm 1 line 8.
n <- nrow(sim.data)
s=15; b=333 #s*b ≈ n  
shuffle <- sample(1:n)

# divide into s bags (each row of the following matrix is one bag)
indices <- matrix(head(shuffle, b * s), nrow = s, byrow = TRUE)
# resample the first bag so that it will have the original size
rboot1 <- as.vector(rmultinom(1, size = n, prob = rep(1/b,b)))
subsample <- sim.data[indices[1, ], ] # the first bag
index<-rep(1:b,rboot1)
index<-sample(index)
newdata<-subsample[index,] # the first bag with size n after resampling

nfolds<-10
cv.idx<-foldsGenerator(newdata,nfolds=nfolds)

# calculate the MCV for the BLB repeated data
registerDoParallel(cores=3)

Lambda<-seq(0.001,4,by=0.02)
MCV<-rep(NA,length(Lambda))

for(m in 1:length(Lambda)){
  ## loop through, holding out one set of data at each point
  Fid<-foreach(i = 1:nfolds,.combine="c")%dopar%{
    hold.out.idx=which(cv.idx==i)
    keep.in.data<-newdata[-hold.out.idx,]
    mod <- rqss(y ~ qss(x, lambda =Lambda[m]), tau = tau, data = keep.in.data)
    new.data<-newdata[hold.out.idx,]
    extrapolated <- new.data$x > max(keep.in.data$x) | new.data$x < min(keep.in.data$x)
    
    y.pred<-as.vector(predict(mod,newdata=new.data[!extrapolated,]))
    mean(rho(tau=tau,x=new.data$y[!extrapolated]-y.pred))
  }
  
  MCV[m]<-mean(Fid)
}

# Figure 4
plot(Lambda, MCV,type='l',xlab=expression(lambda),ylab='MCV')





###-------------------5. Apply our BLB_CVB methodology-------------------
## A grid of lambda values for which we calculate the CVB values.
## The resolution decreases for larger lambda values.
Lambda_candidates = c(0,0.05,0.1,0.15,0.2,0.3,seq(1, 30, by = 0.2),seq(30,100,by=1.5))

## A function that implements the methodology (takes around 4 mins).
## help(BLB_rqss) for the help page to the function
re<-BLB_rqss(cores=3,data=sim.data,tau=0.5,b=333,s=15,D=200,lambda = Lambda_candidates)

## Or just reload the results if you want to save time.
data(re)

## Figure 6 (CVB)
plot(re$Lambda[-(1:5)],re$Fid[-(1:5)],xlab=expression(lambda),ylab='CVB',type='l')
abline(v=6.4,lty=2,col="blue")
text(12,0.206,expression(paste(lambda,"*=6.4")))

## Figure 6 (BLB confidence bands)
library(scales)
sc<-1:201
plot(sim.data,col="grey",xlab='x',ylab='y')
ci.poly<-cbind(c(re$x0[sc],rev(re$x0[sc])),c(re$CI_average[sc, 1],rev(re$CI_average[sc, 2])))
polygon(ci.poly, col=alpha("dodgerblue3",0.6), lty=0)
lines(re$x0,re$x0^2,lwd=2)
legend("bottomleft",lwd=c(2,8),col=c("black",alpha("dodgerblue3",0.6)),legend=c("True median","Confidence bands"))

