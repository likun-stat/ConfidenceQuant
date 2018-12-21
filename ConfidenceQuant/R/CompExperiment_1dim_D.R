#'Comparative Experiment (Discretized).
#'
#'\code{CompExperiment_1dim_D} compares the confidence bands at a given quantile
#'for two different datasets, one related to a treatment and the other to a
#'control. It applies the BLB algorithm to each dataset to get confidence bands
#'using quantile smoothing splines for 1-dim covariate. What's special about
#'this function is that it discretizes the data set to decrease the sample size,
#'and it utilizes ALL the observations.
#'
#'This function runs \code{BLB} twice, once for each dataset. It is based on
#'\code{\link{BLB_Discretize}}, which implements BLB for quantile smoothing splines
#'with a one dimensional covariate dataset. It performs parallelization to speed
#'up the calculation.
#'
#'\if{html}{\figure{comp.png}{options: width=100 alt="Image output"}}
#'\if{latex}{\figure{comp.png}{options: width=3in}}
#'
#'\code{\link{CompPlot_1dim}} takes the results and use ggplot/plotly to visualize
#'them, in which different colors represent different scenarios. See figure
#'above.
#'
#'@import quantreg
#'@import doParallel
#'@import foreach
#'
#'@param cores The number of cores to use for parallel execution. If not
#'  specified, the number of cores is set to the value of
#'  \code{options("cores")}, if specified, or to one-half the number of cores
#'  detected by the \code{\link{parallel}} package.
#'@param treatment A 2-dim optional data frame for the treatment group (or
#'  object coercible by \code{\link{as.data.frame}} to a data frame) containing
#'  the variables in the model. The last column of this data frame must be the
#'  response for the experiment.
#'@param control A 2-dim optional data frame for the control group (or object
#'  coercible by \code{\link{as.data.frame}} to a data frame) containing the
#'  variables in the model. The last column of this data frame must be the
#'  response for the experiment.
#'@param alpha The confidence level required. The default is 0.05.
#'@param tau A specific quantile to be estimated. Must be a number between 0 and
#'  1.
#'@param lambda The smoothing parameter used for \code{treatment} &
#'  \code{control} if \eqn{Search=FALSE}, which governs the tradeoff between
#'  fidelity and the penalty component for the triogram term.
#'@param D A number that specifies for how many x values you want to compute
#'  confidence intervals; the confidence band is made of pointwise intervals at
#'  varius x values. If specified, it will look at \eqn{D} equidistant points.
#'@param b1 The subsample size in the BLB algorithm for \code{treatment}.
#'  Kleiner et al. suggest that the size should be around \eqn{n1^0.6}, where
#'  \eqn{n1} is the data size for \code{treatment}.
#'@param b2 The subsample size in the BLB algorithm for \code{control}. It is
#'  also suggested that the size should be around \eqn{n2^0.6}, where \eqn{n2}
#'  is the data size for \code{control}.
#'@param s The number of subsamples used in the BLB algorithm. Kleiner et al.
#'  suggest that \eqn{s} should be 10~20.
#'@param r The number of bootstrap iterations (samples with with replacement).
#'  \eqn{r=100} is suggested.
#'@param Range A vector with 2 values that specifys the range of the data set
#'  where the user wants to perform BLB over. It is defined using lower and
#'  upper quantile. The default value is \eqn{(0.001,0.999)}.
#'@param M A numeric value that controls how fine that data set should be
#'  discretized.
#'@param Search If \code{TRUE} (which is recommended), then the function will
#'  first search for an optimum smoothing parameter \eqn{\lambda}.
#'
#'@return A list with three parts - \code{result1}, \code{result2}, and
#'  \code{Diff}, which respectively return confidence bands for
#'  \code{treatment}, \code{control} and the difference between the two dataset.
#'  Each part includes the following:
#'
#'  1. \code{x0} and \code{CI_average}, where \code{x0} contains the x values at
#'  which the confidence intervals are evaluated, and \code{CI_average} is 2-dim
#'  matrix which contains the corresponding lower and upper bounds.
#'
#'  2. \code{lambda}, which is the optimum smoothing parameter selected by
#'  \code{BLB_Discretize}. If it is done automatically, the function also returns
#'  \code{Lambda} and \code{Fid}, which respectively stand for a vector of lambda
#'  values and their corresponding cross-validation MCV values.
#'
#' @examples
#' data(treatment)
#' data(control)
#'
#' #alpha=0.05;tau=0.5
#' all<-CompExperiment_1dim_D(cores=7,treatment,control,tau=0.5,Search=TRUE)
#'
#' plot<-CompPlot_1dim(treatment = treatment, control=control, all = all, xlab = 'x', ylab = 'y')
#'
#'@references Kleiner, I. J et al. JRSS B, 2012. \eqn{A Scalable Bootstrap for
#'  Massive Data}.
#'@references Akima, H. (1978). \eqn{A Method of Bivariate Interpolation and
#'  Smooth Surface Fitting for Irregularly Distributed Data Points}. ACM
#'  Transactions on Mathematical Software 4, 148-164.
#'@seealso \code{\link{contour},\link{image}}
#'@seealso \code{\link{BLB_rqss}} for BLB with one dataset that has 1-dim
#'  covariate.
#'@seealso \code{\link{CompExperiment_2dim}} for comparative experiments with
#'  2-dim covariate data sets.
#'
#'@export

CompExperiment_1dim_D<-function(cores=NULL,treatment,control,alpha=0.05,tau=0.25,lambda=2,D=100,s=15,b1=floor(nrow(treatment)/s),b2=floor(nrow(control)/s),r=100, M=80,Range=c(0.001,0.999),Search=FALSE){
  if(is.null(cores)) {cores=round(0.5*detectCores())}
  registerDoParallel(cores=cores)

  if(!all(colnames(treatment)==c("x","y"))) colnames(treatment)<-c("x","y")
  if(!all(colnames(control)==c("x","y"))) colnames(control)<-c("x","y")

  indicator<-(Search|length(lambda)>1)

  ###1. Sample s disjoint subsets with size b
  n<-nrow(treatment)
  Tmp1<-Discretize(dats=treatment,m=M)
  refer<-rep(1:nrow(Tmp1),Tmp1$Weights)
  indices1<-matrix(head(refer,b1*s),nrow=b1)
  Rboot1<-rmultinom(s,size=n,prob=rep(1/b1,b1))

  m<-nrow(control)
  Tmp2<-Discretize(dats=control,m=M)
  refer<-rep(1:nrow(Tmp2),Tmp2$Weights)
  indices2<-matrix(head(refer,b2*s),nrow=b2)
  Rboot2<-rmultinom(s,size=m,prob=rep(1/b2,b2))

  ###2. λ Search
  r1<-list(lambda=lambda)
  r2<-list(lambda=lambda)

  if(indicator){
    r1<-searchLambda_1dim_D(lambda=lambda,Tmp=Tmp1,s=s,indices=indices1,tau=tau,Rboot=Rboot1)
    r2<-searchLambda_1dim_D(lambda=lambda,Tmp=Tmp2,s=s,indices=indices2,tau=tau,Rboot=Rboot2)
  }


  ###3. Bag of Little Bootstrap
  xrange<-sort(c(quantile(Tmp1$x,probs=Range),quantile(Tmp2$x,probs=Range)))[2:3]
  x0 <- seq(xrange[1], xrange[2], by = diff(xrange)/D)

  response_trt <- rep(NA,length(x0))   #Prediction for one bootstrap sample_treatment
  Response_trt<-matrix(nrow=length(x0),ncol=r)

  response_ctr <- rep(NA,length(x0))   #Prediction for one bootstrap sample_control
  Response_ctr<-matrix(nrow=length(x0),ncol=r)

  Difference<-matrix(nrow=length(x0),ncol=r)

  CIs<-foreach (j=1:s, .inorder = FALSE) %dopar% {
    Ind1<-indices1[,j]
    subsample1<-Tmp1[unique(Ind1),]

    Ind2<-indices2[,j]
    subsample2<-Tmp2[unique(Ind2),]

    extrapolated1 = x0 > max(subsample1$x) | x0 < min(subsample1$x)
    extrapolated2 = x0 > max(subsample2$x) | x0 < min(subsample2$x)

    for(k in 1:r){
      response_trt[1:length(x0)] <- NA
      response_ctr[1:length(x0)] <- NA

      rboot1<-as.vector(rmultinom(1,size=n,prob=rep(1/b1,b1)))
      rboot2<-as.vector(rmultinom(1,size=m,prob=rep(1/b2,b2)))
      bootInd1<-aggregate(list(w=rboot1),list(Ind=Ind1),sum)
      bootInd2<-aggregate(list(w=rboot2),list(Ind=Ind2),sum)
      bootData1<-Tmp1[bootInd1$Ind,]; bootData1$Weights<-NULL
      bootData2<-Tmp2[bootInd2$Ind,]; bootData2$Weights<-NULL

      mod1 <- rqss_new(y ~ qss(x, lambda = r1$lambda), tau = tau, data = bootData1, weights = as.vector(bootInd1$w))
      response_trt[!extrapolated1] = predict(mod1, newdata = data.frame(x=x0[!extrapolated1]))
      Response_trt[,k]<-response_trt

      mod2 <- rqss_new(y ~ qss(x, lambda = r2$lambda), tau = tau, data = bootData2, weights = as.vector(bootInd2$w))
      response_ctr[!extrapolated2] = predict(mod2, newdata = data.frame(x=x0[!extrapolated2]))
      Response_ctr[,k]<-response_ctr

      Difference[,k]<-response_trt-response_ctr
    }
    trt<-t(apply(Response_trt,1,function(x) quantile(x, probs=c(alpha/2,1-alpha/2),na.rm=TRUE)))
    ctr<-t(apply(Response_ctr,1,function(x) quantile(x, probs=c(alpha/2,1-alpha/2),na.rm=TRUE)))
    dif<-t(apply(Difference,1,function(x) quantile(x, probs=c(alpha/2,1-alpha/2),na.rm=TRUE)))
    return(list(trt=trt, ctr=ctr, dif=dif))
  }

  #Average CI bounds computed for different data subsets
  IsNull<-rep(TRUE,length(CIs))
  for(j in 1:length(CIs)){
    if(is.null(CIs[[j]])) IsNull[j]<-FALSE
  }
  CIs<-CIs[IsNull]

  #1. treatment
  CI_average<-matrix(NA,nrow = length(x0), ncol = 2)
  for(i in 1:length(x0)){
    times<-0
    ci<-rep(0,2)
    for (j in 1:length(CIs)){
      if(!is.na(CIs[[j]][[1]][i,1]))  {times<-times+1;ci<-ci+CIs[[j]][[1]][i,]}
    }
    CI_average[i,]<-ci/times
  }

  if(indicator) result1<-list(x0=x0, CI_average=CI_average, lambda=r1$lambda, Lambda=r1$Lambda, Fid=r1$Fidelity)
  else result1<-list(x0=x0, CI_average=CI_average, lambda=lambda)

  #2. control
  CI_average<-matrix(NA,nrow = length(x0), ncol = 2)
  for(i in 1:length(x0)){
    times<-0
    ci<-rep(0,2)
    for (j in 1:length(CIs)){
      if(!is.na(CIs[[j]][[2]][i,1]))  {times<-times+1;ci<-ci+CIs[[j]][[2]][i,]}
    }
    CI_average[i,]<-ci/times
  }

  if(indicator) result2<-list(x0=x0, CI_average=CI_average, lambda=r2$lambda, Lambda=r2$Lambda, Fid=r2$Fidelity)
  else result2<-list(x0=x0, CI_average=CI_average, lambda=lambda)

  #3. difference
  CI_average<-matrix(NA,nrow = length(x0), ncol = 2)
  for(i in 1:length(x0)){
    times<-0
    ci<-rep(0,2)
    for (j in 1:length(CIs)){
      if(!is.na(CIs[[j]][[3]][i,1]))  {times<-times+1;ci<-ci+CIs[[j]][[3]][i,]}
    }
    CI_average[i,]<-ci/times
  }

  Diff<-list(x0=x0, CI_average=CI_average)

  return(list(result1=result1, result2=result2, Diff=Diff))
}

#'\code{searchLambda_1dim_D} is a wrapper function to calculate the optimum lambda.
#' @name searchLambda_1dim_D
#' @rdname CompExperiment_1dim_D


searchLambda_1dim_D<-function(lambda,Tmp,s,indices,tau,Rboot){
  if(length(lambda)>1) Lambda<-lambda else Lambda<-c(10,20,100,500,1000,2000,3000,5000,9000,seq(10000,70000,by=10000),100000)
  fid<-rep(NA,length(Lambda))

  Fid<-foreach(j=1:s,.combine = "cbind") %dopar% {
    Ind<-indices[,j]
    train<-aggregate(list(w=Rboot[,j]),list(Ind=Ind),sum)
    Train<-Tmp[train$Ind,]; Train$Weights<-NULL

    for(f in 1:length(Lambda)){
      mod <- rqss_new(y ~ qss(x, lambda = Lambda[f]), tau = tau, data = Train, weights = as.vector(train$w))

      fid_rec<-0
      for(p in (1:s)[-j]){
        Ind.test<-indices[,p]
        test <-aggregate(list(w=Rboot[,p]),list(Ind=Ind.test),sum)
        Test<-Tmp[test$Ind,]; Test$Weights<-NULL

        extrapolated<-(Test$x > max(Train$x) | Test$x < min(Train$x))
        y.pred<-predict(mod, newdata = Test[!extrapolated,])
        temp<-as.vector(test$w)[!extrapolated]
        fid_rec<-fid_rec+sum(temp*rho(tau=tau, Test$y[!extrapolated]-y.pred))/sum(temp)
      }
      fid[f]<-fid_rec/(s-1)
    }

    return(fid)
  }
  Fidelity<-rowMeans(Fid)
  tmp<-abs(diff(Fidelity))/Fidelity[1]
  Ins<-median(tmp)
  which<-which(tmp<Ins)
  if(length(which)>0){lambda<-Lambda[min(which+1,length(Fidelity))]} else lambda<-Lambda[which.min(Fidelity)]

  return(list(lambda=lambda,Lambda=Lambda,Fidelity=Fidelity))
}
