#' Bag of Little Bootstraps (Discretized).
#'
#' \code{BLB_Discretize}  Bag of Litle Bootstraps (BLB) for data sets with 1-dim
#' covariate. Used to generate  confidence bands for a quantile smoothing
#' splines fitted with the rqss function from package \code{\link{quantreg}}.
#' What's special about this function is that it discretizes the data set to
#' decrease the sample size, and it utilizes ALL the observations.
#'
#' This function is based on the "Bag of Little Bootstraps" (BLB) method by
#' Kleiner et al.(2012). It performs parallelization to speed up the
#' calculation.
#'
#' @import quantreg
#' @import doParallel
#' @import foreach
#'
#' @param cores The number of cores to use for parallel execution. If not
#'   specified, the number of cores is set to the value of
#'   \code{options("cores")}, if specified, or to one-half the number of cores
#'   detected by the \code{\link{parallel}} package.
#' @param data A 2-dim optional data frame (or object coercible by
#'   \code{\link{as.data.frame}} to a data frame) containing the variables in
#'   the model. The column names should be specified in a way that “x” is for
#'   the predictor, and “y” for the response.
#' @param alpha The confidence level required. The default is 0.05.
#' @param tau A specific quantile to be estimated. Must be a number between 0
#'   and 1.
#' @param m A numeric value that controls how fine that data set should be
#'   discretized.
#' @param lambda The smoothing parameter governing the tradeoff between fidelity
#'   and the penalty component for the triogram term. If \code{Search=TRUE},
#'   there is no need for users to specify a value.
#' @param D A number that specifies for how many x values you want to compute
#'   confidence intervals; the confidence band is made of pointwise intervals at
#'   varius x values. If specified, it will look at \eqn{D} equidistant points.
#' @param b The subsample size in the BLB algorithm. Kleiner et al. suggest that
#'   the size should be around \eqn{n^0.6}, where \eqn{n} is the data size. It
#'   is better to have \eqn{b*s=nrow(data)} so the asymptotic property is
#'   achieved.
#' @param s The number of subsamples used in BLB algorithm. The default value is
#'   the number of cores the user assigned for the function.
#' @param r The number of bootstrap iterations (samples with with replacement).
#'   \eqn{r=100} is suggested.
#' @param Range A vector with 2 values that specifys the range of the data set
#'   where the user wants to perform BLB over. It is defined using lower and
#'   upper quantile. The default value is \eqn{(0.001,0.999)}.
#' @param Search If \code{TRUE} (which is recommended), then the function will
#'   first search for an optimum smoothing parameter \eqn{\lambda}.
#' @return A list with three parts:
#'
#'   1. \code{x0} and \code{CI_average}, where \code{x0} contains the x values
#'   at which the confidence intervals are evaluated, and \code{CI_average} is
#'   2-dim matrix which contains the corresponding lower and upper bounds.
#'
#'   2. \code{lambda}, which is the optimum smoothing parameter selected by
#'   \code{BLB_Discretize}. If it is done automatically, the function also
#'   returns \code{Lambda} and \code{Fid}, which respectively stand for a vector
#'   of lambda values and their corresponding cross-validation MCV values.
#'
#' @examples
#' data(one)
#' result<-BLB_Discretize(cores=7,data=one, alpha=0.05, tau=0.5, Search=TRUE)
#'
#' plot<-FitPlot_1dim(data=one, result=result, xlab='x',ylab='y')
#' plot
#'
#'
#' @references Micheal, I. J et al.(2012). \eqn{A Scalable Bootstrap for Massive
#'   Data}.
#' @seealso \code{\link{BLB_Discretize_2dim}} for BLB with 2-dim covariate
#'   dataset.
#'
#'@export


BLB_Discretize<-function(cores=NULL,data,m=80,alpha=0.05,tau=0.25,lambda=2,D=100,s=cores, b=floor(nrow(data)/s), r=100, Range=c(0.001,0.999),Search=FALSE){
  if(is.null(cores)) {cores=round(0.5*detectCores())}
  registerDoParallel(cores=cores)

  if(!all(colnames(data)==c("x","y"))) colnames(data)<-c("x","y")

  ###1. Sample s disjoint subsets with size b
  n<-nrow(data)
  Tmp<-Discretize(dats=data, m=m)
  refer<-rep(1:nrow(Tmp),Tmp$Weights)
  indices<-matrix(head(refer,b*s),nrow=b)
  Rboot<-rmultinom(s,size=n,prob=rep(1/b,b))
  indicator<-(Search|length(lambda)>1)


  ###2. λ Search
  if(indicator){
    {if(length(lambda)>1) Lambda<-lambda else Lambda<-c(10,20,100,500,1000,2000,3000,5000,9000,seq(10000,70000,by=10000),100000)
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

    #plot(Lambda[-1],rowMeans(Fid)[-1],type='l',ylab="New Metric",xlab="λ",main="MCV for BLB repeated data")
    Fidelity<-rowMeans(Fid)
    tmp<-abs(diff(Fidelity))/Fidelity[1]
    Ins<-median(tmp)
    which<-which(tmp<Ins)
    if(length(which)>0){lambda<-Lambda[min(which+1,length(Fidelity))]} else lambda<-Lambda[which.min(Fidelity)]}
  }

  ###3. Bag of Little Bootstrap
  xrange<-quantile(Tmp$x,probs=Range)
  x0 <- seq(xrange[1], xrange[2], by = diff(xrange)/D)
  response <- rep(NA,length(x0))   #Prediction for one bootstrap sample
  Response<-matrix(nrow=length(x0),ncol=r)

  CIs<-foreach (j=1:s) %dopar% {
    Ind<-indices[,j]
    subsample<-Tmp[unique(Ind),]
    extrapolated = x0 > max(subsample$x) | x0 < min(subsample$x)

    for(k in 1:r){
      response[1:length(x0)] <- NA
      rboot<-as.vector(rmultinom(1,size=n,prob=rep(1/b,b)))
      bootInd<-aggregate(list(w=rboot),list(Ind=Ind),sum)
      bootData<-Tmp[bootInd$Ind,]; bootData$Weights<-NULL

      mod <- rqss_new(y ~ qss(x, lambda = lambda), tau = tau, data = bootData, weights = as.vector(bootInd$w))
      response[!extrapolated] = predict(mod, newdata = data.frame(x=x0[!extrapolated]))
      Response[,k]<-response
    }
    return(t(apply(Response,1,function(x) quantile(x, probs=c(alpha/2,1-alpha/2),na.rm=TRUE))))
  }

  #Average CI bounds computed for different data subsets
  IsNull<-rep(TRUE,length(CIs))
  for(j in 1:length(CIs)){
    if(is.null(CIs[[j]])) IsNull[j]<-FALSE
  }
  CIs<-CIs[IsNull]

  CI_average<-matrix(NA,nrow = length(x0), ncol = 2)
  for(i in 1:length(x0)){
    times<-0
    ci<-rep(0,2)
    for (j in 1:length(CIs)){
      if(!is.na(CIs[[j]][i,1]))  {times<-times+1;ci<-ci+CIs[[j]][i,]}
    }
    CI_average[i,]<-ci/times
  }

  if(indicator) return(list(x0=x0, CI_average=CI_average, lambda=lambda, Lambda=Lambda, Fid=Fidelity))
  else return(list(x0=x0, CI_average=CI_average, lambda=lambda))
}

#'\code{Discretize} is a function to discretize the data set.
#' @name Discretize
#' @rdname BLB_Discretize
#' @export

Discretize<-function(dats, m=30){
  myIQR<-function (x, na.rm = FALSE, type = 7) diff(quantile(as.numeric(x), c(0.002, 0.998), na.rm = na.rm, names = FALSE, type = type))

  scale<-apply(dats,2,myIQR)/m
  tmp<-round(mapply('/',dats,scale))
  tmp<-data.frame(apply(tmp,2,as.integer))

  Tmp<-aggregate(list(Weights=rep(1,nrow(tmp))),tmp,length)
  Tmp<-data.frame(mapply("*",Tmp,c(scale,1)))
  Tmp<-Tmp[sample(1:nrow(Tmp)),]
  return(Tmp)
}
