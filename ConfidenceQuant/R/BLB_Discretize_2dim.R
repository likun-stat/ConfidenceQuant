#' Bag of Little Bootstraps (Discretized).
#'
#' \code{BLB_Discretize_2dim}  Bag of Litle Bootstraps (BLB) for data sets with
#' 2-dim covariate. Used to generate confidence bands for a quantile smoothing
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
#' @import tripack
#'
#' @param cores The number of cores to use for parallel execution. If not
#'   specified, the number of cores is set to the value of
#'   \code{options("cores")}, if specified, or to one-half the number of cores
#'   detected by the \code{\link{parallel}} package.
#' @param data A 3-dim optional data frame (or object coercible by
#'   \code{\link{as.data.frame}} to a data frame) containing the variables in
#'   the model. The column names should be specified in a way that “x1” and "x2"
#'   are for the two predictors, and “y” for the response.
#' @param m A numeric value that controls how fine that data set should be
#'   discretized.
#' @param alpha The confidence level required. The default is 0.05.
#' @param tau A specific quantile to be estimated. Must be a number between 0
#'   and 1.
#' @param D A number that specifies for how many x values you want to compute
#'   confidence intervals; the confidence surface is made of pointwise intervals
#'   at varius x values. If specified, it will compute intervals on an
#'   equidistant \eqn{D×D} grid.
#' @param lambda The smoothing parameter governing the tradeoff between fidelity
#'   and the penalty component for the triogram term. If \code{Search=TRUE},
#'   there is no need for users to specify a value.
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
#'
#' @return A list with three parts:
#'
#'   1. \code{x0} and \code{CI_average}, where \code{x0} contains the x values
#'   at which the confidence intervals are evaluated, and \code{CI_average} is
#'   2-dim matrix which contains the corresponding lower and upper bounds.
#'
#'   2. \code{lambda}, which is the optimum smoothing parameter chosen by
#'   \code{BLB_Discretize_2dim}. If it is done automatically, the function also
#'   returns \code{Lambda} and \code{Fid}, which respectively stand for a vector
#'   of lambda values and their corresponding cross-validation MCV values.
#'
#'
#' @examples
#' data(two)
#' result<-BLB_Discretize_2dim(cores=7, data=two, alpha=0.05, tau=0.5, Search=TRUE)
#'
#'
#' #1. 3D plot
#' plot<-FitPlot_2dim(data=two, result=result, xlab='x', ylab='y', zlab='z')
#' plot
#'
#' #2. Contour plot
#' #Plot the confidence bands
#' library(akima)
#' valid<-!is.na(result$CI_average[,1])
#' X<-result$x0$x1[valid];Y<-result$x0$x2[valid]
#' lower<-result$CI_average[valid,1]; upper<-result$CI_average[valid,2]
#' x1range<-range(X);x2range<-range(Y)
#' akima1<-interp(X,Y,lower,xo=seq(x1range[1],x1range[2],length=200),yo=seq(x2range[1],x2range[2],length=200))
#' akima2<-interp(X,Y,upper,xo=seq(x1range[1],x1range[2],length=200),yo=seq(x2range[1],x2range[2],length=200))
#'
#' plot(X,Y,main="Prediction Contour Plot",type="n",xlab="x1",ylab="x2")
#' image(akima1,add=TRUE,col = terrain.colors(12))  ###YOU CAN CHANGE THE NUMBER OF COLOR LEVELS!!!!
#' contour(akima1,add=TRUE,nlevels=12)
#'
#'
#' @references Kleiner, I. J et al. JRSS B, 2012. \eqn{A Scalable Bootstrap for
#' Massive Data}.
#' @references Roger Koenker.(2005). \eqn{Quantile regression}. Cambridge
#'   university press.
#' @seealso \code{\link{BLB_rqss_2dim}} for BLB without the discretization.
#'
#'@export



BLB_Discretize_2dim<-function(cores=NULL,data,m=30,alpha=0.05,tau=0.25,lambda=2,D=50, s=cores, b=floor(nrow(data)/s), r=100, Range=c(0.001,0.999), Search=FALSE){
  if(is.null(cores)) {cores=round(0.5*detectCores())}
  registerDoParallel(cores=cores)

  if(!all(colnames(data)==c("x1","x2","y"))) colnames(data)<-c("x1","x2","y")

  ###1. Sample s disjoint subsets with size b
  n<-nrow(data)
  Tmp<-Discretize(dats=data, m=m)
  refer<-rep(1:nrow(Tmp),Tmp$Weights)
  indices<-matrix(head(refer,b*s),nrow=b)
  indicator<-(Search|length(lambda)>1)

  ###2. λ Search
  if(indicator){
    if(length(lambda)>1) Lambda<-lambda else Lambda<-c(10,20,100,500,1000,2000,3000,5000,9000,seq(10000,70000,by=10000),100000)
    fid<-rep(NA,length(Lambda))

    Fid<-foreach(j=1:s,.combine = "cbind") %dopar% {
      rboot1<-as.vector(rmultinom(1,size=n,prob=rep(1/b,b)))
      rboot2<-as.vector(rmultinom(1,size=n,prob=rep(1/b,b)))

      Ind<-indices[,j]
      if(j==s) Ind.test<-indices[,1] else Ind.test<-indices[,j+1]

      train<-aggregate(list(w=rboot1),list(Ind=Ind),sum)
      test <-aggregate(list(w=rboot2),list(Ind=Ind.test),sum)
      Train<-Tmp[train$Ind,]; Train$Weights<-NULL
      Test<-Tmp[test$Ind,]; Test$Weights<-NULL

      mod <- rqss(y ~ qss(cbind(x1,x2), lambda = lambda[1], ndum=0), tau = tau, data = Train)
      tri<-tripack::tri.mesh(mod$qss[[1]][[1]][,1], mod$qss[[1]][[1]][,2])
      NONextrapolated<-tripack::in.convex.hull(tri, Test$x1, Test$x2)

      for(f in 1:length(Lambda)){
        mod <- rqss_new(y ~ qss(cbind(x1,x2), lambda = Lambda[f], ndum=0), tau = tau, data = Train, weights = as.vector(train$w))
        y.pred<-predict(mod, newdata = Test[NONextrapolated,])
        temp<-as.vector(test$w)[NONextrapolated]
        fid[f]<-sum(temp*rho(tau=tau, Test$y[NONextrapolated]-y.pred))/sum(temp)
      }
      return(fid)
    }
    Fidelity<-rowMeans(Fid)
    tmp<-abs(diff(Fidelity))/Fidelity[1]
    Ins<-median(tmp)
    which<-which(tmp<Ins)
    if(length(which)>0){lambda<-Lambda[min(which+1,length(Fidelity))]} else lambda<-Lambda[which.min(Fidelity)]
  }

  ###3. Bag of Little Bootstrap
  x1range<-quantile(Tmp$x1,probs=Range); x2range<-quantile(Tmp$x2,probs=Range)
  x0 <- expand.grid(x1=seq(x1range[1], x1range[2], by = diff(x1range)/D),x2=seq(x2range[1], x2range[2], by = diff(x2range)/D))   #Changed from 100 to 50*50
  response <- rep(NA,nrow(x0))
  Response<-matrix(nrow=nrow(x0),ncol=r)

  CIs<-foreach (j=1:s) %dopar% {
    Ind<-indices[,j]
    subsample<-Tmp[unique(Ind),]

    #x0 should not be outside of the convex hull
    mod <- rqss(y ~ qss(cbind(x1,x2), lambda = lambda[1], ndum=0), tau = tau, data = subsample)
    tri<-tripack::tri.mesh(mod$qss[[1]][[1]][,1], mod$qss[[1]][[1]][,2])
    NONextrapolated<-tripack::in.convex.hull(tri, x0$x1, x0$x2)

    for(k in 1:r){
      response[1:nrow(x0)] <- NA
      rboot<-as.vector(rmultinom(1,size=n,prob=rep(1/b,b)))
      bootInd<-aggregate(list(w=rboot),list(Ind=Ind),sum)
      bootData<-Tmp[bootInd$Ind,]; bootData$Weights<-NULL

      mod <- rqss_new(y ~ qss(cbind(x1,x2), lambda = lambda, ndum=0), tau = tau, data = bootData, weights = as.vector(bootInd$w))
      response[NONextrapolated] = predict(mod, newdata = x0[NONextrapolated,])
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

  CI_average<-matrix(NA,nrow = nrow(x0), ncol = 2)
  for(i in 1:nrow(x0)){
    times<-0
    ci<-rep(0,2)
    for (j in 1:length(CIs)){
      if(!is.na(CIs[[j]][i,1]))  {times<-times+1;ci<-ci+CIs[[j]][i,]}
    }
    if(times>0) CI_average[i,]<-ci/times
    else CI_average[i,]<-rep(NA,2)
  }

  if(indicator) return(list(x0=x0, CI_average=CI_average, lambda=lambda, Lambda=Lambda, Fid=Fidelity))
  else return(list(x0=x0, CI_average=CI_average, lambda=lambda))
}
