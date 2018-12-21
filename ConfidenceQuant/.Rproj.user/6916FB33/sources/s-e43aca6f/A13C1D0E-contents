#' Regular Bootstrap Method.
#'
#' \code{Bootstrap_rqss}  Regular bootstrap for data sets with 1-dim covariate.
#' Used to generate the confidence bands for a quantile smoothing spline fitted
#' with the rqss function from package \code{\link{quantreg}}.
#'
#' This function is based on a regular bootstrap method, which calculates
#' confidence bands for one quantile. It performs parallelization to speed up
#' the calculation.
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
#' @param lambda The smoothing parameter governing the tradeoff between fidelity
#'   and the penalty component for the triogram term. If \code{Search=TRUE},
#'   there is no need for users to specify a value.
#' @param B The number of Monte Carlo iterations using bootstrap with
#'   replacement. \eqn{B=100} is by default.
#' @param D A number that specifies how many x values you want to examine the
#'   confidence bands at. If specified, it will look at \eqn{D} equal-distanced
#'   points.
#' @param Search If \code{TRUE} (which is recommended), then the function will
#'   first search for an optimum smoothing parameter \eqn{\lambda}.
#' @return A list with two parts: \code{x0} and \code{CIs}, where \code{x0}
#'   contains the x values that we are examining the confidence intervals at,
#'   and \code{CIs} is 2-dim matrix which contains the corresponding lower bound
#'   and upper bound.
#' @examples
#' data(one)
#' result<-Bootstrap_rqss(data=one, alpha=0.05, tau=0.5, B=100, Search=TRUE)
#'
#'
#' plot<-FitPlot_1dim(data=one, result=result, xlab='x',ylab='y')
#' plot
#'
#' @references Micheal, I. J et al.(2012). \eqn{A Scalable Bootstrap for Massive
#'   Data}.
#' @seealso \code{\link{BLB_rqss}} for BLB with 1-dim covariate data set.
#' @seealso \code{\link{Bootstrap_rqss_2dim}} for regular Bootstrap method with
#'   2-dim covariate data set.
#'
#'@export




Bootstrap_rqss<-function(cores=NULL,data=parent.frame(),alpha=0.05,tau=0.25,lambda=2,B=100,D=100, Search=FALSE, x0=NULL, warning.catch=FALSE){
  if(is.null(cores)) {cores=3}
  registerDoParallel(cores=cores)

  if(missing(data)) stop("Missing data input. Need to give values.")

  if(!all(colnames(data)==c("x","y"))) colnames(data)<-c("x","y")

  n<-nrow(data)
  subsample<-data.frame(matrix(nrow=n,ncol=ncol(data)))   #Store the subsample
  colnames(subsample)<-colnames(data)

  if(is.null(x0)) {xrange<-range(data$x)
  x0 <- seq(xrange[1], xrange[2], by = diff(xrange)/D)}
  response <- rep(NA,length(x0))   #Prediction for one bootstrap sample
  indicator<-(Search|length(lambda)>1)


  #λ Search
  if(indicator){
    if(length(lambda)>1) Lambda<-lambda else Lambda<-c(seq(0,400,by=10),500,1000,1500,2000,3000)
    folds<-foldsGenerator(sim.data=data,nfolds=10)

    Fid<-foreach(j=1:length(Lambda),.combine = "c") %dopar% {
      multifoldCV(cores=cores,lambda=Lambda[j],tau,sim.data=data,cv.idx=folds)
    }

    Fidelity<-Fid
    lambda<-Lambda[which.min(Fidelity)]
  }

  #Bootstrap
  Response<-foreach (j=1:B,.combine='rbind') %dopar% {
    #Subsample the data
    ind<-FALSE
    subindex<-sample(1:n,n,replace=TRUE)
    subsample[1:n,]<-data[subindex,]
    extrapolated = x0 > max(subsample$x) | x0 < min(subsample$x)

    response[1:length(x0)] <- NA
    if(warning.catch){
      tt <- tryCatch(rqss(y ~ qss(x, lambda = lambda), tau = tau, data = subsample) ,error=function(e) e, warning=function(w) w)
      if(is(tt,"warning")) ind<-TRUE
      if(!ind){
        mod<-rqss(y ~ qss(x, lambda = lambda), tau = tau, data = subsample)
        response[!extrapolated] = predict(mod, newdata = data.frame(x=x0[!extrapolated]))
      }}
    else{
      mod<-rqss(y ~ qss(x, lambda = lambda), tau = tau, data = subsample)
      response[!extrapolated] = predict(mod, newdata = data.frame(x=x0[!extrapolated]))
    }
    c(response,ind)
  }

  #p<-length(x0)
  p<-1
  if(warning.catch){
    filter<-sum(Response[,102])
    cat(paste("WARNING CATCH: ",filter," iterations were disregarded out of ", B,"\n"))
    CIs<-t(apply(Response[!Response[,102],1:101],2,function(x) quantile(x, probs=c(alpha/(2*p),1-alpha/(2*p)),na.rm=TRUE)))}
  else {CIs<-t(apply(Response[,1:101],2,function(x) quantile(x, probs=c(alpha/(2*p),1-alpha/(2*p)),na.rm=TRUE)))}

  if(indicator) return(list(x0=x0, CI_average=CIs, lambda=lambda, Lambda=Lambda, Fid=Fidelity))
  else return(list(x0=x0, CI_average=CIs, lambda=lambda))
}
