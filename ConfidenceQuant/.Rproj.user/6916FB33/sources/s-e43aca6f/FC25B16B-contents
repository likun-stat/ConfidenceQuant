#' Regular Bootstrap Method.
#'
#' \code{Bootstrap_rqss_2dim}  Regualr bootstrap for data sets with 2-dim
#' covariate. Used to generate the confidence bands for a quantile smoothing
#' splines fitted with the rqss function from package \code{\link{quantreg}}.
#'
#' This function is based on a regular bootstrap method, which calculates
#' confidence bands for one quantile. It performs parallelization to speed up
#' the calculation.
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
#'   are for the two predictors respectively, and “y” for the response.
#' @param alpha The confidence level required. The default is 0.05.
#' @param tau A specific quantile to be estimated. Must be a number between 0
#'   and 1.
#' @param D A number that determines the density of a grid of x values over
#' which pointwise confidence intervals will be computed. If specified,  a
#' \eqn{D×D} grid of equidistant points on the plane is created.
#' @param lambda The smoothing parameter governing the tradeoff between fidelity
#'   and the penalty component for the triogram term. If \code{Search=TRUE},
#'   there is no need for users to specify a value.
#' @param B The number of Monte Carlo iterations using bootstrap with
#'   replacement. \eqn{B=100} is by default.
#' @param Search If \code{TRUE} (which is recommended), then the function will
#'   first search for an optimum smoothing parameter \eqn{\lambda}.
#' @return A list with two parts: \code{x0} and \code{CIs}, where \code{x0}
#'   contains the x values that we are examining the confidence intervals at,
#'   and \code{CIs} is 2-dim matrix which contains the corresponding lower bound
#'   and upper bound.
#' @examples
#' data(two)
#' result<-Bootstrap_rqss_2dim(data=two, alpha=0.05, tau=0.5,B=100,Search=TRUE)
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
#' @references Micheal, I. J et al.(2012). \eqn{A Scalable Bootstrap for Massive
#'   Data}.
#' @references Roger Koenker.(2005). \eqn{Quantile regression}. Cambridge
#'   university press.
#' @seealso \code{\link{BLB_rqss_2dim}} for BLB with 2-dim covariate data set.
#' @seealso \code{\link{Bootstrap_rqss}} for regular Bootstrap method with
#'   1-dim covariate data set
#'
#'@export






Bootstrap_rqss_2dim<-function(cores=NULL,data=parent.frame(),alpha=0.05,tau=0.25,lambda=2,B=100,D=50,x0=NULL,Search=TRUE){
  if(is.null(cores)) {cores=3}
  registerDoParallel(cores=cores)

  if(missing(data)) stop("Missing data input. Need to give values.")

  if(!all(colnames(data)==c("x1","x2","y"))) colnames(data)<-c("x1","x2","y")

  n<-nrow(data)
  subsample<-data.frame(matrix(nrow=n,ncol=ncol(data)))   #Store the subsample
  colnames(subsample)<-colnames(data)

  x1range<-range(data$x1); x2range<-range(data$x2)
  if(is.null(x0)) x0 <- expand.grid(x1=seq(x1range[1], x1range[2], by = diff(x1range)/D),x2=seq(x2range[1], x2range[2], by = diff(x2range)/D))   #Changed from 100 to 50*50
  response <- rep(NA,nrow(x0))   #Prediction for one bootstrap sample
  #Response<-matrix(nrow=nrow(x0),ncol=r)
  indicator<-(Search|length(lambda)>1)


  #λ Search
  if(indicator){
    if(length(lambda)>1) Lambda<-lambda else  Lambda<-c(seq(0,150,by=15),seq(170,400,by=20),400,500,1000,1500,2000,3000)
    folds<-foldsGenerator(sim.data=data,nfolds=10)

    Fid<-foreach(j=1:length(Lambda),.combine = "c") %dopar% {
      multifoldCV_2dim(cores=cores,lambda=Lambda[j],tau,sim.data=data,cv.idx=folds)
    }

    Fidelity<-Fid
    lambda<-Lambda[which.min(Fidelity)]
  }

  #Bag of Little Bootstrap
  Response<-foreach (j=1:B,.combine='rbind') %dopar% {
    #Subsample the data
    ind<-FALSE
    subindex<-sample(1:n,n,replace=TRUE)
    subsample[1:n,]<-data[subindex,]

    #x0 should not be outside of the convex hull
    mod <- rqss(y ~ qss(cbind(x1,x2), lambda = lambda, ndum=0), tau = tau, data = subsample)
    tri<-tripack::tri.mesh(mod$qss[[1]][[1]][,1], mod$qss[[1]][[1]][,2])
    NONextrapolated<-tripack::in.convex.hull(tri, x0$x1, x0$x2)
    response[NONextrapolated] = predict(mod, newdata = x0[NONextrapolated,])
    response
  }

  #p<-nrow(x0)
  p<-1
  CIs<-t(apply(Response,2,function(x) quantile(x, probs=c(alpha/(2*p),1-alpha/(2*p)),na.rm=TRUE)))

  if(indicator) return(list(x0=x0, CI_average=CIs, lambda=lambda, Lambda=Lambda, Fid=Fidelity))
  else return(list(x0=x0, CI_average=CIs, lambda=lambda))
}
