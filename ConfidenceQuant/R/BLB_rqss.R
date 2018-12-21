#' Bag of Little Bootstraps.
#'
#' \code{BLB_rqss}  Bag of Litle Bootstraps (BLB) for data sets with 1-dim
#' covariate. Used to generate  confidence bands for a quantile smoothing
#' splines fitted with the rqss function from package \code{\link{quantreg}}.
#'
#' This function is based on the "Bag of Little Bootstraps" (BLB) method by
#' Kleiner et al.(2012). It performs parallelization to speed up the
#' calculation.
#'
#' @import quantreg
#' @import doParallel
#' @import foreach
#' @importFrom ggplot2 ggplot geom_point geom_ribbon geom_line labs theme aes element_text element_blank
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
#' @param D A number that specifies for how many x values you want to compute
#'   confidence intervals; the confidence band is made of pointwise intervals at
#'   varius x values. If specified, it will look at \eqn{D} equidistant points.
#' @param b The subsample size in the BLB algorithm. Kleiner et al. suggest that
#'   the size should be around \eqn{n^0.6}, where \eqn{n} is the data size.
#' @param s The number of subsamples used in BLB algorithm. Jordan et al.
#'   suggest that \eqn{s} should be 10~20.
#' @param r The number of bootstrap iterations (samples with with replacement).
#'   \eqn{r=100} is suggested.
#' @param Search If \code{TRUE} (which is recommended), then the function will
#'   first search for an optimum smoothing parameter \eqn{\lambda}.
#'@param xlab,ylab Titles for x and y axis. See \code{\link{title}}.
#' @return A list with three parts:
#'
#'   1. \code{x0} and \code{CI_average}, where \code{x0} contains the x values
#'   at which the confidence intervals are evaluated, and \code{CI_average} is
#'   2-dim matrix which contains the corresponding lower and upper bounds.
#'
#'   2. \code{lambda}, which is the optimum smoothing parameter selected by
#'   \code{BLB_rqss}. If it is done automatically, the function also returns
#'   \code{Lambda} and \code{Fid}, which respectively stand for a vector of lambda
#'   values and their corresponding cross-validation MCV values.
#'
#' @examples
#' data(one)
#' result<-BLB_rqss(cores=7, data=one, alpha=0.05, tau=0.5, Search=TRUE)
#'
#' plot<-FitPlot_1dim(data=one, result=result, xlab='x',ylab='y')
#' plot
#'
#' @references Micheal, I. J et al.(2012). \eqn{A Scalable Bootstrap for Massive
#'   Data}.
#' @seealso \code{\link{BLB_rqss_2dim}} for BLB with 2-dim covariate dataset.
#'
#'@export


BLB_rqss<-function(cores=NULL,data,alpha=0.05,tau=0.25,lambda=2,D=100,b=ceiling(nrow(data)^0.6), s=15, r=100, Search=FALSE){
  set.seed(100)
  if(is.null(cores)) {cores=round(0.5*detectCores())}
  registerDoParallel(cores=cores)

  if(!all(colnames(data)==c("x","y"))) colnames(data)<-c("x","y")

  ###1. Sample s disjoint subsets with size b
  n<-nrow(data)
  shuffle<-sample(1:n)
  indices<-matrix(head(shuffle,b*s),nrow=s,byrow=TRUE)
  Rboot<-rmultinom(s,size=n,prob=rep(1/b,b))
  indicator<-(Search|length(lambda)>1)


  ###2. λ Search
  if(indicator){
    if(length(lambda)>1) Lambda<-lambda else Lambda<-c(seq(0,400,by=10),500,1000,1500,2000,3000)
    fid<-rep(NA,length(Lambda))

    Fid<-foreach(j=1:s,.combine = "cbind") %dopar% {

      subsample<-data[indices[j,],]

      for(f in 1:length(Lambda)){
        mod <- rqss_new(y ~ qss(x, lambda = Lambda[f]), tau = tau, data = subsample, weights = as.vector(Rboot[,j]))
        fid_rec<-0
        for(p in (1:s)[-j]){
          testdata<-data[indices[p,],]
          extrapolated<-(testdata$x > max(subsample$x) | testdata$x < min(subsample$x))
          y.pred<-predict(mod, newdata = testdata[!extrapolated,])
          temp<-as.vector(Rboot[,p])[!extrapolated]
          fid_rec<-fid_rec+sum(temp*rho(tau=tau, testdata$y[!extrapolated]-y.pred))/sum(temp)
        }
        fid[f]<-fid_rec/(s-1)
      }

      return(fid)
    }

    #plot(Lambda[-1],rowMeans(Fid)[-1],type='l',ylab="New Metric",xlab="λ",main="MCV for BLB repeated data")
    Fidelity<-rowMeans(Fid)
    lambda<-Lambda[which.min(Fidelity)]
  }


  ###3. Bag of Little Bootstrap
  xrange<-range(data$x)
  x0 <- seq(xrange[1], xrange[2], by = diff(xrange)/D)
  response <- rep(NA,length(x0))   #Prediction for one bootstrap sample
  Response<-matrix(nrow=length(x0),ncol=r)

  CIs<-foreach (j=1:s) %dopar% {
    subsample<-data[indices[j,],]
    extrapolated = x0 > max(subsample$x) | x0 < min(subsample$x)
    for(k in 1:r){
      response[1:length(x0)] <- NA
      rboot<-as.vector(rmultinom(1,size=n,prob=rep(1/b,b)))
      mod <- rqss_new(y ~ qss(x, lambda = lambda), tau = tau, data = subsample, weights = rboot)
      response[!extrapolated] = predict(mod, newdata = data.frame(x=x0[!extrapolated]))
      Response[,k]<-response
    }
    return(t(apply(Response,1,function(x) quantile(x, probs=c(alpha/2,1-alpha/2),na.rm=TRUE))))
    #t(apply(Response,1,function(x) quantile(x, probs=c(alpha/(2*p),1-alpha/(2*p)),na.rm=TRUE)))
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
  #R CMD Rd2pdf ConfidenceQuant
}



#'\code{FitPlot_1dim} is a function to plot the results.
#' @name FitPlot_1dim
#' @rdname BLB_rqss
#' @export

FitPlot_1dim<-function(data,result,xlab='x',ylab='y'){
  if(!all(colnames(data)==c("x","y"))) colnames(data)<-c("x","y")

  n <- nrow(data)
  yo <- range(c(quantile(data$y, 0.001), quantile(data$y,0.999)))
  xo <- range(result$x0[!is.na(result$CI_average[,1])])


  dat<-data.frame(x0=result$x0,lower=result$CI_average[,1],upper=result$CI_average[,2])
  lower <- approxfun(result$x0, result$CI_average[, 1])
  upper<- approxfun(result$x0, result$CI_average[, 2])


  h <- ggplot(data[sample(n, max(min(n * 0.001, 5000),2999)), ])
  h <- h + geom_point(aes(x = x, y = y), alpha = 0.1, colour = "grey") +
    geom_ribbon(data = dat, aes(x = x0, ymin = lower, ymax = upper), alpha = 0.7, fill = "dodgerblue3") +
    geom_line(data = dat, aes(x = x0, y = upper), alpha = 0.7, colour = "dodgerblue3")+ coord_cartesian(xlim = xo, ylim = yo, expand = FALSE)  +
    labs(title = "Confidence Bands", x = xlab, y = ylab) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.title = element_blank(), legend.position = c(0.99,  0.99), legend.justification = c("right", "top"))
  return(h)
}
