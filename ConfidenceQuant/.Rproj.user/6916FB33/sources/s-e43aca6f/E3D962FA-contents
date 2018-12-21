#'Comparative Experiment.
#'
#'\code{CompExperiment_2dim_boot} compares the confidence regions for a given
#'quantile for two different datasets, one related to a treatment and the other
#'to a control. It applies the regular bootstrap algorithm to each dataset to get confidence
#'regions using quantile smoothing splines for 2-dim covariate.
#'
#'This function runs \code{Bootstrap} twice, once for each dataset. It is based
#'on \code{\link{Bootstrap_rqss_2dim}}, which implements regular bootstrap for
#'quantile smoothing splines with two dimensional covariate dataset. It
#'performs parallelization to speed up the calculation.
#'
#'\if{html}{\figure{comp2.png}{options: width=100 alt="Image output"}}
#'\if{latex}{\figure{comp2.png}{options: width=3in}}
#'
#'This function can also return a colored contour plot, in which different
#'colors represent different scenarios. See figure above.
#'
#'@import quantreg
#'@import doParallel
#'@import foreach
#'@import akima
#'
#'@param cores The number of cores to use for parallel execution. If not
#'  specified, the number of cores is set to the value of
#'  \code{options("cores")}, if specified, or to one-half the number of cores
#'  detected by the \code{\link{parallel}} package.
#'@param treatment A 3-dim optional data frame for the treatment group (or
#'  object coercible by \code{\link{as.data.frame}} to a data frame) containing
#'  the variables in the model. The last column of this data frame must be the
#'  response for the experiment.
#'@param control A 3-dim optional data frame for the control group (or object
#'  coercible by \code{\link{as.data.frame}} to a data frame) containing the
#'  variables in the model. The last column of this data frame must be the
#'  response for the experiment.
#'@param alpha The confidence level required. The default is 0.05.
#'@param tau A specific quantile to be estimated. Must be a number between 0 and
#'  1.
#'@param lambda The smoothing parameter used for \code{treatment} &
#'  \code{control} if \eqn{Search=FALSE}, which governs the tradeoff between
#'  fidelity and the penalty component for the triogram term.
#'@param D A number that determines the density of your grid of x values that
#'  you want to predict on. If specified, we will examine the confidence bands
#'  on a \eqn{D×D} grid.
#' @param B The number of Monte Carlo iterations using bootstrap with
#'   replacement. \eqn{B=100} is by default.
#' @param Search If \code{TRUE} (which is recommended), then the function will
#'   first search for an optimum smoothing parameter \eqn{\lambda}.
#'
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
#'  \code{BLB_rqss_2dim}. If it is done automatically, the function also returns
#'  \code{Lambda} and \code{Fid}, which respectively stand for a vector of lambda
#'  values and their corresponding cross-validation MCV values.

#' @examples
#' data(treatment_2dim)
#' data(control_2dim)
#'
#' #alpha=0.05;tau=0.5
#' all<-CompExperiment_2dim_boot(cores=7, treatment_2dim, control_2dim, tau=0.5, B=100, Search=TRUE)
#'
#' plot<-CompPlot_2dim(control = control_2dim,treatment = treatment_2dim,all = all,xlab="x1",
#' ylab="x2",zlab="z")
#'
#'
#'@references Micheal, I. J et al.(2012). \eqn{A Scalable Bootstrap for Massive
#'  Data}.
#'@references Akima, H. (1978). \eqn{A Method of Bivariate Interpolation and
#'  Smooth Surface Fitting for Irregularly Distributed Data Points}. ACM
#'  Transactions on Mathematical Software 4, 148-164.
#'@seealso \code{\link{contour},\link{image}}
#'@seealso \code{\link{Bootstrap_rqss_2dim}} for simple bootstrap method with
#'  one dataset that has 1-dim covariate.
#'@seealso \code{\link{CompExperiment_1dim_boot}} for comparative experiments
#'  with 1-dim covariate data sets.
#'
#'@export



CompExperiment_2dim_boot<-function(cores=NULL,treatment,control,alpha=0.05,tau=0.25,lambda=2,D=50, B=100, Search=FALSE){
  if(is.null(cores)) {cores=round(0.5*detectCores())}
  registerDoParallel(cores=cores)

  if(!all(colnames(treatment)==c("x1","x2","y"))) colnames(treatment)<-c("x1","x2","y")
  if(!all(colnames(control)==c("x1","x2","y"))) colnames(control)<-c("x1","x2","y")

  indicator<-(Search|length(lambda)>1)

  ### λ Search
  r1<-list(lambda=lambda)
  r2<-list(lambda=lambda)

  if(indicator){
    if(length(lambda)>1) Lambda<-lambda else Lambda<-c(seq(0,150,by=15),seq(170,400,by=20),400,500,1000,1500,2000,3000)
    folds1<-foldsGenerator(sim.data=treatment,nfolds=10)
    folds2<-foldsGenerator(sim.data=control,nfolds=10)


    Fid1<-foreach(j=1:length(Lambda),.combine = "c") %dopar% {
      multifoldCV_2dim(cores=cores,lambda=Lambda[j],tau,sim.data=treatment,cv.idx=folds1)
    }

    Fid2<-foreach(j=1:length(Lambda),.combine = "c") %dopar% {
      multifoldCV_2dim(cores=cores,lambda=Lambda[j],tau,sim.data=control,cv.idx=folds2)
    }

    r1<-list(lambda=Lambda[which.min(Fid1)],Lambda=Lambda,Fidelity=Fid1)
    r2<-list(lambda=Lambda[which.min(Fid2)],Lambda=Lambda,Fidelity=Fid2)
  }

  x1range<-sort(c(range(treatment$x1), range(control$x1)))[2:3]
  x2range<-sort(c(range(treatment$x2), range(control$x2)))[2:3]
  x0 <- expand.grid(x1=seq(x1range[1], x1range[2], by = diff(x1range)/D),x2=seq(x2range[1], x2range[2], by = diff(x2range)/D))


  temp<-Bootstrap_rqss_2dim(cores=cores,data=treatment,alpha=alpha,tau=tau,lambda=r1$lambda,B=B,D=D,x0=x0,Search=FALSE)
  if(indicator) result1<-list(x0=temp$x0, CI_average=temp$CI_average, lambda=r1$lambda, Lambda=r1$Lambda, Fid=r1$Fidelity)
  else result1<-list(x0=temp$x0, CI_average=temp$CI_average, lambda=r1$lambda)

  temp<-Bootstrap_rqss_2dim(cores=cores,data=control,alpha=alpha,tau=tau,lambda=r2$lambda,B=B,D=D,x0=x0,Search=FALSE)
  if(indicator) result2<-list(x0=temp$x0, CI_average=temp$CI_average, lambda=r2$lambda, Lambda=r2$Lambda, Fid=r2$Fidelity)
  else result2<-list(x0=temp$x0, CI_average=temp$CI_average, lambda=r2$lambda)

  Diff<-list(x0=x0,CI_average=result1$CI_average-result2$CI_average)

  return(list(result1=result1, result2=result2, Diff=Diff))
}
