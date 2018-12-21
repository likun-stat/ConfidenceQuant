#'Bag of Little Bootstraps.
#'
#'\code{BLB_rqss_2dim}  Bag of Litle Bootstraps (BLB) for data sets with 2-dim
#'covariate. Used to generate confidence bands for a quantile smoothing splines
#'fitted with the rqss function from package \code{\link{quantreg}}.
#'
#'This function is based on the "Bag of Little Bootstraps" (BLB) method by
#'Kleiner et al.(2012). It performs parallelization to speed up the calculation.
#'
#'@import quantreg
#'@import doParallel
#'@import foreach
#'@import tripack
#'@import akima
#'@import plotly
#'
#'@param cores The number of cores to use for parallel execution. If not
#'  specified, the number of cores is set to the value of
#'  \code{options("cores")}, if specified, or to one-half the number of cores
#'  detected by the \code{\link{parallel}} package.
#'@param data A 3-dim optional data frame (or object coercible by
#'  \code{\link{as.data.frame}} to a data frame) containing the variables in the
#'  model. The column names should be specified in a way that “x1” and "x2" are
#'  for the two predictors, and “y” for the response.
#'@param alpha The confidence level required. The default is 0.05.
#'@param tau A specific quantile to be estimated. Must be a number between 0 and
#'  1.
#'@param D A number that specifies for how many x values you want to compute
#'  confidence intervals; the confidence surface is made of pointwise intervals
#'  at varius x values. If specified, it will compute intervals on an
#'  equidistant \eqn{D×D} grid.
#'@param lambda The smoothing parameter governing the tradeoff between fidelity
#'  and the penalty component for the triogram term. If \code{Search=TRUE},
#'  there is no need for users to specify a value.
#'@param b The subsample size in the BLB algorithm. Kleiner et al. suggest that
#'  the size should be around \eqn{n^0.6}, where \eqn{n} is the data size.
#'@param s The number of subsamples used in BLB algorithm. Kleiner et al.
#'  suggest that \eqn{s} should be 10~20.
#'@param r The number of bootstrap iterations (samples with with replacement).
#'  \eqn{r=100} is suggested.
#'@param Search If \code{TRUE} (which is recommended), then the function will
#'  first search for an optimum smoothing parameter \eqn{\lambda}.
#'@param xlab,ylab,zlab Titles for the axises. See \code{\link{title}}.
#'
#'@return A list with three parts:
#'
#'  1. \code{x0} and \code{CI_average}, where \code{x0} contains the x values at
#'  which the confidence intervals are evaluated, and \code{CI_average} is 2-dim
#'  matrix which contains the corresponding lower and upper bounds.
#'
#'  2. \code{lambda}, which is the optimum smoothing parameter selected by
#'  \code{BLB_rqss_2dim}. If it is done automatically, the function also returns
#'  \code{Lambda} and \code{Fid}, which respectively stand for a vector of lambda
#'  values and their corresponding cross-validation MCV values.
#'
#'
#' @examples
#' data(two)
#' result<-BLB_rqss_2dim(cores=7, data=two, alpha=0.05, tau=0.5, Search=TRUE)
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
#'@references Kleiner, I. J et al. JRSS B, 2012. \eqn{A Scalable Bootstrap for
#'  Massive Data}.
#'@references Roger Koenker.(2005). \eqn{Quantile regression}. Cambridge
#'  university press.
#'@seealso \code{\link{BLB_rqss}} for BLB with 1-dim covariate dataset.
#'
#'@export






BLB_rqss_2dim<-function(cores=NULL,data,alpha=0.05,tau=0.25,lambda=2,D=50,b=ceiling(nrow(data)^0.6), s=15, r=100, Search=FALSE){
  if(is.null(cores)) {cores=round(0.5*detectCores())}
  registerDoParallel(cores=cores)

  if(!all(colnames(data)==c("x1","x2","y"))) colnames(data)<-c("x1","x2","y")

  ###1. Sample s disjoint subsets with size b
  n<-nrow(data)
  shuffle<-sample(1:n)
  indices<-matrix(head(shuffle,b*s),nrow=s,byrow=TRUE)
  indicator<-(Search|length(lambda)>1)

  ###2. λ Search
  if(indicator){
    if(length(lambda)>1) Lambda<-lambda else Lambda<-c(seq(0,150,by=15),seq(170,400,by=20),400,500,1000,1500,2000,3000)
    fid<-rep(NA,length(Lambda))

    Fid<-foreach(j=1:s,.combine = "cbind") %dopar% {
      rboot1<-as.vector(rmultinom(1,size=n,prob=rep(1/b,b)))
      rboot2<-as.vector(rmultinom(1,size=n,prob=rep(1/b,b)))
      subsample<-data[indices[j,],]
      if(j==s) testdata<-data[indices[1,],] else testdata<-data[indices[j+1,],]

      mod <- rqss(y ~ qss(cbind(x1,x2), lambda = lambda[1], ndum=0), tau = tau, data = subsample)
      tri<-tripack::tri.mesh(mod$qss[[1]][[1]][,1], mod$qss[[1]][[1]][,2])
      NONextrapolated<-tripack::in.convex.hull(tri, testdata$x1, testdata$x2)


      for(f in 1:length(Lambda)){
        mod <- rqss_new(y ~ qss(cbind(x1,x2), lambda = Lambda[f], ndum=0), tau = tau, data = subsample, weights = rboot1)
        y.pred<-predict(mod, newdata = testdata[NONextrapolated,])
        temp<-rboot2[NONextrapolated]
        fid[f]<-sum(temp*rho(tau=tau, testdata$y[NONextrapolated]-y.pred))/sum(temp)
      }
      return(fid)
    }

    #plot(Lambda[-1],rowMeans(Fid)[-1],type='l',ylab="New Metric",xlab="λ",main="MCV for BLB repeated data")
    Fidelity<-rowMeans(Fid)
    lambda<-Lambda[which.min(Fidelity)]
  }


  ###3. Bag of Little Bootstrap
  x1range<-range(data$x1); x2range<-range(data$x2)
  x0 <- expand.grid(x1=seq(x1range[1], x1range[2], by = diff(x1range)/D),x2=seq(x2range[1], x2range[2], by = diff(x2range)/D))   #Changed from 100 to 50*50
  response <- rep(NA,nrow(x0))   #Prediction for one bootstrap sample
  Response<-matrix(nrow=nrow(x0),ncol=r)

  CIs<-foreach (j=1:s) %dopar% {
    subsample<-data[indices[j,],]

    #x0 should not be outside of the convex hull
    mod <- rqss(y ~ qss(cbind(x1,x2), lambda = lambda[1], ndum=0), tau = tau, data = subsample)
    tri<-tripack::tri.mesh(mod$qss[[1]][[1]][,1], mod$qss[[1]][[1]][,2])
    NONextrapolated<-tripack::in.convex.hull(tri, x0$x1, x0$x2)


    for(k in 1:r){
      response[1:nrow(x0)] <- NA
      rboot<-as.vector(rmultinom(1,size=n,prob=rep(1/b,b)))
      mod <- rqss_new(y ~ qss(cbind(x1,x2), lambda = lambda, ndum=0), tau = tau, data = subsample, weights = rboot)
      response[NONextrapolated] = predict(mod, newdata = x0[NONextrapolated,])
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



#'\code{FitPlot_2dim} is a function to plot the results.
#' @name FitPlot_2dim
#' @rdname BLB_rqss_2dim
#' @export


FitPlot_2dim<-function(data, result, Range=c(0.002,0.998), xlab="log(avtp_med)", ylab="log(niqr_med)", zlab="log(playdelay_med)"){
  if(!all(colnames(data)==c("x1","x2","y"))) colnames(data)<-c("x1","x2","y")
  n<-nrow(data)
  x0<-result$x0
  CI_average<-result$CI_average
  X<-x0$x1;Y<-x0$x2;Z<-CI_average[,1];Z_1<-CI_average[,2]
  X<-X[!is.na(Z)];Y<-Y[!is.na(Z)];Z_1<-Z_1[!is.na(Z)];Z<-Z[!is.na(Z)]


  m<-300
  xo <- range(X)
  yo <- range(Y)
  xo1<-quantile(data$x1,probs=Range)
  yo1<-quantile(data$x2,probs=Range)
  zo<-quantile(data$y,probs=Range)
  xo <- sort(c(xo,xo1))[2:3]
  yo <- sort(c(yo,yo1))[2:3]

  akima_lower<-interp(X,Y,Z,xo=seq(xo[1],xo[2],length=m),yo=seq(yo[1],yo[2],length=m+1))
  akima_higher<-interp(X,Y,Z_1,xo=seq(xo[1],xo[2],length=m),yo=seq(yo[1],yo[2],length=m+1))

  f2 <- list(
    family = "Old Standard TT, serif",
    size = 14,
    color="grey")

  scene<-list(aspectmode="manual",aspectratio=list(x=1,y=1,z=0.95),
              xaxis=list(title=xlab,titlefont=f2,range=xo),
              yaxis=list(title=ylab,titlefont=f2,range=yo),
              zaxis=list(title=zlab,titlefont=f2,range=zo,showline=TRUE),
              camera = list(eye = list(x = 2.25, y = -2.25, z = 0.9)))

  a <- list(
    text = "Confidence Bands",
    xref = "paper",
    yref = "paper",
    yanchor = "bottom",
    xanchor = "center",
    align = "center",
    x = 0.5,
    y = 0.92,
    showarrow = FALSE
  )


  zo<-range(result$CI_average,na.rm=TRUE)
  zdiff<-diff(zo)
  col<-colorRampPalette(c("red","green"))(3)

  Tmp<-data[sample(n, max(min(n * 0.001, 5000),2999)), ]
  trace<-plot_ly(Tmp,x=~x1, y=~x2, z=~y, showlegend=FALSE)%>%
    add_markers(color=I("grey"), size=I(2.5), symbol=I(20), opacity=0.2,hoverinfo=FALSE)%>%
    add_surface(x=akima_lower$x,y=akima_lower$y,z=t(akima_lower$z),opacity=0.7, showscale=FALSE) %>%
    add_surface(x=akima_higher$x,y=akima_higher$y,z=t(akima_higher$z),opacity=0.7, showscale=FALSE)%>%
    layout(annotations = a, scene=scene)
  return(trace)
}
