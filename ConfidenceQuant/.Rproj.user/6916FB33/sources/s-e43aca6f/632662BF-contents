#'Comparative Experiment.
#'
#'\code{CompExperiment_2dim} compares the confidence regions for a given
#'quantile for two different datasets, one related to a treatment and the other
#'to a control. It applies the BLB algorithm to each dataset to get confidence
#'regions using quantile smoothing splines for 2-dim covariate.
#'
#'This function runs \code{BLB} twice, once for each dataset. It is based on
#'\code{\link{BLB_rqss_2dim}}, which implements BLB for quantile smoothing splines
#'with a two-dimensional covariate dataset. It performs parallelization to speed
#'up the calculation.
#'
#'\if{html}{\figure{comp2.png}{options: width=100 alt="Image output"}}
#'\if{latex}{\figure{comp2.png}{options: width=3in}}
#'
#'\code{CompPlot_2dim} takes the results and use ggplot/plotly to visualize
#'them, in which different colors represent different scenarios. See figure
#'above.
#'
#'@import quantreg
#'@import doParallel
#'@import foreach
#'@import akima
#'@import plotly
#'@import grDevices
#'@importFrom ggplot2 ggplot geom_point geom_ribbon geom_line labs theme geom_hline geom_rect coord_cartesian scale_fill_manual scale_fill_gradientn
#'       geom_raster geom_tile geom_contour aes element_rect element_line
#'
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
#'@param D A number that determines the density of a grid of x values at which
#'  the quantile function will be predicted. If specified, it will evaluate a
#'  confidence surface on a \eqn{D×D} grid.
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
#'@param xlab,ylab,zlab Titles for the axises. See \code{\link{title}}.
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
#'
#' @examples
#' data(treatment_2dim)
#' data(control_2dim)
#'
#' #alpha=0.05;tau=0.5
#' all<-CompExperiment_2dim(cores=7, treatment_2dim, control_2dim, tau=0.5, Search=TRUE)
#'
#' plot<-CompPlot_2dim(control = control_2dim,treatment = treatment_2dim,all = all,xlab="x1",
#' ylab="x2",zlab="z")
#'
#'
#' @references Kleiner, I. J et al. JRSS B, 2012. \eqn{A Scalable Bootstrap for
#' Massive Data}.
#'@references Akima, H. (1978). \eqn{A Method of Bivariate Interpolation and Smooth Surface Fitting for Irregularly Distributed Data Points}. ACM Transactions on Mathematical Software 4, 148-164.
#'@seealso \code{\link{contour},\link{image}}
#'@seealso \code{\link{BLB_rqss_2dim}} for BLB with one dataset that has 1-dim
#'  covariate.
#'@seealso \code{\link{CompExperiment_1dim}} for comparative experiments with
#'  1-dim covariate data sets.
#'
#'
#'@export
#'
CompExperiment_2dim<-function(cores=NULL,treatment,control,alpha=0.05,tau=0.25,lambda=2,D=50,b1=ceiling(nrow(treatment)^0.6),b2=ceiling(nrow(control)^0.6), s=15, r=100, Search=FALSE){
  if(is.null(cores)) {cores=round(0.5*detectCores())}
  registerDoParallel(cores=cores)

  if(!all(colnames(treatment)==c("x1","x2","y"))) colnames(treatment)<-c("x1","x2","y")
  if(!all(colnames(control)==c("x1","x2","y"))) colnames(control)<-c("x1","x2","y")

  indicator<-(Search|length(lambda)>1)

  ###1. Sample s disjoint subsets with size b
  n<-nrow(treatment)
  shuffle<-sample(1:n)
  indices1<-matrix(head(shuffle,b1*s),nrow=s,byrow=TRUE)

  m<-nrow(control)
  shuffle<-sample(1:m)
  indices2<-matrix(head(shuffle,b2*s),nrow=s,byrow=TRUE)


  ###2. λ Search
  r1<-list(lambda=lambda)
  r2<-list(lambda=lambda)

  if(indicator){
    r1<-searchLambda_2dim(lambda=lambda,data=treatment,s=s,indices=indices1,tau=tau,n=n,b=b1)
    r2<-searchLambda_2dim(lambda=lambda,data=control,s=s,indices=indices2,tau=tau,n=m,b=b2)
    #plot(r$Lambda[-1],r$Fidelity[-1],type='l',ylab="New Metric",xlab="λ",main="MCV for BLB repeated data")
  }


  ###3. Bag of Little Bootstrap
  x1range<-sort(c(range(treatment$x1), range(control$x1)))[2:3]
  x2range<-sort(c(range(treatment$x2), range(control$x2)))[2:3]
  x0 <- expand.grid(x1=seq(x1range[1], x1range[2], by = diff(x1range)/D),x2=seq(x2range[1], x2range[2], by = diff(x2range)/D))

  response_trt <- rep(NA,nrow(x0))   #Prediction for one bootstrap sample_treatment
  Response_trt<-matrix(nrow=nrow(x0),ncol=r)

  response_ctr <- rep(NA,nrow(x0))   #Prediction for one bootstrap sample_control
  Response_ctr<-matrix(nrow=nrow(x0),ncol=r)

  Difference<-matrix(nrow=nrow(x0),ncol=r)

  CIs<-foreach (j=1:s, .inorder = FALSE) %dopar% {
    subsample1<-treatment[indices1[j,],]
    subsample2<-control[indices2[j,],]

    #x0 should not be outside of the convex hull
    mod1 <- rqss(y ~ qss(cbind(x1,x2), lambda = lambda[1], ndum=0), tau = tau, data = subsample1)
    tri1<-tripack::tri.mesh(mod1$qss[[1]][[1]][,1], mod1$qss[[1]][[1]][,2])
    NONextrapolated1<-tripack::in.convex.hull(tri1, x0$x1, x0$x2)

    mod2 <- rqss(y ~ qss(cbind(x1,x2), lambda = lambda[1], ndum=0), tau = tau, data = subsample2)
    tri2<-tripack::tri.mesh(mod2$qss[[1]][[1]][,1], mod2$qss[[1]][[1]][,2])
    NONextrapolated2<-tripack::in.convex.hull(tri2, x0$x1, x0$x2)

    for(k in 1:r){
      response_trt[1:nrow(x0)] <- NA
      response_ctr[1:nrow(x0)] <- NA
      rboot1<-as.vector(rmultinom(1,size=n,prob=rep(1/b1,b1)))
      rboot2<-as.vector(rmultinom(1,size=m,prob=rep(1/b2,b2)))

      mod1 <- rqss_new(y ~ qss(cbind(x1,x2), lambda = r1$lambda, ndum=0), tau = tau, data = subsample1, weights = rboot1)
      response_trt[NONextrapolated1] = predict(mod1, newdata = x0[NONextrapolated1,])
      Response_trt[,k]<-response_trt

      mod2 <- rqss_new(y ~ qss(cbind(x1,x2), lambda = r2$lambda, ndum=0), tau = tau, data = subsample2, weights = rboot2)
      response_ctr[NONextrapolated2] = predict(mod2, newdata = x0[NONextrapolated2,])
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
  CI_average<-matrix(NA,nrow = nrow(x0), ncol = 2)
  for(i in 1:nrow(x0)){
    times<-0
    ci<-rep(0,2)
    for (j in 1:length(CIs)){
      if(!is.na(CIs[[j]][[1]][i,1]))  {times<-times+1;ci<-ci+CIs[[j]][[1]][i,]}
    }
    if(times>0) CI_average[i,]<-ci/times
    else CI_average[i,]<-rep(NA,2)
  }

  if(indicator) result1<-list(x0=x0, CI_average=CI_average, lambda=r1$lambda, Lambda=r1$Lambda, Fid=r1$Fidelity)
  else result1<-list(x0=x0, CI_average=CI_average, lambda=lambda)

  #2. control
  CI_average<-matrix(NA,nrow = nrow(x0), ncol = 2)
  for(i in 1:nrow(x0)){
    times<-0
    ci<-rep(0,2)
    for (j in 1:length(CIs)){
      if(!is.na(CIs[[j]][[2]][i,1]))  {times<-times+1;ci<-ci+CIs[[j]][[2]][i,]}
    }
    if(times>0) CI_average[i,]<-ci/times
    else CI_average[i,]<-rep(NA,2)
  }

  if(indicator) result2<-list(x0=x0, CI_average=CI_average, lambda=r2$lambda, Lambda=r2$Lambda, Fid=r2$Fidelity)
  else result2<-list(x0=x0, CI_average=CI_average, lambda=lambda)

  #3. difference
  CI_average<-matrix(NA,nrow = nrow(x0), ncol = 2)
  for(i in 1:nrow(x0)){
    times<-0
    ci<-rep(0,2)
    for (j in 1:length(CIs)){
      if(!is.na(CIs[[j]][[3]][i,1]))  {times<-times+1;ci<-ci+CIs[[j]][[3]][i,]}
    }
    if(times>0) CI_average[i,]<-ci/times
    else CI_average[i,]<-rep(NA,2)
  }

  Diff<-list(x0=x0, CI_average=CI_average)

  return(list(result1=result1, result2=result2, Diff=Diff))
}


#'\code{searchLambda_2dim} is a wrapper function to calculate the optimum lambda.
#' @name searchLambda_2dim
#' @rdname CompExperiment_2dim
#' @export

searchLambda_2dim<-function(lambda,data,s,indices,tau,n,b){
  if(length(lambda)>1) Lambda<-lambda else Lambda<-c(seq(0,150,by=15),seq(170,400,by=20),400,500,1000,1500,2000,3000)
  fid<-rep(NA,length(Lambda))

  Fid<-foreach(j=1:s,.combine = "cbind", .inorder = FALSE) %dopar% {
    rboot1<-as.vector(rmultinom(1,size=n,prob=rep(1/b,b)))
    rboot2<-as.vector(rmultinom(1,size=n,prob=rep(1/b,b)))
    subsample<-data[indices[j,],]
    if(j==s) testdata<-data[indices[1,],] else testdata<-data[indices[j+1,],]

    mod <- rqss(y ~ qss(cbind(x1,x2), lambda = lambda, ndum=0), tau = tau, data = subsample)
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

  return(list(lambda=lambda,Lambda=Lambda,Fidelity=Fidelity))
}


#'\code{CompPlot_2dim} is a function to plot the results.
#' @name CompPlot_2dim
#' @rdname CompExperiment_2dim
#' @export

CompPlot_2dim<-function(control, treatment, all, Range=c(0.002,0.998) ,xlab="log(avtp_med)",ylab="log(niqr_med)",zlab="log(playdelay_med)"){
  if(!all(colnames(treatment)==c("x1","x2","y"))) colnames(treatment)<-c("x1","x2","y")
  if(!all(colnames(control)==c("x1","x2","y"))) colnames(control)<-c("x1","x2","y")
  result1<-all$result1
  result2<-all$result2
  x0_1<-result1$x0
  CI_average_1<-result1$CI_average
  X1<-x0_1$x1;Y1<-x0_1$x2;Z1<-CI_average_1[,1];Z1_1<-CI_average_1[,2]
  X1<-X1[!is.na(Z1)];Y1<-Y1[!is.na(Z1)];Z1_1<-Z1_1[!is.na(Z1)];Z1<-Z1[!is.na(Z1)]

  x0_2<-result2$x0
  CI_average_2<-result2$CI_average
  X2<-x0_2$x1;Y2<-x0_2$x2;Z2<-CI_average_2[,1];Z2_1<-CI_average_2[,2]
  X2<-X2[!is.na(Z2)];Y2<-Y2[!is.na(Z2)];Z2_1<-Z2_1[!is.na(Z2)];Z2<-Z2[!is.na(Z2)]


  m<-300

  xo <- sort(c(range(X1),range(X2)))[2:3]
  yo <- sort(c(range(Y1),range(Y2)))[2:3]
  xo1<-range(quantile(treatment$x1,probs=Range),quantile(control$x1,probs=Range))
  yo1<-range(quantile(treatment$x2,probs=Range),quantile(control$x2,probs=Range))
  xo <- sort(c(xo,xo1))[2:3]
  yo <- sort(c(yo,yo1))[2:3]
  akima1_lower<-interp(X1,Y1,Z1,xo=seq(xo[1],xo[2],length=m),yo=seq(yo[1],yo[2],length=m+1))
  akima1_higher<-interp(X1,Y1,Z1_1,xo=seq(xo[1],xo[2],length=m),yo=seq(yo[1],yo[2],length=m+1))

  akima2_lower<-interp(X2,Y2,Z2,xo=seq(xo[1],xo[2],length=m),yo=seq(yo[1],yo[2],length=m+1))
  akima2_higher<-interp(X2,Y2,Z2_1,xo=seq(xo[1],xo[2],length=m),yo=seq(yo[1],yo[2],length=m+1))



  col<-colorRampPalette(c("red", "green"))(11)
  colorscale=list(c(0.0, col[1]),
                  list(0.1, col[2]),
                  list(0.2, col[3]),
                  list(0.3, col[4]),
                  list(0.4, col[5]),
                  list(0.5, col[6]),
                  list(0.6, col[7]),
                  list(0.7, col[8]),
                  list(0.8, col[9]),
                  list(0.9, col[10]),
                  list(1.0, col[11]))



  f2 <- list(
    family = "Old Standard TT, serif",
    size = 14,
    color="grey")

  scene1<-list(aspectmode="manual",aspectratio=list(x=1,y=1,z=0.95),
               xaxis=list(title=xlab,titlefont=f2,range=xo),
               yaxis=list(title=ylab,titlefont=f2,range=yo),
               zaxis=list(title=zlab,titlefont=f2,showline=TRUE),
               camera = list(eye = list(x = 2.25, y = -2.25, z = 0.9)))

  scene2<-list(aspectmode="manual",aspectratio=list(x=1,y=1,z=0.95),
               xaxis=list(title=xlab,titlefont=f2),
               yaxis=list(title=ylab,titlefont=f2),
               zaxis=list(title="Difference",titlefont=f2,showline=TRUE),
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

  b <- list(
    text = "Difference Plot",
    xref = "paper",
    yref = "paper",
    yanchor = "bottom",
    xanchor = "center",
    align = "center",
    x = 0.5,
    y = 0.92,
    showarrow = FALSE
  )

  ###Figure 1
  zo<-range(c(range(result1$CI_average,na.rm=TRUE),range(result2$CI_average,na.rm=TRUE)))
  zdiff<-diff(zo)
  col<-colorRampPalette(c("red","green"))(3)

  trace1<-plot_ly(showlegend=FALSE,scene="scene1")%>%
    add_surface(x=akima1_lower$x,y=akima1_lower$y,z=t(akima1_lower$z),colors=col,opacity=0.5,cmin=zo[1]-4*zdiff,cmax=zo[2],showscale=FALSE) %>%
    add_surface(x=akima1_higher$x,y=akima1_higher$y,z=t(akima1_higher$z),opacity = 0.5,cmin=zo[1]-4*zdiff,cmax=zo[2],showscale=FALSE)%>%
    add_fun(function(plot){plot%>%add_surface(x=akima2_lower$x,y=akima2_lower$y,z=t(akima2_lower$z),cmin=zo[1],cmax=zo[2]+4*zdiff,opacity=0.5,showscale=FALSE) %>%
        add_surface(x=akima2_higher$x,y=akima2_higher$y,z=t(akima2_higher$z),cmin=zo[1],cmax=zo[2]+4*zdiff,opacity = 0.5,showscale=FALSE)})%>%
    layout(annotations = a)



  ###Figure 2
  Diff<-all$Diff
  x0<-Diff$x0
  CI_average<-Diff$CI_average
  X<-x0$x1;Y<-x0$x2;Z<-CI_average[,1];Z2<-CI_average[,2]
  X<-X[!is.na(Z)];Y<-Y[!is.na(Z)];Z2<-Z2[!is.na(Z)];Z<-Z[!is.na(Z)]

  akima_lower<-interp(X,Y,Z,xo=seq(xo[1],xo[2],length=m),yo=seq(yo[1],yo[2],length=m+1))
  akima_higher<-interp(X,Y,Z2,xo=seq(xo[1],xo[2],length=m),yo=seq(yo[1],yo[2],length=m+1))
  x<-akima_lower$x
  y<-akima_lower$y
  z<-(akima_lower$z+akima_higher$z)/2

  #Identify where the insignificant area is
  which<-akima_lower$z<0 & akima_higher$z>0
  which[is.na(which)]<-FALSE
  z[which]<-0


  #gg-plot
  tmp<-expand.grid(x=x,y=y)
  pts<-tmp[as.vector(which),]
  temp<-sum(which)                                   #For the grey area
  tmp<-cbind(tmp,z=as.vector(z))
  tmp<-tmp[complete.cases(tmp),]                     #For geom_raster

  #dotted in grey area
  dotted<-expand.grid(x=seq(xo[1],xo[2],length=100),y=seq(yo[1],yo[2],length=100))
  tmp0<-tmp[sample(nrow(tmp),0.4*nrow(tmp)),]
  find<-interp(tmp0$x,tmp0$y,tmp0$z,xo=seq(xo[1],xo[2],length=100),yo=seq(yo[1],yo[2],length=100))
  which<-as.vector(find$z==0)
  which[is.na(which)]<-FALSE
  dotted<-dotted[which,]

  binwidth<-diff(range(tmp$z))/10
  vals<-abs(tmp$z)
  tmp2<-tmp;tmp2$z[vals>0]<-1                        #For contour plot
  lim<-max(vals)
  myPalette<-colorRampPalette(c("red","gray85","green"))
  sc <- scale_fill_gradientn(colours = myPalette(100), limits=c(-lim, lim))
  g<-ggplot() +
    geom_raster(data=tmp, aes(x,y,fill = z), interpolate = TRUE, na.rm=TRUE)+
    geom_tile(data=pts, aes(x=x,y=y),fill="grey70")+
    geom_point(data=dotted, aes(x=x,y=y),shape=46)+
    geom_contour(data=tmp2, aes(x, y, z = z), binwidth = 1, colour='black')+
    geom_contour(data=tmp, aes(x, y, z = z), binwidth = binwidth, colour='grey30')+
    coord_cartesian(expand=FALSE, xlim=xo,ylim=yo)+
    sc + labs(x=xlab,y=ylab,fill=paste("Significance\n   Legend"))+
    theme(panel.background = element_rect(fill = NA,colour = "grey50"),
          legend.position = "bottom", panel.grid.major.x = element_line(colour = "grey80"),
          panel.grid.major.y = element_line(colour = "grey80"))



  #Figure 3
  cmin<-min(akima_lower$z,na.rm=TRUE)
  cmax<-max(akima_higher$z,na.rm=TRUE)
  C<-max(c(abs(cmin),abs(cmax)))
  getYellow<-colorRampPalette(c("green","yellow"))(11)[8]
  col<-myPalette(100)

  zoffset<-cmin-(cmax-cmin)*0.1
  temp<-as.vector(akima_higher$z)
  temp[!is.na(temp)]<-zoffset

  trace2= plot_ly(showscale=FALSE,scene="scene2")%>%add_surface(x=akima_lower$x,y=akima_lower$y,z=t(akima_lower$z),colors=col,opacity = 0.7,cmin=-C,cmax=C)%>%
    add_surface(x=akima_higher$x,y=akima_higher$y,z=t(akima_higher$z),opacity = 0.7,colors=col,cmin=-C,cmax=C)%>%
    add_surface(x=akima_higher$x,y=akima_higher$y,z=t(matrix(0,nrow=nrow(akima_higher$z),ncol=ncol(akima_higher$z))),colorscale=colorscale,opacity=0.7)%>%
    add_surface(x=akima_higher$x,y=akima_higher$y,z=t(matrix(temp,nrow=nrow(akima_higher$z),ncol=ncol(akima_higher$z))),surfacecolor=t((akima_higher$z+akima_lower$z)/2),cmin=-C,cmax=C)%>%
    layout(annotations = b)


  trace<-subplot(trace1,trace2, titleX = TRUE, titleY = TRUE)%>%layout(title="Comparative Experiment",titlefont=list(family="Impact",size=16), scene=scene1,scene2=scene2)

  return(list(trace=trace, g=g))
}



CompExperiment_2dim_old<-function(cores=NULL,treatment=parent.frame(),control=parent.frame(),alpha=0.05,tau=0.5,lambda_trt=2,lambda_ctr=2,
                              b1=ceiling(nrow(treatment)^0.6), b2=ceiling(nrow(control)^0.6),s=20, r=200,figure=TRUE,
                              col_trt="green",col_ctr="red",col_mid="yellow",x1="x1",x2="x2",D=50){

  if(!all(colnames(treatment)==c("x1","x2","y"))|!all(colnames(control)==c("x1","x2","y"))) {colnames(treatment)<-c("x1","x2","y");colnames(control)<-c("x1","x2","y")}

  temp<-BLB_rqss_2dim(cores=cores,data=treatment,alpha=alpha,tau=tau,lambda=lambda_trt,b=b1,s=s,r=r,D=D)
  x0_1<-temp$x0
  CI_average_1<-temp$CI_average
  X1<-x0_1$x1;Y1<-x0_1$x2;Z1<-CI_average_1[,1];Z1_1<-CI_average_1[,2]
  X1<-X1[!is.na(Z1)];Y1<-Y1[!is.na(Z1)];Z1_1<-Z1_1[!is.na(Z1)];Z1<-Z1[!is.na(Z1)]
  cat("Run BLB for treatment sample\n")

  temp<-BLB_rqss_2dim(cores=cores,data=control,alpha=alpha,tau=tau,lambda=lambda_ctr,b=b1,s=s,r=r,D=D)
  x0_2<-temp$x0
  CI_average_2<-temp$CI_average
  X2<-x0_2$x1;Y2<-x0_2$x2;Z2<-CI_average_2[,1];Z2_1<-CI_average_2[,2]
  X2<-X2[!is.na(Z2)];Y2<-Y2[!is.na(Z2)];Z2_1<-Z2_1[!is.na(Z2)];Z2<-Z2[!is.na(Z2)]
  cat("Run BLB for control sample\n")

  xo<-range(range(treatment$x1),range(control$x1))
  yo<-range(range(treatment$x2),range(control$x2))
  akima1_lower<-interp(X1,Y1,Z1,xo=seq(xo[1],xo[2],length=200),yo=seq(yo[1],yo[2],length=200))
  akima1_higher<-interp(X1,Y1,Z1_1,xo=seq(xo[1],xo[2],length=200),yo=seq(yo[1],yo[2],length=200))
  # plot(X1,Y1,main="Sample 1")
  # image(akima1_lower,col=terrain.colors(12),add=TRUE)

  akima2_lower<-interp(X2,Y2,Z2,xo=seq(xo[1],xo[2],length=200),yo=seq(yo[1],yo[2],length=200))
  akima2_higher<-interp(X2,Y2,Z2_1,xo=seq(xo[1],xo[2],length=200),yo=seq(yo[1],yo[2],length=200))
  # plot(X2,Y2,main="Sample 2")
  # image(akima2_lower,col=terrain.colors(12),add=TRUE)

  I<-rep(NA,200*200)
  I[which(as.vector(akima1_lower$z-akima2_higher$z)>0)]<-1
  I[which(as.vector(akima2_lower$z-akima1_higher$z)>0)]<--1
  temp<-c(which(as.vector(akima1_lower$z-akima2_higher$z)>0),which(as.vector(akima2_lower$z-akima1_higher$z)>0))
  # all(!duplicated(temp))
  whole<-which(!is.na(as.vector(akima1_lower$z-akima2_higher$z)))
  rest<-whole[!(whole %in% temp)]
  I[rest]<-0
  temp<-sort(unique(I))
  I<-matrix(I,200,200)

  col<-rep(col_mid,length(temp))
  for(i in 1:length(temp)){
    if(temp[i]==1)  col[i]<-col_trt
    else if(temp[i]==-1) col[i]<-col_ctr
  }
  new_akima<-list(x=akima1_lower$x,y=akima1_lower$y,z=I)

  if(figure) {
    plot(X1,Y1,main="Comparative Experiment",type="n",xlab="x1",ylab="x2",xlim=xo,ylim=yo)
    image(new_akima,add=TRUE,col=col)
    contour(new_akima,add=TRUE,nlevels=3)
    legend("topleft",legend=c("treatment","Control","In-between"),pch=15,col=c(col_trt,col_ctr,col_mid))
  }
  return(new_akima)
}
