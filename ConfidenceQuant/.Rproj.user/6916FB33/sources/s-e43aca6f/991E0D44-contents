#'Comparative Experiment.
#'
#'\code{CompExperiment_1dim} compares the confidence bands at a given quantile
#'for two different datasets, one related to a treatment and the other to a
#'control. It applies the BLB algorithm to each dataset to get confidence bands
#'using quantile smoothing splines for 1-dim covariate.
#'
#'This function runs \code{BLB} twice, once for each dataset. It is based on
#'\code{\link{BLB_rqss}}, which implements BLB for quantile smoothing splines
#'with a one dimensional covariate dataset. It performs parallelization to speed
#'up the calculation.
#'
#'\if{html}{\figure{comp.png}{options: width=100 alt="Image output"}}
#'\if{latex}{\figure{comp.png}{options: width=3in}}
#'
#'\code{CompPlot_1dim} takes the results and use ggplot/plotly to visualize
#'them, in which different colors represent different scenarios. See figure
#'above.
#'
#'@import quantreg
#'@import doParallel
#'@import foreach
#'@import rootSolve
#'@import scales
#'@importFrom ggplot2 ggplot geom_point geom_ribbon geom_line labs theme geom_hline
#'    geom_rect coord_cartesian scale_fill_manual aes element_text element_blank element_rect element_line
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
#'
#'  In \code{CompPlot_1dim} this stands for how many percentage of data you want
#'  to throw away for the final plot.
#'@param tau A specific quantile to be estimated. Must be a number between 0 and
#'  1.
#'@param lambda The smoothing parameter governing the tradeoff between fidelity
#'  and the penalty component for the triogram term. If \code{Search=TRUE},
#'  there is no need for users to specify a value.
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
#'@param xlab,ylab Titles for x and y axis. See \code{\link{title}}.
#'@param Search If \code{TRUE} (which is recommended), then the function will
#'  first search for an optimum smoothing parameter \eqn{\lambda}.
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
#'  \code{BLB_rqss}. If it is done automatically, the function also returns
#'  \code{Lambda} and \code{Fid}, which respectively stand for a vector of lambda
#'  values and their corresponding cross-validation MCV values.
#'
#' @examples
#' data(treatment)
#' data(control)
#'
#' #alpha=0.05;tau=0.5
#' all<-CompExperiment_1dim(cores=7, treatment, control, tau=0.5, Search=TRUE)
#'
#' plot<-CompPlot_1dim(treatment = treatment, control=control, all = all, xlab = 'x', ylab = 'y')
#'
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

CompExperiment_1dim<-function(cores=NULL,treatment,control,alpha=0.05,tau=0.25,lambda=2,D=100,b1=ceiling(nrow(treatment)^0.6),b2=ceiling(nrow(control)^0.6), s=15, r=100, Search=FALSE){
  if(is.null(cores)) {cores=round(0.5*detectCores())}
  registerDoParallel(cores=cores)

  if(!all(colnames(treatment)==c("x","y"))) colnames(treatment)<-c("x","y")
  if(!all(colnames(control)==c("x","y"))) colnames(control)<-c("x","y")

  indicator<-(Search|length(lambda)>1)

  ###1. Sample s disjoint subsets with size b
  n<-nrow(treatment)
  shuffle<-sample(1:n)
  indices1<-matrix(head(shuffle,b1*s),nrow=s,byrow=TRUE)
  Rboot1<-rmultinom(s,size=n,prob=rep(1/b1,b1))

  m<-nrow(control)
  shuffle<-sample(1:m)
  indices2<-matrix(head(shuffle,b2*s),nrow=s,byrow=TRUE)
  Rboot2<-rmultinom(s,size=m,prob=rep(1/b2,b2))


  ###2. λ Search
  r1<-list(lambda=lambda)
  r2<-list(lambda=lambda)

  if(indicator){
    r1<-searchLambda(lambda=lambda,data=treatment,s=s,indices=indices1,tau=tau,Rboot=Rboot1)
    r2<-searchLambda(lambda=lambda,data=control,s=s,indices=indices2,tau=tau,Rboot=Rboot2)
    #plot(r$Lambda[-1],r$Fidelity[-1],type='l',ylab="New Metric",xlab="λ",main="MCV for BLB repeated data")
  }


  ###3. Bag of Little Bootstrap
  xrange<-sort(c(range(treatment$x), range(control$x)))[2:3]
  x0 <- seq(xrange[1], xrange[2], by = diff(xrange)/D)

  response_trt <- rep(NA,length(x0))   #Prediction for one bootstrap sample_treatment
  Response_trt<-matrix(nrow=length(x0),ncol=r)

  response_ctr <- rep(NA,length(x0))   #Prediction for one bootstrap sample_control
  Response_ctr<-matrix(nrow=length(x0),ncol=r)

  Difference<-matrix(nrow=length(x0),ncol=r)

  CIs<-foreach (j=1:s, .inorder = FALSE) %dopar% {
    subsample1<-treatment[indices1[j,],]
    subsample2<-control[indices2[j,],]

    extrapolated1 = x0 > max(subsample1$x) | x0 < min(subsample1$x)
    extrapolated2 = x0 > max(subsample2$x) | x0 < min(subsample2$x)

    for(k in 1:r){
      response_trt[1:length(x0)] <- NA
      response_ctr[1:length(x0)] <- NA
      rboot1<-as.vector(rmultinom(1,size=n,prob=rep(1/b1,b1)))
      rboot2<-as.vector(rmultinom(1,size=m,prob=rep(1/b2,b2)))

      mod1 <- rqss_new(y ~ qss(x, lambda = r1$lambda), tau = tau, data = subsample1, weights = rboot1)
      response_trt[!extrapolated1] = predict(mod1, newdata = data.frame(x=x0[!extrapolated1]))
      Response_trt[,k]<-response_trt

      mod2 <- rqss_new(y ~ qss(x, lambda = r2$lambda), tau = tau, data = subsample2, weights = rboot2)
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

#'\code{searchLambda} is a wrapper function to calculate the optimum lambda.
#' @name searchLambda
#' @rdname CompExperiment_1dim
#' @export


searchLambda<-function(lambda,data,s,indices,tau,Rboot){
  if(length(lambda)>1) Lambda<-lambda else Lambda<-c(seq(0,400,by=10),500,1000,1500,2000,3000)
  fid<-rep(NA,length(Lambda))

  Fid<-foreach(j=1:s,.combine = "cbind", .inorder = FALSE) %dopar% {

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
  return(list(lambda=lambda,Lambda=Lambda,Fidelity=Fidelity))
}


#'\code{CompPlot_1dim} is a function to plot the results.
#' @name CompPlot_1dim
#' @rdname CompExperiment_1dim
#' @export

CompPlot_1dim<-function(treatment,control,all,xlab="x",ylab="y",alpha=0.05){
  if(!all(colnames(treatment)==c("x","y"))) colnames(treatment)<-c("x","y")
  if(!all(colnames(control)==c("x","y"))) colnames(control)<-c("x","y")
  result1<-all$result1
  result2<-all$result2
  n1<-nrow(treatment)
  n2<-nrow(control)
  yo<-range(c(quantile(treatment$y,0.001),quantile(treatment$y,0.999)),c(quantile(control$y,0.001),quantile(control$y,0.999)))

  #To solve the two functions, we need to extrapolate.
  xo_inner<-sort(c(range(result1$x0[!is.na(result1$CI_average[,1])]),range(result2$x0[!is.na(result2$CI_average[,1])])))[2:3]
  lower1<-approxfun(result1$x0,result1$CI_average[,1], yleft=result1$CI_average[min(which(!is.na(result1$CI_average[,1]))),1], yright=result1$CI_average[max(which(!is.na(result1$CI_average[,1]))),1])
  upper1<-approxfun(result1$x0,result1$CI_average[,2], yleft=result1$CI_average[min(which(!is.na(result1$CI_average[,1]))),2], yright=result1$CI_average[max(which(!is.na(result1$CI_average[,1]))),2])
  lower2<-approxfun(result2$x0,result2$CI_average[,1], yleft=result2$CI_average[min(which(!is.na(result2$CI_average[,1]))),1], yright=result2$CI_average[max(which(!is.na(result2$CI_average[,1]))),1])
  upper2<-approxfun(result2$x0,result2$CI_average[,2], yleft=result2$CI_average[min(which(!is.na(result2$CI_average[,1]))),2], yright=result2$CI_average[max(which(!is.na(result2$CI_average[,1]))),2])

  ##1. Get the roots

  roots_1<-uniroot.all(function(x) upper1(x)-lower2(x),xo_inner,n=1000)
  roots_2<-uniroot.all(function(x) upper2(x)-lower1(x),xo_inner,n=1000)

  #--Change it back without extropolate.
  lower1<-approxfun(result1$x0,result1$CI_average[,1])
  upper1<-approxfun(result1$x0,result1$CI_average[,2])
  lower2<-approxfun(result2$x0,result2$CI_average[,1])
  upper2<-approxfun(result2$x0,result2$CI_average[,2])


  ##2. Decide who wins
  interests<-sort(c(xo_inner[1],roots_1,roots_2,xo_inner[2]))
  grid<-rep("Tie",length(interests)-1)
  for(i in 1:(length(interests)-1)){
    temp<-mean(c(interests[i],interests[i+1]))
    if(lower1(temp)-upper2(temp)>0) grid[i]<-"Treatment>Control"
    if(lower2(temp)-upper1(temp)>0) grid[i]<-"Control>Treatment"
  }
  grid<-factor(grid,levels=c("Treatment>Control","Tie","Control>Treatment" ))
  color<-c("dodgerblue3","grey40","indianred2")


  ##3. gg-plot
  dat1<-data.frame(x0=result1$x0,lower=result1$CI_average[,1],upper=result1$CI_average[,2])
  dat2<-data.frame(x0=result2$x0,lower=result2$CI_average[,1],upper=result2$CI_average[,2])

  #--Define xlim, ylim
  valid<-!is.na(all$Diff$CI_average[,1])
  xo<-range(result1$x0[valid])
  ylim<-c(yo[1]-0.035*diff(yo),yo[2])

  #--Define Rectangle
  xmin<-c(xo[1],sort(c(roots_1,roots_2)))
  xmax<-c(sort(c(roots_1,roots_2)),xo[2])
  ymax<-rep(yo[1],length(xmin))
  ymin<-rep(yo[1]-0.035*diff(yo),length(xmin))
  rec<-data.frame(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax)

  cole<-factor(c("Treatment>Control","Tie","Control>Treatment"))
  Grid<-as.character(grid)
  for(i in 1:length(cole)){
    add<-tail(rec,1)
    add[1:2]<-add[1:2]+10000
    rec<-rbind(rec,add)
    Grid<-c(Grid,as.character(cole[i]))
  }
  Grid<-factor(Grid,levels=c("Treatment>Control","Tie","Control>Treatment" ))

  #--Define plot
  h<-ggplot(treatment[sample(n1,max(min(n1*0.001,5000),3000)),])
  h<-h+geom_point(aes(x=x,y=y),alpha=0.09,colour="lightblue")+
    geom_point(data=control[sample(n2,max(min(n2*0.001,5000),3000)),],aes(x=x,y=y),alpha=0.045,colour="red")+
    geom_ribbon(data=dat1,aes(x=x0, ymin = lower, ymax = upper), alpha=0.7, fill = "dodgerblue3")+
    geom_ribbon(data=dat2,aes(x=x0, ymin = lower, ymax = upper), alpha=0.7, fill = "red")+
    geom_line(data=dat1,aes(x=x0,y=lower),alpha=0.7, colour="dodgerblue3")+
    geom_line(data=dat1,aes(x=x0,y=upper),alpha=0.7, colour="dodgerblue3")+
    geom_line(data=dat2,aes(x=x0,y=lower),alpha=0.7, colour="red")+
    geom_line(data=dat2,aes(x=x0,y=upper),alpha=0.7, colour="red")+
    geom_rect(data=rec, aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,fill = Grid), colour = "grey70")+
    coord_cartesian(xlim=xo,ylim = ylim,expand=FALSE)+
    scale_fill_manual(values=color)+
    labs(title="Comparative Experiment",x=xlab,y=ylab)+
    theme(plot.title = element_text(hjust = 0.5, face="bold"),legend.title=element_blank(),legend.position = c(.99, .99), legend.justification = c("right","top"))

  ##4. Difference plot
  dat<-data.frame(x0=all$Diff$x0,ymin=all$Diff$CI_average[,1],ymax=all$Diff$CI_average[,2])
  dat<-dat[valid,]
  Ylim<-range((dat$ymin+dat$ymax)/2)
  Ylim[1]<-Ylim[1]-0.01*diff(yo)
  Ylim[2]<-Ylim[2]+0.01*diff(yo)

  temp<-c(treatment$x,control$x)
  left<-quantile(temp,alpha/2)
  right<-quantile(temp,1-alpha/2)

  d1<-dat[which(dat$x0<=left),]
  d2<-dat[which(dat$x0>left&dat$x0<=right),]
  d3<-dat[which(dat$x0>=right),]

  temp1<-tail(d1,1)
  temp2<-head(d2,1)
  conceited<-diff(c(temp1$x0,temp2$x0))/10

  med1<-data.frame(x0=left-conceited, ymin=approx(x=c(temp1$x0,temp2$x0),y=c(temp1$ymin,temp2$ymin),xout=left-conceited,rule=2)$y, ymax=approx(x=c(temp1$x0,temp2$x0),y=c(temp1$ymax,temp2$ymax),xout=left-conceited,rule=2)$y)
  med2<-data.frame(x0=left, ymin=approx(x=c(temp1$x0,temp2$x0),y=c(temp1$ymin,temp2$ymin),xout=left)$y, ymax=approx(x=c(temp1$x0,temp2$x0),y=c(temp1$ymax,temp2$ymax),xout=left)$y)

  temp3<-tail(d2,1)
  temp4<-head(d3,1)
  conceited<-diff(c(temp3$x0,temp4$x0))/10

  med3<-data.frame(x0=right, ymin=approx(x=c(temp3$x0,temp4$x0),y=c(temp3$ymin,temp4$ymin),xout=right)$y, ymax=approx(x=c(temp3$x0,temp4$x0),y=c(temp3$ymax,temp4$ymax),xout=right)$y)
  med4<-data.frame(x0=right+conceited, ymin=approx(x=c(temp3$x0,temp4$x0),y=c(temp3$ymin,temp4$ymin),xout=right+conceited,rule=2)$y, ymax=approx(x=c(temp3$x0,temp4$x0),y=c(temp3$ymax,temp4$ymax),xout=right+conceited,rule=2)$y)

  if(med1$x0>temp1$x0) d1<-rbind(d1,med1)
  d2<-rbind(med2,d2,med3)
  if(med4$x0<temp4$x0) d3<-rbind(med4,d3)

  ##Futher divide d2
  Which <- interests < right
  d2$col <- grid[max(which(Which))]
  myList<-list()
  myList[[1]]<-d2

  which <- interests > left & interests < right
  r <- interests[which]
  if (length(r) > 1) {
    which[max(which(which)) + 1] <- TRUE
    which <- which[-1]
    grid <- grid[which]
    tmp1 <- approx(x = d2$x0, y = d2$ymin, xout = r)
    tmp2 <- approx(x = d2$x0, y = d2$ymax, xout = r)
    myList<-list()

    old.ind <- 1
    for (i in 1:length(r)) {
      ind <- max(which(d2$x0 < r[i]))
      add <- data.frame(x0 = r[i], ymin = tmp1$y[i], ymax = tmp2$y[i],col = tail(grid, 1))
      d2 <- rbind(d2[1:ind, ], add, add, d2[(ind + 1):nrow(d2),])
      d2$col[old.ind:(ind + 1)] <- grid[i]
      myList[[i]]<-d2[old.ind:(ind + 1),]
      old.ind <- ind + 2
    }
    myList[[length(r)+1]]<-d2[old.ind:nrow(d2),]
  }

  color2 <- c("blue", "grey20", "red")
  g <- ggplot()+
    geom_ribbon(data = d1, aes(x = x0,ymin = ymin, ymax = ymax), alpha = 0.7, fill = "grey70")+
    geom_ribbon(data = d3, aes(x = x0, ymin = ymin, ymax = ymax), alpha = 0.7, fill = "grey70")+
    geom_line(data=d1,aes(x=x0,y=ymin),alpha=0.7,colour="grey70")+geom_line(data=d1,aes(x=x0,y=ymax),alpha=0.7,colour="grey70")+
    geom_line(data=d2,aes(x=x0,y=ymin),alpha=0.7,colour="grey70")+geom_line(data=d2,aes(x=x0,y=ymax),alpha=0.7,colour="grey70")

  cols<-rep(NA,length(myList))
  for(i in 1:length(myList)){
    temp<-myList[[i]]
    cols[i]<-color2[grep(temp$col[1],levels(temp$col[1]))]
  }

  for(i in 1:length(myList)){
    g<-g+geom_ribbon(data=myList[[i]], aes(x = x0, ymin = ymin, ymax = ymax), fill = cols[i], alpha = 0.7)+
      geom_line(data=myList[[i]],aes(x=x0,y=ymin),alpha=0.7,colour=cols[i])+geom_line(data=myList[[i]],aes(x=x0,y=ymax),alpha=0.7,colour=cols[i])
  }

  g<-g+geom_hline(yintercept = 0, color = "darkgreen", linetype = 1, size = 1) +
    coord_cartesian(xlim = xo, ylim=Ylim,expand = FALSE) + labs(y = "Difference") +
    theme(axis.title.x = element_blank(), panel.background = element_rect(fill = NA,colour = "grey50"),
          legend.position = "none", panel.grid.major.x = element_line(colour = "grey80"),
          panel.grid.minor = element_blank(), panel.grid.major.y = element_blank())

  return(list(h=h,g=g))
}







CompExperiment_1dim_old<-function(cores=NULL,treatment=parent.frame(),control=parent.frame(),alpha=0.05,tau=0.5,lambda_trt=2,lambda_ctr=2,
                              b1=ceiling(nrow(treatment)^0.6), b2=ceiling(nrow(control)^0.6),s=20, r=100,figure=TRUE,
                              col_trt="green",col_ctr="red",col_mid="yellow",xlab="Log(throughput[scaled])",ylab="Log(playdelay[scaled]",D=100){

  if(!all(colnames(treatment)==c("x","y"))|!all(colnames(control)==c("x","y"))) {colnames(treatment)<-c("x","y");colnames(control)<-c("x","y")}

  temp<-BLB_rqss(cores=cores, treatment, alpha=alpha, tau=tau, b=b1, s=s, r=r, lambda=lambda_trt,D=D)
  x0_1<-temp$x0
  CI_average_1<-temp$CI_average
  cat("Run BLB for treatment sample\n")

  temp<-BLB_rqss(cores=cores, control, alpha=alpha, tau=tau, b=b2, s=s, r=r, lambda=lambda_ctr,D=D)
  x0_2<-temp$x0
  CI_average_2<-temp$CI_average
  cat("Run BLB for control sample\n")

  xo<-range(range(treatment$x),range(control$x))
  yo<-range(range(treatment$y),range(control$y))
  xo_inner<-sort(c(range(x0_1[!is.na(CI_average_1[,1])]),range(x0_2[!is.na(CI_average_2[,1])])))[2:3]

  #To solve the two functions, we need to extrapolate.
  lower1<-approxfun(x0_1,CI_average_1[,1], yleft=CI_average_1[min(which(!is.na(CI_average_1[,1]))),1], yright=CI_average_1[max(which(!is.na(CI_average_1[,1]))),1])
  upper1<-approxfun(x0_1,CI_average_1[,2], yleft=CI_average_1[min(which(!is.na(CI_average_1[,1]))),2], yright=CI_average_1[max(which(!is.na(CI_average_1[,1]))),2])
  lower2<-approxfun(x0_2,CI_average_2[,1], yleft=CI_average_2[min(which(!is.na(CI_average_2[,1]))),1], yright=CI_average_2[max(which(!is.na(CI_average_2[,1]))),1])
  upper2<-approxfun(x0_2,CI_average_2[,2], yleft=CI_average_2[min(which(!is.na(CI_average_2[,1]))),2], yright=CI_average_2[max(which(!is.na(CI_average_2[,1]))),2])

  #Get the roots
  roots_1<-uniroot.all(function(x) upper1(x)-lower2(x),xo_inner,n=1000)
  roots_2<-uniroot.all(function(x) upper2(x)-lower1(x),xo_inner,n=1000)

  #Change it back without extropolate.
  lower1<-approxfun(x0_1,CI_average_1[,1])
  upper1<-approxfun(x0_1,CI_average_1[,2])
  lower2<-approxfun(x0_2,CI_average_2[,1])
  upper2<-approxfun(x0_2,CI_average_2[,2])


  interests<-sort(c(xo_inner[1],roots_1,roots_2,xo_inner[2]))
  grid<-rep(0,length(interests)-1)
  for(i in 1:(length(interests)-1)){
    temp<-mean(c(interests[i],interests[i+1]))
    if(lower1(temp)-upper2(temp)>0) grid[i]<-1
    if(lower2(temp)-upper1(temp)>0) grid[i]<--1
  }

  result<-stepfun(tail(interests,length(interests)-1),c(grid,tail(grid,1)))
  # curve(result,from=-3,to=3)

  if(figure){
    temp<-max(quantile(control$y,0.99),quantile(treatment$y,0.99))
    id1<-sample(nrow(control),0.4*nrow(control))
    id2<-sample(nrow(treatment),0.4*nrow(treatment))
    plot(control$x[id1], control$y[id1], pch=20, col=alpha("red",0.02), xlab=xlab, ylab=ylab,ylim=c(yo[1],temp), main="Comparative Experiment")
    points(treatment$x[id2], treatment$y[id2], pch=20, col=alpha("green",0.02))

    curve(lower2,add=TRUE,lwd=2,lty=2,from=min(x0_2))
    curve(upper2,add=TRUE,lwd=2,lty=2,from=min(x0_2))
    curve(lower1,add=TRUE,lwd=2,from=min(x0_1))
    curve(upper1,add=TRUE,lwd=2,from=min(x0_1))

    if(length(roots_1)!=0){
      for(i in 1:length(roots_1)) lines(x=rep(roots_1[i],2),y=c(yo[1]-10,upper1(roots_1[i])),lty=3)
    }

    if(length(roots_2)!=0){
      for(i in 1:length(roots_2)) lines(x=rep(roots_2[i],2),y=c(yo[1]-10,upper2(roots_2[i])),lty=3)
    }

    abline(h=yo[1],col="grey")

    temp<-sort(unique(grid))
    col<-rep(col_mid,length(temp))
    for(i in 1:length(temp)){
      if(temp[i]==1)  col[i]<-col_trt
      else if(temp[i]==-1) col[i]<-col_ctr
    }

    grid=matrix(grid)
    temp1<-0.1*diff(xo)
    temp2<-0.1*diff(yo)
    image(c(xo[1]-temp1, sort(c(roots_1,roots_2)), xo[2]+temp1), c(yo[1]-temp1,yo[1]), grid, add=TRUE, col=col)
    legend("topright",legend=c("Treatment","Control","In-between"),pch=15,col=c(col_trt,col_ctr,col_mid))
  }
  return(result)
}
