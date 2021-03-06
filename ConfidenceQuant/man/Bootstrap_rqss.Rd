% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Bootstrap_rqss.R
\name{Bootstrap_rqss}
\alias{Bootstrap_rqss}
\title{Regular Bootstrap Method.}
\usage{
Bootstrap_rqss(cores = NULL, data = parent.frame(), alpha = 0.05,
  tau = 0.25, lambda = 2, B = 100, D = 100, Search = FALSE,
  x0 = NULL, warning.catch = FALSE)
}
\arguments{
\item{cores}{The number of cores to use for parallel execution. If not
specified, the number of cores is set to the value of
\code{options("cores")}, if specified, or to one-half the number of cores
detected by the \code{\link{parallel}} package.}

\item{data}{A 2-dim optional data frame (or object coercible by
\code{\link{as.data.frame}} to a data frame) containing the variables in
the model. The column names should be specified in a way that “x” is for
the predictor, and “y” for the response.}

\item{alpha}{The confidence level required. The default is 0.05.}

\item{tau}{A specific quantile to be estimated. Must be a number between 0
and 1.}

\item{lambda}{The smoothing parameter governing the tradeoff between fidelity
and the penalty component for the triogram term. If \code{Search=TRUE},
there is no need for users to specify a value.}

\item{B}{The number of Monte Carlo iterations using bootstrap with
replacement. \eqn{B=100} is by default.}

\item{D}{A number that specifies how many x values you want to examine the
confidence bands at. If specified, it will look at \eqn{D} equal-distanced
points.}

\item{Search}{If \code{TRUE} (which is recommended), then the function will
first search for an optimum smoothing parameter \eqn{\lambda}.}
}
\value{
A list with two parts: \code{x0} and \code{CIs}, where \code{x0}
  contains the x values that we are examining the confidence intervals at,
  and \code{CIs} is 2-dim matrix which contains the corresponding lower bound
  and upper bound.
}
\description{
\code{Bootstrap_rqss}  Regular bootstrap for data sets with 1-dim covariate.
Used to generate the confidence bands for a quantile smoothing spline fitted
with the rqss function from package \code{\link{quantreg}}.
}
\details{
This function is based on a regular bootstrap method, which calculates
confidence bands for one quantile. It performs parallelization to speed up
the calculation.
}
\examples{
data(one)
result<-Bootstrap_rqss(data=one, alpha=0.05, tau=0.5, B=100, Search=TRUE)


plot<-FitPlot_1dim(data=one, result=result, xlab='x',ylab='y')
plot

}
\references{
Micheal, I. J et al.(2012). \eqn{A Scalable Bootstrap for Massive
  Data}.
}
\seealso{
\code{\link{BLB_rqss}} for BLB with 1-dim covariate data set.

\code{\link{Bootstrap_rqss_2dim}} for regular Bootstrap method with
  2-dim covariate data set.
}
