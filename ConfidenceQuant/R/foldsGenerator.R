#'Fold Generator for Cross-Validation.
#'
#'\code{foldsGenerator} is used to generate folds for a specific data set when
#'performing crossvalidation.
#'
#'This function first identifies the number of observations in the data set.
#'With the pre-specified number of folds \code{nfolds}, the function generates a
#'sequence of randoms numbers which vary from \eqn{1} to \eqn{nfolds}. Each
#'number indicates which fold an observation belongs to.
#'
#'
#'
#'@param sim.data A data frame (or object coercible by
#'  \code{\link{as.data.frame}} to a data frame) on which the cross validation
#'  is performed.
#'@param nfolds The number of folds for the cross validation.
#'@param n The number of rows if the data frame not provided.
#'
#'
#'@return A vector of random numbers which vary from \eqn{1} to \eqn{nfolds}.
#' @examples
#' data(one)
#'
#' folds<-foldsGenerator(one,nfolds=10)
#' @export

foldsGenerator<-function(sim.data,nfolds=10,n=NULL){
  if(!is.null(n)) {data(sim.data)} else {n<-nrow(sim.data);data(sim.data1)}
  
  cv.idx<-as.numeric(matrix(1:nfolds,nrow=n,ncol=1))
  cv.idx<-sample(cv.idx)
  
  return(cv.idx)
}

