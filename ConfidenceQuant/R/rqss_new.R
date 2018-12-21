#' Additive Quantile Regression Smoothing
#'
#' \code{rqss_new} is based on \code{\link{rqss}} function in
#' \code{\link{quantreg}} package. The main difference is that this new function
#' implements the weights option, which can be quite useful. It is used to fit
#' function for additive quantile regression models with possible univariate
#' and/or bivariate nonparametric terms estimated by total variation
#' regularization.
#'
#' Total variation regularization for univariate and bivariate nonparametric
#' quantile smoothing is described in Koenker, Ng and Portnoy (1994) and Koenker
#' and Mizera(2003) respectively.  The additive model extension of this approach
#' depends crucially on the sparse linear algebra implementation for R described
#' in Koenker and Ng (2003).  There are extractor methods \code{\link{logLik}}
#' and \code{\link{AIC}} that is relevant to lambda selection.  A more detailed
#' description of some recent developments of these methods is available from
#' within the package with \code{vignette("rqss")}.  Since this function uses
#' sparse versions of the interior point algorithm it may also prove to be
#' useful for fitting linear models without \code{\link{qss}} terms when the
#' design has a sparse structure, as for example when there is a complicated
#' factor structure.
#'
#' If the \pkg{MatrixModels} and \pkg{Matrix} packages are both loadable then
#' the linear in parameters portion of the design matrix is made in sparse
#' matrix form, this is helpful in large applications with many factor variables
#' for which dense formation of the design matrix would take too much space.
#'
#' @import quantreg
#' @import tripack
#'
#' @param formula a formula object, with the response on the left of a `~'
#'   operator,  and terms, separated by `+' operators, on the right. The terms
#'   may include \code{qss} terms that represent additive nonparametric
#'   components.  These terms can be univariate or bivariate.  See
#'   \code{\link{qss}} for details on how to specify these terms.
#'
#' @param tau the quantile to be estimated, this must be a number between 0 and
#'   1.
#' @param data a data.frame in which to interpret the variables named in the
#'   formula, or in the subset and the weights argument.
#' @param weights vector of observation weights. The length of weights must be
#'   the same as the number of observations.  The weights must be positive
#'   integers which basically indicates the number of repeats.
#' @param na.action a function to filter missing data. This is applied to the
#'   model.frame after any subset argument has been used. The default (with
#'   \code{na.fail}) is to create an error if any missing values are found.  A
#'   possible alternative is \code{na.omit}, which deletes observations that
#'   contain one or more missing values.
#' @param method the algorithmic method used to compute the fit.  There are
#'   currently two options.   Both are implementations of the Frisch--Newton
#'   interior point method described in detail in Portnoy and Koenker(1997).
#'   Both are implemented using sparse Cholesky decomposition as described in
#'   Koenker and Ng (2003).
#'
#'   Option \code{"sfnc"} is used if the user specifies inequality constraints.
#'   Option \code{"sfn"} is used if there are no inequality constraints. Linear
#'   inequality constraints on the fitted coefficients are specified by a matrix
#'   \code{R} and a vector \code{r}, specified inside the \code{qss} terms,
#'   representing the constraints in the form \eqn{Rb \ge r}{Rb >= r}.
#'
#'   The option \code{method = "lasso"} allows one to penalize the coefficients
#'   of the covariates that have been entered linearly as in
#'   \code{\link{rq.fit.lasso}}; when this is specified then there should be an
#'   additional \code{lambda} argument specified that determines the amount of
#'   shrinkage.
#' @param lambda can be either a scalar, in which case all the slope
#'   coefficients are assigned this value, or alternatively, the user can
#'   specify a vector of length equal to the number of linear covariates plus
#'   one (for the intercept) and these values will be used as coordinate
#'   dependent shrinkage factors.
#' @param contrasts  a list giving contrasts for some or all of the factors
#'   default = \code{NULL} appearing in the model formula. The elements of the
#'   list should have the same name as the variable and should be either a
#'   contrast matrix (specifically, any full-rank matrix with as many rows as
#'   there are levels in the factor), or else a function to compute such a
#'   matrix given the number of levels.
#' @param ztol A zero tolerance parameter used to determine the number of
#' zero residuals in the fitted object which in turn determines the effective
#' dimensionality of the fit.
#' @param control control argument for the fitting routines
#' (see \code{\link{sfn.control}}
#' @param ... Other arguments passed to fitting routines
#'
#'
#' @return The function returns a fitted object representing the estimated
#' model specified in the formula.  See \code{\link{rqss.object}}
#' for further details on this object, and references to methods
#' to look at it.
#'
#'
#'
#' @references Koenker, R. and S. Portnoy (1997)
#' The Gaussian Hare and the Laplacean
#' Tortoise:  Computability of Squared-error vs Absolute Error Estimators,
#' (with discussion).
#' \emph{Statistical Science} \bold{12}, 279--300.
#' @references Koenker, R., P. Ng and S. Portnoy, (1994)
#' Quantile Smoothing Splines;
#' \emph{Biometrika} \bold{81}, 673--680.
#' @references Koenker, R. and I. Mizera, (2003)
#' Penalized Triograms: Total Variation Regularization for Bivariate Smoothing;
#' \emph{JRSS(B)} \bold{66}, 145--163.
#' @references Koenker, R. and P. Ng (2003)
#' SparseM:  A Sparse Linear Algebra Package for R,
#' \emph{J. Stat. Software}.
#' @references Esmond G. Ng and Barry W. Peyton, "Block sparse Cholesky algorithms
#' on advanced uniprocessor computers". SIAM J. Sci. Stat. Comput.
#' 14  (1993), pp. 1034-1056.
#' @references John R. Gilbert, Esmond G. Ng, and Barry W. Peyton, "An efficient
#' algorithm to compute row and column counts for sparse Cholesky
#' factorization". SIAM J. Matrix Anal. Appl. 15 (1994), pp. 1075-1091.
#'
#' @author Roger Koenker, Modified by Likun Zhang
#'
#' @seealso \code{\link{rqss}}, \code{\link{qss}}
#' @examples
#' One<-one[sample(1:2999,2000),]
#' tau<-0.5
#' fit0<-rqss(y ~ qss(x, lambda = 1), tau = tau, data = One)
#'
#'
#' x0<-seq(2.3,5.1,length=200)
#' response<-rep(NA,200)
#' extrapolated <- x0 > max(One$x) | x0 < min(One$x)
#' response[!extrapolated] = predict(fit0, newdata = data.frame(x=x0[!extrapolated]))
#'
#' plot(One,col="grey",pch=20,main="Implement Weights in rqss()")
#' lines(x0,response,lty=2,lwd=3)
#'
#'
#' times<-sample(1:30,2000,replace=TRUE)
#' index<-rep(1:2000,times=times)
#' One_new<-One[index,]
#' fit_rep<-rqss(y ~ qss(x, lambda = 1), tau = tau, data = One_new)
#'
#' response1<-rep(NA,200)
#' extrapolated <- x0 > max(One_new$x) | x0 < min(One_new$x)
#' response1[!extrapolated] = predict(fit_rep, newdata = data.frame(x=x0[!extrapolated]))
#'
#' lines(x0,response1,col="red",lwd=2)
#' legend("topright",lty=c(2,1),lwd=c(3,2),col=c("black","red"),legend=c("Original Data","Repeated Data"))
#' @keywords regression, smooth, robust
#'
#'@export


rqss_new<-function (formula, tau = 0.5, data = parent.frame(), weights,
          na.action, method = "sfn", lambda = NULL, contrasts = NULL,
          ztol = 1e-05, control = sfn.control(), ...)
{
  call <- match.call()
  m <- match.call(expand.dots = FALSE)
  temp <- c("", "formula", "data", "weights", "na.action")
  m <- m[match(temp, names(m), nomatch = 0)]
  m[[1]] <- as.name("model.frame")
  special <- "qss"
  Terms <- if (missing(data))
    terms(formula, special) else terms(formula, special, data = data)
  qssterms <- attr(Terms, "specials")$qss
  dropx <- NULL
  if (length(tau) > 1) {
    tau <- tau[1]
    warning("multiple taus not supported, using first element")
  }
  if (length(qssterms)) {
    tmpc <- untangle.specials(Terms, "qss")
    ord <- attr(Terms, "order")[tmpc$terms]
    if (any(ord > 1))
      stop("qss can not be used in an interaction")
    dropx <- tmpc$terms
    if (length(dropx))
      Terms <- Terms[-dropx]
    attr(Terms, "specials") <- tmpc$vars
    fnames <- function(x) {
      fy <- all.names(x[[2]])
      if (fy[1] == "cbind")
        fy <- fy[-1]
      fy
    }
    fqssnames <- unlist(lapply(parse(text = tmpc$vars), fnames))
    qssnames <- unlist(lapply(parse(text = tmpc$vars), function(x) deparse(x[[2]])))
  }
  m$formula <- Terms
  ff <- delete.response(terms(formula(Terms)))
  if (exists("fqssnames")) {
    ffqss <- paste(fqssnames, collapse = "+")
    ff <- paste(deparse(formula(ff)), "+", ffqss)
  }
  m <- eval(m, parent.frame())


  weights <- model.extract(m, weights)   ###Verified the weights are still here!
  process <- (tau < 0 || tau > 1)
  Y <- model.extract(m, "response")
  if (requireNamespace("MatrixModels") && requireNamespace("Matrix")) {
    X <- MatrixModels::model.Matrix(Terms, m, contrasts=NULL,
                                    sparse = TRUE)
    vnames <- dimnames(X)[[2]]
    X <- as(X, "matrix.csr")
  }else {
    X <- model.matrix(Terms, m, contrasts)
    vnames <- dimnames(X)[[2]]
  }

  p <- ncol(X)   ###X is the parametric part so far
  pf <- environment(formula)
  nrL <- 0


  if (method == "lasso") {
    if (!length(lambda))
      stop("No lambda specified for lasso constraint")
    if (length(lambda) == 1)
      lambda <- c(0, rep(lambda, p - 1))
    if (length(lambda) != p)
      stop("lambda must be either of length p, or length 1")
    if (any(lambda < 0))
      stop("negative lambdas disallowed")
    L <- diag(lambda, nrow = length(lambda))
    L <- L[which(lambda != 0), , drop = FALSE]
    L <- as.matrix.csr(L)
    nrL <- nrow(L)
    ncL <- ncol(L)
  }


  ####X is still the eparametric part
  if (length(qssterms) > 0) {
    F <- as.matrix.csr(X)

    qss <- lapply(tmpc$vars, function(u) eval(parse(text = u),
                                              data, enclos = pf))     ####Very important

    mqss <- length(qss)
    ncA <- rep(0, mqss + 1)
    nrA <- rep(0, mqss + 1)
    nrR <- rep(0, mqss + 1)
    for (i in 1:mqss) {
      F <- cbind(F, qss[[i]]$F)
      ncA[i + 1] <- ncol(qss[[i]]$A)
      nrA[i + 1] <- nrow(qss[[i]]$A)
      nrR[i + 1] <- ifelse(is.null(nrow(qss[[i]]$R)), 0,
                           nrow(qss[[i]]$R))
      vnames <- c(vnames, paste(qssnames[i], 1:ncA[i + 1], sep = ""))
    }

    A <- as.matrix.csr(0, sum(nrA), sum(ncA))
    if (sum(nrR) > 0) {
      R <- as.matrix.csr(0, sum(nrR), sum(ncA))
      nrR <- cumsum(nrR)
    }
    ncA <- cumsum(ncA)
    nrA <- cumsum(nrA)
    lambdas <- rep(0, mqss)
    for (i in 1:mqss) {
      lambdas[i] <- qss[[i]]$lambda
      Arows <- (1 + nrA[i]):nrA[i + 1]
      Acols <- (1 + ncA[i]):ncA[i + 1]
      A[Arows, Acols] <- qss[[i]]$lambda * qss[[i]]$A
      if (nrR[i] < nrR[i + 1])
        R[(1 + nrR[i]):nrR[i + 1], (1 + ncA[i]):ncA[i + 1]] <- qss[[i]]$R
    }
    A <- cbind(as.matrix.csr(0, nrA[mqss + 1], p), A)
    if (nrR[mqss + 1] > 0) {
      R <- cbind(as.matrix.csr(0, nrR[mqss + 1], p), R)
      r <- rep(0, nrR[mqss + 1])
    }
    else {
      R <- NULL
      r <- NULL
    }
    if (method == "lasso") {
      A <- rbind(cbind(L, as.matrix.csr(0, nrL, ncol(F) -
                                          ncL)), A)
    }

   if (length(weights)) {
    wF<-F*weights
    wY<-Y*weights
    X <- rbind(wF, A)
    Y <- c(wY, rep(0, nrow(A)))
    rhs <- selectMethod("t",signature = "matrix.csr")(as.matrix.csr(rbind((1 - tau) * wF, 0.5 * A))) %*% rep(1, nrow(X))
    }
    else{
      X <- rbind(F, A)
      Y <- c(Y, rep(0, nrow(A)))
      rhs <- selectMethod("t",signature = "matrix.csr")(rbind((1 - tau) * F, 0.5 * A)) %*% rep(1, nrow(X))
   }

    # if(mqss!=1) XpX <- t(X) %*% X
    # if(mqss==1){m<-ncol(X);nnzdmax<-m+9*2+6*(m-5)} else{nnzdmax <- XpX@ia[length(XpX@ia)] - 1}
    XpX <-selectMethod("t",signature = "matrix.csr")(X) %*% X
    nnzdmax <- XpX@ia[length(XpX@ia)] - 1
    if (is.null(control[["nsubmax"]]))
      control[["nsubmax"]] <- max(nnzdmax, floor(1000 + exp(-1.6) * nnzdmax^1.2))
    if (is.null(control[["nnzlmax"]]))
      control[["nnzlmax"]] <- floor(2e+05 - 2.8 * nnzdmax +  7e-04 * nnzdmax^2)
    if (is.null(control[["tmpmax"]]))
      control[["tmpmax"]] <- floor(1e+05 + exp(-12.1) * nnzdmax^2.35)


    fit <- if (length(r) > 0)
      rqss.fit(X, Y, tau = tau, rhs = rhs, method = "sfnc",
               R = R, r = r, control = control, ...)
    else rqss.fit(X, Y, tau = tau, rhs = rhs, method = "sfn",
                  control = control, ...)


    for (i in 1:mqss) {
      ML <- p + 1 + ncA[i]
      MU <- p + ncA[i + 1]
      qss[[i]] <- list(xyz = cbind(qss[[i]]$x$x, qss[[i]]$x$y,
                                   c(0, fit$coef[ML:MU])), dummies = qss[[i]]$dummies)
      if (ncol(qss[[i]]$xyz) == 2)
        class(qss[[i]]) <- "qss1"
      else class(qss[[i]]) <- "qss2"
    }
    names(qss) <- qssnames
    fit$qss <- qss
  }
  else {
    X <- as.matrix.csr(X)
    nrA <- 0
    if (method == "lasso") {
      rhs <- t(rbind((1 - tau) * X, 0.5 * L)) %*% rep(1, nrow(X) + nrow(L))
      X <- rbind(X, L)
      Y <- c(Y, rep(0, nrL))
    }
    else rhs <- NULL
    if (length(weights)) {
      if (any(weights < 0))
        stop("negative weights not allowed")
      X <- X * weights
      Y <- Y * weights
    }
    fit <- rqss.fit(X, Y, tau = tau, rhs = rhs, control = control, ...)
    fit$nrA <- nrA
  }


  names(fit$coef) <- vnames
  n <- length(fit$resid) - nrL - nrA[length(nrA)]
  uhat <- fit$resid[1:n]
  Rho <- function(u, tau) sum(u * (tau - (u < 0)))
  fit$fidelity <- Rho(uhat, tau)
  fit$edf <- sum(abs(uhat) < ztol)
  fit$X <- X
  fit$n <- n
  fit$nrL <- nrL
  fit$terms <- Terms
  fit$fake.formula <- ff
  fit$formula <- formula
  fit$method <- method
  fit$call <- call
  fit$tau <- tau
  if (length(qssterms)) {
    fit$lambdas <- lambdas
    fit$qssnames <- qssnames
    fit$nrA <- nrA
    fit$ncA <- cumsum(c(p, diff(ncA)))
  }
  else fit$ncA <- p
  attr(fit, "na.message") <- attr(m, "na.message")
  class(fit) <- "rqss"
  fit
}

