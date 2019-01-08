# install installr
if (!"installr" %in% rownames(installed.packages()))
  install.packages("installr")
library(installr)
# require2 will automatically install nonexistent packages
require2(Rcpp)
require2(RcppArmadillo)
require2(RcppParallel)
require2(glmnet)

logisticFunc <- function(t) pmin(pmax(1/(1+exp(-t)), 10*.Machine$double.eps), 1 - 10*.Machine$double.eps)

getRandomDataFunc <- function(n, p) {
  X <- matrix(rnorm(n*p), n)
  trueBeta <- as.matrix(sample(-5:5, p+1L, TRUE)) / 10
  y <- sapply(logisticFunc(cbind(1, X) %*% trueBeta), function(prob) rbinom(1, 1, prob))
  return(list(X = cbind(1, X), y = y))
}

# cpp file
Rcpp::sourceCpp("armaLogisticRegL1Func.cpp")

# R
rLogisticRegFunc <- function(X, y, epsilon = 1e-8, maxit = 25L) {
  # get degree of freedom
  dfResidual <- nrow(X) - ncol(X)
  # initialize parameters
  it <- 0L
  converged <- TRUE
  coefOld <- rep(0, ncol(X))
  pHatOld <- y / 2 + 0.25
  devOld <- Inf
  # iterations to update coefficients
  repeat {
    # update coefficients
    coefNew <- coefOld + solve(t(X) %*% diag(as.vector(pHatOld * (1-pHatOld))) %*% X, t(X) %*% (y - pHatOld))
    # calculate new p hat
    pHatNew <- logisticFunc(X %*% coefNew)
    # calculate new deviance
    devNew <- sum(-2*c(y*log(pHatNew), (1-y)*log((1-pHatNew))))
    # check whether converged
    if (abs(devNew - devOld) / (abs(devNew) + 0.1) < epsilon)
      break
    # check whether reaches maximum iteration
    if (it >= maxit) {
      warning("Exceed maximum iteration, it is not converged!")
      converged <- FALSE
      break
    }
    # update parameters for next iteration
    it <- it + 1L
    coefOld <- coefNew
    pHatOld <- pHatNew
    devOld <- devNew
  }
  # retrun
  return(list(
    coefficients = as.vector(coefNew),
    se = sqrt(diag(solve(t(X) %*% diag(as.vector(pHatNew * (1-pHatNew))) %*% X))),
    deviance = devNew,
    dfResidual = dfResidual,
    iter = it,
    converged = converged
  ))
}
