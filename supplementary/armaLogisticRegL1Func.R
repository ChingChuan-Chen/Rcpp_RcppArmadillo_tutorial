# install installr
if (!"installr" %in% rownames(installed.packages()))
  install.packages("installr")
library(installr)
# require2 will automatically install nonexistent packages
require2(Rcpp)
require2(RcppArmadillo)
require2(RcppParallel)
require2(data.table)
require2(pipeR)
require2(ggplot2)

logisticFunc <- function(t) pmin(pmax(1/(1+exp(-t)), 10*.Machine$double.eps), 1 - 10*.Machine$double.eps)

getRandomDataFunc <- function(n, p) {
  X <- scale(matrix(rnorm(n*p), n))
  trueBeta <- as.matrix(sample(-5:5, p+1L, TRUE)) / 10
  y <- sapply(logisticFunc(cbind(1, X) %*% trueBeta), function(prob) rbinom(1, 1, prob))
  return(list(X = cbind(1, X), y = y))
}

# cpp file
Rcpp::sourceCpp("armaLogisticRegL1Func.cpp")
armaLogisticRegL1Func <- function(X, y, nlambda = 100L, lambda = NULL, 
                                  lambdaMinRatio = ifelse(nrow(X) < ncol(X), 1e-2, 1e-4),
                                  intercept = TRUE, epsilon = 1e-8, maxit = 1e5L) {
  if (is.null(lambda)) {
    if (intercept) {
      idx <- 2L:ncol(X)
    } else {
      idx <- 1L:ncol(X)
    }
    logLabmdaMax <- log(max(abs(t(X[ , idx]) %*% y) / nrow(X)))
    lambda <- exp(seq(logLabmdaMax, logLabmdaMax + log(lambdaMinRatio), length.out = nlambda))
  } else {
    stopifnot(all(lambda >= 0))
    nlambda <- length(lambda)
  }
  fit <- armaLogisticRegL1FuncCpp(X, y, nlambda, lambda, lambdaMinRatio, intercept, epsilon, maxit)
  # colnames(fit$results) <- c("Lambda", "%Dev", "df")
  return(fit)
}

# R
rLogisticL1RegFunc <- function(X, y, nlambda = 100L, lambda = NULL, 
                               lambdaMinRatio = ifelse(nrow(X) < ncol(X), 1e-2, 1e-4),
                               intercept = TRUE, epsilon = 1e-8, maxit = 1e5L) {
  if (is.null(lambda)) {
    if (intercept) {
      idx <- 2L:ncol(X)
    } else {
      idx <- 1L:ncol(X)
    }
    logLabmdaMax <- log(max(abs(t(X[ , idx]) %*% y) / nrow(X)))
    lambda <- exp(seq(logLabmdaMax, logLabmdaMax + log(lambdaMinRatio), length.out = nlambda))
  } else {
    stopifnot(all(lambda >= 0))
    nlambda <- length(lambda)
  }
  
  # get degree of freedom
  dfResidual <- nrow(X) - ncol(X)
  # initialize parameters
  it <- 1L
  converged <- TRUE
  lastConverged <- FALSE
  # initialize coefs
  coefs <- rep(0, ncol(X))
  startIdx <- 1L
  if (intercept) {
    coefs[1] <- log(mean(y) / (1-mean(y)))
    startIdx <- 2L
  }
  coefMat <- matrix(NA_real_, ncol(X), nlambda)
  coefMat[ , 1] <- coefs
  # get deviance and objective value
  pHat <- logisticFunc(X %*% coefs)
  deviance <- sum(-2*c(y*log(pHat), (1-y)*log(1-pHat)))
  deviances <- vector("double", nlambda)
  deviances[1] <- deviance
  nulldev <- deviance
  devianceRatio <- vector("double", nlambda)
  objOld <- deviance / 2
  # other parameters
  xColMeans <- colMeans(X)
  XWX <- 0.25 * colMeans(X^2)
  residuals <- (y - pHat) / 0.25

  # start pathwise coordiante descent
  for (l in 2L:nlambda) {
    repeat {
      # active set
      locUpdate <- seq.int(startIdx, ncol(X))
      if (!converged)
        locUpdate <- which(abs(coefs[startIdx:length(coefs)]) > 0) + (startIdx - 1L)
      
      # update intercept
      if (intercept)
        coefs[1] <- coefs[1] + crossprod(X[ , 1], residuals) / nrow(X)
      
      # pathwise coordinate descent
      for (k in seq_along(locUpdate)) {
        uv <- XWX[locUpdate[k]] * coefs[locUpdate[k]] + 0.25 * crossprod(X[ , locUpdate[k]], residuals) / nrow(X)
        uv <- ifelse(abs(uv) < lambda[l], 0, uv)
        if (uv != 0.0) {
          coefs[locUpdate[k]] <- (uv - sign(uv) * lambda[l]) / (XWX[locUpdate[k]] + lambda[l])
          pHat <- logisticFunc(X %*% coefs)
          residuals <- (y - pHat) / 0.25
        }
      }
      
      # update deviance
      deviance <- sum(-2*c(y*log(pHat), (1-y)*log((1-pHat))))
      # update objective value
      objNew <- deviance / 2 + lambda[l] * sum(abs(coefs[startIdx:length(coefs)]))
      
      # update converging status
      lastConverged <- converged
      # check whether converged
      converged <- (objOld - objNew) / (objNew - 0.1) < epsilon
      if (converged && lastConverged)
        break
      # update parameters for next iteration
      it <- it + 1L
      objOld <- objNew
    }
    
    coefMat[ , l] <- coefs
    deviances[l] <- deviance
    devianceRatio[l] <- 1 - deviance / nulldev 
    if (((devianceRatio[l] - devianceRatio[l-1]) / devianceRatio[l] < 1e-5) || (devianceRatio[l] > 0.999))
      break
    if (it > maxit)
      stop("Exceed maximum iteration")
  }
  
  # retrun
  return(list(
    results = cbind(Lambda = lambda[1:l], "%Dev" = c(NA_real_, -diff(deviances[1:l]) / deviances[2:l]),
                    df = colSums(abs(coefMat[ , 1:l]) > 0) - 1),
    coefficients = coefMat[ , 1:l],
    nulldev = nulldev,
    deviances = deviances[1:l],
    dfResidual = dfResidual,
    iter = it,
    converged = converged
  ))
}

# check correctness
set.seed(100)
zz <- getRandomDataFunc(100L, 20L)
fit <- with(zz, rLogisticL1RegFunc(X, y))

fitRes <- data.table(
  lambda = rep(fit$results[ , 1], each = nrow(fit$coefficients) - 1L), 
  variable = paste0("V", 1L:(nrow(fit$coefficients)-1L)) %>>%
    factor(levels = .),
  coef = as.vector(fit$coefficients[-1L, ])
)

ggplot(fitRes, aes(lambda, coef, colour = variable)) + geom_point() + geom_line() +
  guides(colour = guide_legend(ncol=2))


fit2 <- with(zz, armaLogisticRegL1Func(X, y))
fitRes <- data.table(
  lambda = rep(fit2$results, each = nrow(fit2$coefficients) - 1L), 
  variable = paste0("V", 1L:(nrow(fit2$coefficients)-1L)) %>>%
    factor(levels = .),
  coef = as.vector(fit2$coefficients[-1L, ])
)

ggplot(fitRes, aes(lambda, coef, colour = variable)) + geom_point() + geom_line() +
  guides(colour = guide_legend(ncol=2))
