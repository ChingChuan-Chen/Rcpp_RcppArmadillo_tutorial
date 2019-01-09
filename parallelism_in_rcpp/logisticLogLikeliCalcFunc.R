# install installr
if (!"installr" %in% rownames(installed.packages()))
  install.packages("installr")
library(installr)
# require2 will automatically install nonexistent packages
require2(Rcpp)
require2(RcppArmadillo)
require2(RcppParallel)
require2(microbenchmark)

logisticFunc <- function(t) pmin(pmax(1/(1+exp(-t)), 1e-15), 1 - 1e-15)

getRandomDataFunc <- function(n, p) {
  X <- matrix(rnorm(n*p), n)
  trueBeta <- as.matrix(sample(-5:5, p+1L, TRUE)) / 10
  y <- sapply(logisticFunc(cbind(1, X) %*% trueBeta), function(prob) rbinom(1, 1, prob))
  return(list(X = cbind(1, X), y = y))
}

# cpp file
Rcpp::sourceCpp("logisticLogLikeliCalcFunc.cpp")

# R
rLogisticLogLikeliCalcFunc <- function(X, y, beta) {
  pHat <- logisticFunc(X %*% beta)
  return(sum(c(y*log(pHat), (1-y)*log((1-pHat)))))
}

# check whether resutls are correct
set.seed(100)
randomData <- getRandomDataFunc(100L, 20L)
fit <- glm(randomData$y ~ 0 + randomData$X, family = binomial())
all.equal(logLik(fit)[1], rLogisticLogLikeliCalcFunc(randomData$X, randomData$y, coef(fit)))
all.equal(logLik(fit)[1], armaLogisticLogLikeliCalcFunc(randomData$X, randomData$y, coef(fit)))

# benchmark
randomData <- getRandomDataFunc(5000000L, 20L)
fit <- glm(randomData$y~0+randomData$X, family = binomial())
microbenchmark(
  arma = with(randomData, armaLogisticLogLikeliCalcFunc(randomData$X, randomData$y, coef(fit))),
  r = with(randomData, rLogisticLogLikeliCalcFunc(randomData$X, randomData$y, coef(fit))), 
  r2 = with(randomData, logLik(fit)),
  times = 5L
)
