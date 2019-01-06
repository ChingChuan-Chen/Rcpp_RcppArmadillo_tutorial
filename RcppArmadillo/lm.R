# install installr
if (!"installr" %in% rownames(installed.packages()))
  install.packages("installr")
library(installr)
# require2 will automatically install nonexistent packages
require2(Rcpp)
require2(RcppArmadillo)
require2(MASS)
require2(microbenchmark)

getRandomDataFunc <- function(n, p) {
  X <- matrix(rnorm(n*p), n)
  trueBeta <- as.matrix(sample(-5:5, p+1L, TRUE)) / 10
  y <- cbind(1, X) %*% trueBeta + rnorm(n)
  return(list(X = cbind(1, X), y = y))
}

# RcppArmadillo
Rcpp::sourceCpp("armaLmFunc.cpp")

# R
rLmFunc <- function(X, y) {
  dfResidual <- nrow(X) - ncol(X)
  pinvXTX <- MASS::ginv(t(X) %*% X)
  coef <- pinvXTX %*% t(X) %*% y
  res <- y - X %*% coef
  errVar <- as.vector(crossprod(res)) / dfResidual
  se <- sqrt(errVar * diag(pinvXTX))
  return(list(
    coefficients = as.vector(coef),
    errVar = errVar,
    se = se,
    dfResidual = dfResidual
  ))
}

# check that the results of two functions are equal
set.seed(100)

lm.fit(a$X, a$y)$df.residual

# check rLmFunc is equal to lm.fit
with(getRandomDataFunc(100L, 20L),
     all.equal(coef(lm.fit(X, y)), coef(rLmFunc(X, y)), 
               check.attributes = FALSE))
with(getRandomDataFunc(100L, 20L), 
     all.equal(summary(lm(y ~ 0 + X))$sigma, 
               sqrt(rLmFunc(X, y)$errVar), check.attributes = FALSE))
with(getRandomDataFunc(100L, 20L), 
     all.equal(summary(lm(y ~ 0 + X))$coefficients[ , 2], 
               rLmFunc(X, y)$se, check.attributes = FALSE))

# check armaLmFunc is equal to rLmFunc
with(getRandomDataFunc(100L, 20L), all.equal(armaLmFunc(X, y), rLmFunc(X, y)))

# benchmark performance
# n = 100 / p = 20
microbenchmark(
  arma = with(getRandomDataFunc(100L, 20L), armaLmFunc(X, y)),
  r = with(getRandomDataFunc(100L, 20L), rLmFunc(X, y)),
  r2 = with(getRandomDataFunc(100L, 20L), lm.fit(X, y)),
  times = 500L
)

# n = 10000 / p = 20
microbenchmark(
  arma = with(getRandomDataFunc(10000L, 20L), armaLmFunc(X, y)),
  r = with(getRandomDataFunc(10000L, 20L), rLmFunc(X, y)),
  r2 = with(getRandomDataFunc(10000L, 20L), lm.fit(X, y)),
  times = 100L
)

# n = 10000 / p = 200
microbenchmark(
  arma = with(getRandomDataFunc(10000L, 200L), armaLmFunc(X, y)),
  r = with(getRandomDataFunc(10000L, 200L), rLmFunc(X, y)),
  r2 = with(getRandomDataFunc(10000L, 200L), lm.fit(X, y)),
  times = 3L
)
