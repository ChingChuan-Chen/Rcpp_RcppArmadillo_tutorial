# install installr
if (!"installr" %in% rownames(installed.packages()))
  install.packages("installr")
library(installr)
# require2 will automatically install nonexistent packages
require2(Rcpp)
require2(RcppArmadillo)
require2(MASS)
require2(microbenchmark)

logisticFunc <- function(t) 1/(1+exp(-t))

getRandomDataFunc <- function(n, p) {
  X <- matrix(rnorm(n*p), n)
  trueBeta <- as.matrix(sample(-5:5, p+1L, TRUE)) / 10
  y <- sapply(logisticFunc(cbind(1, X) %*% trueBeta), function(prob) rbinom(1, 1, prob))
  return(list(X = cbind(1, X), y = y))
}

# RcppArmadillo
Rcpp::sourceCpp("armaLogisticRegFunc.cpp")

# R
rLogisticRegFunc <- function(X, y, tol = 1e-8, maxit = 25L) {
  df <- nrow(X) - ncol(X)
  yhatInit <- y / 2 + 0.25
  coef_old <- solve(t(X) %*% diag(yhatInit * (1-yhatInit)) %*% X, t(X) %*% (y - yhatInit))
  it <- 1L
  converged <- TRUE
  repeat {
    yhat <- pmax(logisticFunc(X %*% coef_old, tol))
    W <- diag(as.vector(yhat * (1-yhat)))
    coef_new <- coef_old + solve(t(X) %*% W %*% X, t(X) %*% (y - yhat))
    if (all(abs(coef_new - coef_old) < tol))
      break
    if (it > maxit) {
      warning("Exceed maximum iteration, it is not converged!")
      converged <- FALSE
      break
    }
    coef_old <- coef_new
    it <- it + 1L
  }
  
  
  
  return(list(
    coefficients = as.vector(coef_old),
    se = sqrt(diag(solve(t(X) %*% W %*% X))),
    deviance = 0,
    null.deviance = 0,
    df = df,
    iter = it,
    converged = converged
  ))
}

# check that the results of two functions are equal
set.seed(100)

a <- getRandomDataFunc(100L, 20L)
with(a, coef(glm.fit(X, y, family = binomial())))
zz <- with(a, rLogisticRegFunc(X, y))

yhat <- logisticFunc(a$X %*% zz$coefficients)
sum(sqrt(c(-2*log(1-yhat[a$y == 0]), -2*log(yhat[a$y == 1]))))

sqrt(-2*log(1-yhat[2]))
sqrt(-2*log(yhat[2]))

cbind(zz$coefficients, zz$se)
summary(glm(a$y ~ 0 + a$X, family = binomial()))$coefficients[ , 1:2]

summary(glm(a$y ~ 0 + a$X, family = binomial()))

glm.fit(a$X, a$y)$res[2]

aa <- binomial()$linkfun((a$y + 0.5)/(1 + 1))
b <- binomial()$linkinv(aa)
binomial()$dev.resids
binomial()$mu.eta(0)
binomial()$linkinv(1)

residuals <- (y - mu)/mu.eta(eta)


log(0.25/0.75)

diag(summary(glm(a$y ~ 0 + a$X, family = binomial()))$cov.scaled)

# se
summary(glm(a$y ~ 0 + a$X, family = binomial()))$coefficients[ , 2]

with(a, coef(glm.fit(X, y, family = binomial())))

with(getRandomDataFunc(100L, 20L), 
     all.equal(armaLogisticRegFunc(X, y), rLogisticRegFunc(X, y)))

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




