# install installr
if (!"installr" %in% rownames(installed.packages()))
  install.packages("installr")
library(installr)
# require2 will automatically install nonexistent packages
require2(Rcpp)
require2(RcppArmadillo)
require2(RcppParallel)
require2(kernlab)
require2(microbenchmark)

set.seed(100)
N <- 3000L
p <- 100L
b <- 500L
X <- matrix(rnorm(N*p), ncol = p)
center <- X[sample.int(N, b),]

# cpp
sourceCpp("rbfKernelMatrixFunc.cpp")

# R
rRbfKernelMatrixFunc <- function(X, Y, sigma = 1) kernelMatrix(rbfdot(sigma=1/(2*sigma^2)), X, Y)@.Data

# check whether results are matched
all.equal(rRbfKernelMatrixFunc(X, center), armaRbfKernelMatrixFunc(X, center)) # TRUE
all.equal(rRbfKernelMatrixFunc(X, center, 5), armaRbfKernelMatrixFunc(X, center, 5)) # TRUE

# benchmark
microbenchmark(
  kernlab = rRbfKernelMatrixFunc(X, center), 
  armaParallel = armaRbfKernelMatrixFunc(X, center),
  times = 20L
)
