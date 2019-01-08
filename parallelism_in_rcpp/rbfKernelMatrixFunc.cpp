// [[Rcpp::depends(RcppArmadillo, RcppParallel)]]
#include <RcppArmadillo.h>
#include <RcppParallel.h>

// RcppParall Worker
struct RbfKernelCalcWorker: public RcppParallel::Worker {
  // input
  const arma::vec& xSqNormMat;
  const arma::rowvec& centerSqNormMat;
  const double& s2;
  // output
  arma::mat& rbfKernelMat;
  
  // init
  RbfKernelCalcWorker(const arma::vec& xSqNormMat, const arma::rowvec& centerSqNormMat, 
                      const double& s2, arma::mat& rbfKernelMat):
    xSqNormMat(xSqNormMat), centerSqNormMat(centerSqNormMat), 
    s2(s2), rbfKernelMat(rbfKernelMat) {}
  
  // calculation
  void operator()(std::size_t begin, std::size_t end) {
    for (arma::uword k = begin; k < end; ++k) {
      // get location
      arma::uword j = k / rbfKernelMat.n_rows;
      arma::uword i = k - j * rbfKernelMat.n_rows;
      // minus squared norm of x and center
      rbfKernelMat(i, j) -= (xSqNormMat(i) + centerSqNormMat(j));
      rbfKernelMat(i, j) = std::exp(rbfKernelMat(i, j) / s2);
    }
  }
};

// [[Rcpp::export]]
arma::mat armaRbfKernelMatrixFunc(arma::mat x, arma::mat center, double sigma = 1.0) {
  // init return matrix with X^T * C 
  arma::mat rbfKernelMat(x * center.t());
  // get 
  arma::vec xSqNormMat = sum(square(x), 1) / 2.0;
  // get C^T * C / 2
  arma::rowvec centerSqNormMat = (sum(square(center), 1)).t() / 2.0;
  // get sigma^2
  double s2 = pow(sigma, 2.0);
  
  // init worker
  RbfKernelCalcWorker rbfKernelCalcWorkerInstance(xSqNormMat, centerSqNormMat, s2, rbfKernelMat);
  // start parallel computing
  RcppParallel::parallelFor(0, rbfKernelMat.n_elem, rbfKernelCalcWorkerInstance);
  // return
  return rbfKernelMat;
}
