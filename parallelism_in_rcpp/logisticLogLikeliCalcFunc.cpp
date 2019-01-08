// [[Rcpp::depends(RcppArmadillo, RcppParallel)]]
#include <RcppArmadillo.h>
#include <RcppParallel.h>

inline arma::vec logisticFunc(arma::vec t) {
  return(arma::clamp(arma::pow(1.0 + arma::exp(-t), -1.0), 1e-15, 1-1e-15));
}

struct LogisticLogLikeliCalcWorker: public RcppParallel::Worker {
  // inputs
  const arma::mat& X;
  const arma::vec& y;
  const arma::vec& beta;
  // output
  double logisticLogLikeli;
  
  // init
  LogisticLogLikeliCalcWorker(const arma::mat& X, const arma::vec& y, const arma::vec& beta):
    X(X), y(y), beta(beta), logisticLogLikeli(0) {}
  
  LogisticLogLikeliCalcWorker(const LogisticLogLikeliCalcWorker& worker, RcppParallel::Split):
    X(worker.X), y(worker.y), beta(worker.beta), logisticLogLikeli(0) {}

  // calculation
  void operator()(std::size_t begin, std::size_t end) {
    for (arma::uword i = begin; i < end; ++i) {
      double pHat = arma::as_scalar(logisticFunc(X.row(i) * beta));
      logisticLogLikeli += y(i) * std::log(pHat) + (1.0 - y(i)) * std::log(1.0 - pHat);
    }
  }
  
  // add logisticLogLikeli together
  void join(const LogisticLogLikeliCalcWorker& rhs) {
    logisticLogLikeli += rhs.logisticLogLikeli; 
  }
};

// [[Rcpp::export]]
double armaLogisticLogLikeliCalcFunc(const arma::mat& X, const arma::vec& y, const arma::vec& beta) {
  // init worker
  LogisticLogLikeliCalcWorker logisticLogLikeliCalcWorkerInstance(X, y, beta);
  // implementation reduce
  parallelReduce(0, X.n_rows, logisticLogLikeliCalcWorkerInstance);
  // return
  return logisticLogLikeliCalcWorkerInstance.logisticLogLikeli;
}
