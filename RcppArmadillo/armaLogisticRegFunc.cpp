#define RCPP_ARMADILLO_RETURN_ANYVEC_AS_VECTOR
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

inline arma::vec logisticFunc(arma::vec t) {
  return(arma::clamp(arma::pow(1.0 + arma::exp(-t), -1.0), 1e-15, 1-1e-15));
}

// [[Rcpp::export]]
Rcpp::List armaLogisticRegFunc(const arma::mat& X, const arma::vec& y, double epsilon = 1e-8, arma::uword maxit = 25) {
  // get degree of freedom
  arma::uword dfResidual = X.n_rows - X.n_cols;
  // initialize parameters
  arma::uword it = 0;
  bool converged = true;
  arma::vec coefOld = arma::zeros<arma::vec>(X.n_cols), coefNew,
    pHatOld = y / 2.0 + 0.25, pHatNew;
  double devOld = 9999999.0, devNew;
  // iterations to update coefficients
  while (true) {
    // update coefficients
    coefNew = coefOld + arma::solve(arma::trans(X) * arma::diagmat(pHatOld % (1.0 - pHatOld)) * X, arma::trans(X) * (y - pHatOld));
    /* For performance concern, following code is prefered
     * 
     * arma::mat X2 = X;
     * arma::vec w = arma::sqrt(pHatOld % (1.0 - pHatOld));
     * X2.each_col() %= w
     * coefNew = coefOld + arma::solve(arma::trans(X2) * X2, arma::trans(X) * (y - pHatOld));
     *
     */
    // calculate new p hat
    pHatNew = logisticFunc(X * coefNew);
    // calculate new deviance
    devNew = -2.0 * (sum(y % arma::log(pHatNew)) + sum((1.0 - y) % arma::log(1.0-pHatNew)));
    // check whether converged
    if (std::abs(devNew - devOld) / (std::abs(devNew) + 0.1) < epsilon)
      break;
    // check whether reaches maximum iteration
    if (it >= maxit) {
      Rcpp::warning("Exceed maximum iteration, it is not converged!");
      converged = false;
      break;
    }
    // update parameters for next iteration
    ++it;
    coefOld = coefNew;
    pHatOld = pHatNew;
    devOld = devNew;
  }
  
  // std.errors of coefficients
  arma::vec se = arma::sqrt(diagvec((arma::trans(X) * arma::diagmat(pHatOld % (1.0 - pHatOld)) * X).i()));
  /* For performance concern, following code is prefered
   * 
   * arma::mat X2 = X;
   * arma::vec w = arma::sqrt(pHatOld % (1.0 - pHatOld));
   * X2.each_col() %= w
   * arma::vec se = arma::sqrt(diagvec((arma::trans(X2) * X2).i()));
   *
   */
  // return
  return Rcpp::List::create(Rcpp::Named("coefficients") = coefNew,
                            Rcpp::Named("se")           = se,
                            Rcpp::Named("deviance")     = devNew,
                            Rcpp::Named("dfResidual")   = dfResidual,
                            Rcpp::Named("iter")         = it,
                            Rcpp::Named("converged")    = converged);
}
