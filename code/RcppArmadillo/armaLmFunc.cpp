#define RCPP_ARMADILLO_RETURN_ANYVEC_AS_VECTOR
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// [[Rcpp::export]]
Rcpp::List armaLmFunc(const arma::mat& X, const arma::vec& y) {
  arma::uword df = X.n_rows - X.n_cols;
  // fit model y ~ X, extract residuals
  arma::mat pinvXTX = arma::pinv(arma::trans(X)*X);
  arma::vec coef = pinvXTX * arma::trans(X) * y;
  arma::vec res  = y - X*coef;
  double s2 = norm(res, 2) / df;
  // std.errors of coefficients
  arma::vec se = arma::sqrt(s2 * arma::diagvec(pinvXTX));
  return Rcpp::List::create(Rcpp::Named("coefficients") = coef,
                            Rcpp::Named("se")           = se,
                            Rcpp::Named("df")           = df);
}
