#define RCPP_ARMADILLO_RETURN_ANYVEC_AS_VECTOR
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// [[Rcpp::export]]
Rcpp::List armaLmFunc(const arma::mat& X, const arma::vec& y) {
  arma::uword dfResidual = X.n_rows - X.n_cols;
  // fit model y ~ X, extract residuals
  arma::mat pinvXTX = arma::pinv(arma::trans(X)*X);
  arma::vec coef = pinvXTX * arma::trans(X) * y;
  arma::vec res  = y - X*coef;
  double errVar = dot(res, res) / dfResidual;
  // std.errors of coefficients
  arma::vec se = arma::sqrt(errVar * arma::diagvec(pinvXTX));
  return Rcpp::List::create(Rcpp::Named("coefficients") = coef,
                            Rcpp::Named("errVar")       = errVar,
                            Rcpp::Named("se")           = se,
                            Rcpp::Named("dfResidual")  = dfResidual);
}
