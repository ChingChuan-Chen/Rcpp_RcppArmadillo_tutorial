#define RCPP_ARMADILLO_RETURN_ANYVEC_AS_VECTOR
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// [[Rcpp::export]]
Rcpp::List armaLmFunc(const arma::mat& X, const arma::vec& y) {
  // get degree of freedom
  arma::uword dfResidual = X.n_rows - X.n_cols;
  // get pseudo-inverse of matrix XTX
  arma::mat pinvXTX = arma::pinv(arma::trans(X)*X);
  // get coef of linear model
  arma::vec coef = pinvXTX * arma::trans(X) * y;
  // get residuals
  arma::vec res  = y - X*coef;
  // calculate variance of residuals
  double errVar = dot(res, res) / dfResidual;
  // std.errors of coefficients
  arma::vec se = arma::sqrt(errVar * arma::diagvec(pinvXTX));
  // return
  return Rcpp::List::create(Rcpp::Named("coefficients") = coef,
                            Rcpp::Named("errVar")       = errVar,
                            Rcpp::Named("se")           = se,
                            Rcpp::Named("dfResidual")  = dfResidual);
}
