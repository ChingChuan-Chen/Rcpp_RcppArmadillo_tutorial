#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector cpp_rnorm(int n, double mu = 0.0, double sigma = 1.0) {
  NumericVector x = rnorm(n, mu, sigma);
  return x;
}

/*** R
cpp_rnorm(5)
cpp_rnorm(5, 1, 2)
*/
