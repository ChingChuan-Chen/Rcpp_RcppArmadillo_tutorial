#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector multiply(NumericVector x, NumericVector y) {
  return x * y;
}

/*** R
multiply(0:4, 2:6)
*/
