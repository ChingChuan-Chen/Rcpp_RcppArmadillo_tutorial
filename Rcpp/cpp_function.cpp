#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector cpp_function(NumericVector x, Function f) {
  NumericVector res = f(x);
  return res;
}

/*** R
cpp_function(-1:5, abs)
cpp_function(-1:5, function(x) x^2)
*/
