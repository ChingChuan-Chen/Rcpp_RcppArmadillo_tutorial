#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector rnorm_with_opts_cpp(int n, Environment env) {
  double mu = 0.0, sigma = 1.0;
  if (env.exists("mu"))
    mu = env["mu"];
  if (env.exists("sigma"))
    sigma = env["mu"];
  NumericVector res = rnorm(n, mu, sigma);
  return res;
}

/*** R
rnorm_with_opts <- function(n, ...) {
  rnorm_with_opts_cpp(n, as.environment(list(...)))
}
rnorm_with_opts(5)
rnorm_with_opts(5, mu = 1)
rnorm_with_opts(5, mu = 1, sigma = 2)
rnorm_with_opts(5, whatever = 1) # no impact
*/
  