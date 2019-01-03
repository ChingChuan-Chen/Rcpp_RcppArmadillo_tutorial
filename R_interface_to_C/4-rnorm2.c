#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>

SEXP rnorm2(SEXP n_, SEXP mu_, SEXP sigma_) {
  int n = asInteger(n_);
  double mu = asReal(mu_), sigma = asReal(sigma_);

  SEXP result = PROTECT(allocVector(REALSXP, n));
  int i;
  for (i = 0; i < n; ++i)
    REAL(result)[i] = rnorm(mu, sigma);
  UNPROTECT(1);
  return result;
}
