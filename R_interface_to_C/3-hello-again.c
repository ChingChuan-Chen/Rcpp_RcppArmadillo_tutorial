#include <R.h>
#include <Rinternals.h>

SEXP hello_again(SEXP n_) {
  int n = asInteger(n_), i;
  for (i = 0; i < n; ++i)
    Rprintf("Hi ya'll %d times!\n", i);
  return R_NilValue;
}
