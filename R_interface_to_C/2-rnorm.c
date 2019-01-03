#include <R.h>
#include <Rmath.h>

void C_rnorm(int *n, double* mu, double* sigma, double* x) {
  GetRNGstate();
  for (int i = 0; i < *n; ++i)
    x[i] = rnorm(*mu, *sigma);
  PutRNGstate();
}
