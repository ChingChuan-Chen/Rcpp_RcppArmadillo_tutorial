#include <Rcpp.h>
// [[Rcpp::plugins(openmp)]]
#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::export]]
double long_computation(int nb, int threads = 1) {
#ifdef _OPENMP
  if (threads > 0) omp_set_num_threads(threads);
  Rprintf("Number of threads=%i\n", omp_get_max_threads());
#endif
  
  double sum = 0;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
  for (int i = 0; i < nb; ++i) {
    for (int j = 0; j < nb; ++j) {
      sum += R::dlnorm(i+j, 0.0, 1.0, 0);
    }
  }
  return sum + nb;
}
