#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
void hello(int n = 1) {
  for (int i = 0; i < n; ++i)
    Rprintf("Hi ya'll %d times!\n", i);
}

/*** R
hello()
hello(5)
*/
