#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List modified_df(DataFrame df) {
  IntegerVector a = df["a"];
  DataFrame df2 = df;
  df2["c"] = a * 2;
  return List::create(_["original"] = df, 
                      _["modified"] = as<DataFrame>(df2));
}

/*** R
modified_df(data.frame(a = 1:5, b = rnorm(5)))
*/
  