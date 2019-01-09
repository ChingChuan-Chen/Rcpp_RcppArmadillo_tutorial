#define RCPP_ARMADILLO_RETURN_ANYVEC_AS_VECTOR
// [[Rcpp::depends(RcppArmadillo, RcppParallel)]]
#include <RcppArmadillo.h>
#include <RcppParallel.h>

inline arma::vec logisticFunc(arma::vec t) {
  return(arma::clamp(arma::pow(1.0 + arma::exp(-t), -1.0), 1e-15, 1-1e-15));
}

// [[Rcpp::export]]
Rcpp::List armaLogisticRegL1FuncCpp(const arma::mat& X, const arma::vec& y, 
                                    const arma::uword& nlambda, const arma::vec& lambda, 
                                    const double& lambdaMinRatio, const bool& intercept, 
                                    const double& epsilon, arma::uword& maxit) {
  // get degree of freedom
  arma::uword dfResidual = X.n_rows - X.n_cols;
  // initialize parameters
  arma::uword it = 1;
  bool converged = true, lastConverged = false;
  // initialize coefs
  arma::vec coefs = arma::zeros<arma::vec>(X.n_cols);
  arma::uword startIdx = 0L;
  if (intercept) {
    coefs(0) = std::log(arma::mean(y) / (1 - arma::mean(y)));
    ++startIdx;
  }
  arma::mat coefMat = arma::zeros<arma::mat>(X.n_cols, nlambda);
  coefMat.col(0) = coefs;
  // get deviance and objective value
  arma::vec pHat = logisticFunc(X * coefs), deviances = arma::zeros<arma::vec>(nlambda), 
    devianceRatio = arma::zeros<arma::vec>(nlambda);
  double deviance = -2.0 * (sum(y % arma::log(pHat)) + sum((1.0 - y) % arma::log(1.0-pHat))), 
    nulldev = deviance, objOld = deviance / 2.0, objNew;
  deviances(0) = deviance;
  // other parameters
  arma::rowvec xColMeans = arma::mean(X, 0), XWX = 0.25 * arma::mean(arma::square(X), 0);
  arma::vec residuals = (y - pHat) / 0.25;
  arma::uword l;
  arma::uvec locUpdate;
  double uv, sign;
  
  for (l = 1; l < nlambda; ++l) {
    while (true) {
      // active set
      if (converged) {
        locUpdate = arma::linspace<arma::uvec>(startIdx, coefs.n_elem-1, coefs.n_elem);
      } else {
        locUpdate = arma::find(arma::abs(coefs) > 0.0);
      }
      
      // update intercept
      if (intercept)
        coefs(0) = coefs(0) + arma::dot(X.col(0), residuals) / (double) X.n_rows;
      // pathwise coordinate descent
      for (arma::uword k = 0; k < locUpdate.n_elem; ++k) {
        uv = XWX(locUpdate(k)) * coefs(locUpdate(k)) + 0.25 * dot(X.col(locUpdate(k)), residuals) / (double) X.n_rows;
        uv = (std::abs(uv) > lambda(l)) ? 0 : uv;
        if (uv != 0.0) {
          sign = (uv > 0) ? 1.0 : -1.0;
          coefs(locUpdate(k)) = (uv - sign * lambda(l)) / (XWX(locUpdate(k)) + lambda(l));
          pHat = logisticFunc(X * coefs);
          residuals = (y - pHat) / 0.25;
        }
      }
      
      // update deviance
      deviance = -2.0 * (sum(y % arma::log(pHat)) + sum((1.0 - y) % arma::log(1.0-pHat)));;
      // update objective value
      objNew = deviance / 2 + lambda[l] * sum(abs(coefs.rows(startIdx, coefs.n_elem-1)));
      
      // update converging status
      lastConverged = converged;
      // check whether converged
      converged = (objOld - objNew) / (objNew - 0.1) < epsilon;
      if (converged && lastConverged)
        break;
      // update parameters for next iteration
      ++it;
      objOld = objNew;
    }
    coefMat.col(l) = coefs;
    deviances(l) = deviance;
    devianceRatio(l) = 1 - deviance / nulldev;
    if (((devianceRatio(l) - devianceRatio(l-1)) / devianceRatio(l) < 1e-5) || (devianceRatio(l) > 0.999))
      break;
    if (it > maxit)
      Rcpp::stop("Exceed maximum iteration");
  }
  
  // return
  return Rcpp::List::create(
    Rcpp::Named("results") = lambda.rows(0, l),
    Rcpp::Named("coefficients") = coefMat,
    Rcpp::Named("nulldev")    = nulldev,
    Rcpp::Named("deviances")    = deviances,
    Rcpp::Named("dfResidual")   = dfResidual,
    Rcpp::Named("iter")         = it,
    Rcpp::Named("converged")    = converged
  );
}
