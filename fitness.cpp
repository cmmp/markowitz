#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double fitnessCpp(const NumericVector& w, const double lambda,
                  const NumericVector& mean_returns, const NumericMatrix& covs) {
  double riskFactor = 0;
  double returnFactor = 0;

  int n = w.size();

  for(int i = 0; i < n; i++) {
    returnFactor += w[i] * mean_returns[i];
  }

  returnFactor *= -(1 - lambda);

  for(int i = 0; i < (n - 1); i++) {
    for(int j = i + 1; j < n; j++) {
      riskFactor += w[i] * w[j] * covs(i, j);
    }
  }

  riskFactor *= lambda;

  return - (riskFactor + returnFactor);
}
