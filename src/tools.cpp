#include "RcppExtended.hpp"

#include "includes/tools.h"

using namespace Rcpp;

// [[Rcpp::export]]
SEXP cppReglin(const NumericMatrix& x, const NumericMatrix& y, const IntegerVector& sizes) {
  const size_t width = x.rows();
  const size_t height = x.cols();
  std::vector<double8> reglinResults(width);
  std::vector<size_t> sizesVector(sizes.begin(), sizes.end());
  calculateReglin(
    &x[0],
    &y[0],
    &sizesVector[0],
    &reglinResults[0],
    width,
    height
  );
  
  return Rcpp::wrap(reglinResults);
}

// [[Rcpp::export]]
NumericVector cppChisq(
    const NumericMatrix& x,
    const NumericMatrix& y, 
    const NumericMatrix& e,
    const NumericVector& sizes
) {
  const size_t width = x.rows();
  const size_t height = x.cols();
  std::vector<size_t> sizesVector(sizes.begin(), sizes.end());
  std::vector<double> results(width);
  calculateChisq(&x[0],&y[0],&e[0],&sizesVector[0],&results[0],width,height);
  return NumericVector(results.begin(), results.end());
}