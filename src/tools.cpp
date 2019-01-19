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
NumericVector cppReglinSimplified(
    const NumericMatrix& x,
    const NumericMatrix& y,
    const IntegerVector& sizes
) {
  const size_t width = x.rows();
  const size_t height = x.cols();
  std::vector<size_t> colSizes = Rcpp::as<std::vector<size_t> >(sizes);
  NumericVector result(width);
  /*extern "C"
   void calculateReglinSimplified(
   const double *h_x,
   const double *h_y,
   const size_t width,
   const size_t height,
   const size_t *h_cols_sizes,
   double *h_results
   );*/
  calculateReglinSimplified(
    &x[0], &y[0], width, height, &colSizes[0], &result[0]
  );
  return result;
}

// [[Rcpp::export]]
NumericVector cppChisq(
    const NumericMatrix& x,
    const NumericMatrix& y, 
    const NumericMatrix& e,
    const IntegerVector& sizes
) {
  const size_t width = x.rows();
  const size_t height = x.cols();
  std::vector<size_t> sizesVector(sizes.begin(), sizes.end());
  std::vector<double> results(width);
  calculateChisq(&x[0],&y[0],&e[0],&sizesVector[0],&results[0],width,height);
  return NumericVector(results.begin(), results.end());
}