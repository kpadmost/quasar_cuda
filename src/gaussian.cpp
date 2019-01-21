#include "RcppExtended.hpp"
using namespace Rcpp;
#include "includes/gaussian.h"


// [[Rcpp::export]]
SEXP cppFitGaussian(
  const NumericMatrix& x,
  const NumericMatrix& y,
  const IntegerVector& sizes,
  SEXP results
) {
  const size_t width = x.rows();
  const size_t height = x.cols();
  std::vector<size_t> sizesVector = Rcpp::as<std::vector<size_t> >(sizes);
  std::vector<double4> resultsVector = Rcpp::as<std::vector<double4> >(results);
  fitGaussian(
    &x[0],
    &y[0],
    &sizesVector[0],
    &resultsVector[0],
    width,
    height,
    MAX_FITGAUSSIAN_ITER
  );
  return(Rcpp::wrap(resultsVector));
}

// [[Rcpp::export]]
NumericVector cppCalculateFWHM(
  SEXP fitParams
) {
  std::vector<double4> fitParamsVector = Rcpp::as<std::vector<double4> >(fitParams);
  const size_t vectorSize = fitParamsVector.size();
  NumericVector result(vectorSize);
  calculateFwhm(
    &fitParamsVector[0],
    &result[0],
    vectorSize
  );
  return result;
}


// [[Rcpp::export]]
NumericVector cppCalculateGaussianChisq(
  const NumericMatrix& wavelengthsMatrix,
  const NumericMatrix& spectrumsMatrix,
  const NumericVector& errorsMatrix,
  SEXP fitGaussianParams,
  const IntegerVector& sizes
) {
  const uint width = wavelengthsMatrix.rows();
  const uint height = wavelengthsMatrix.cols();
  std::vector<uint> sizesVector = Rcpp::as<std::vector<uint>>(sizes);
  std::vector<double4> fitParamsVector = Rcpp::as<std::vector<double4>>(fitGaussianParams);
  NumericVector result(width);
  calculateGaussianChisq(
    &wavelengthsMatrix[0],
    &spectrumsMatrix[0],
    &errorsMatrix[0],
    &fitParamsVector[0],
    &sizesVector[0],
    &result[0],
    width,
    height
  );
  return result;
}


// [[Rcpp::export]]
NumericMatrix cppCalculateGaussian(
  const NumericMatrix& xMatrix,
  SEXP fitGaussianParams,
  const IntegerVector& sizes 
) {
  const uint width = xMatrix.nrow();
  const uint height = xMatrix.ncol();
  std::vector<uint> sizesVector = Rcpp::as<std::vector<uint>>(sizes);
  std::vector<double4> fitParamsVector = Rcpp::as<std::vector<double4>>(fitGaussianParams);
  NumericMatrix result(width, height);
  calculateGaussian(
    &xMatrix[0],
    &result[0],
    &fitParamsVector[0],
    &sizesVector[0],
    width,
    height
  );
  return result;
}