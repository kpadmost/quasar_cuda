#include <Rcpp.h>

#include "includes/continuum.h"

#include "RcppExtended.hpp"
using namespace Rcpp;

List cppCalculateCfunDcfun(const NumericMatrix& wavelengthsMatrix, SEXP rReglin, SEXP rCReglin) {
  std::vector<double8> reglin = as<std::vector<double8>>(rReglin);
  std::vector<double8> cReglin = as<std::vector<double8>>(rCReglin);
  const size_t width = wavelengthsMatrix.rows();
  const size_t height = wavelengthsMatrix.cols();
  NumericMatrix continuumMatrix(width, height);
  NumericMatrix dContinuumMatrix(width, height);
  calculateDcfun(
    &wavelengthsMatrix[0],
    &dContinuumMatrix[0],
    &continuumMatrix[0],
    &reglin[0],
    &cReglin[0],
    width,
    height
  );
  return List::create(
    Named("cfun") = continuumMatrix,
    Named("dcfun") = dContinuumMatrix
  );
}

/*extern "C"
 void reduceChisqs(
 double *h_chisq,
 const size_t *h_sizes,
 const size_t width
 );*/

NumericVector cppReduceContinuumChisq(NumericVector& chisq, const IntegerVector sizes) {
  const size_t spectrumsNumber = chisq.size();
  std::vector<size_t> sizesVector(sizes.begin(), sizes.end());
  NumericVector chisqCopy(chisq.begin(), chisq.end());
  reduceChisqs(
    &chisqCopy[0],
              &sizesVector[0],
              spectrumsNumber
  );
  
  return chisqCopy;
  
}