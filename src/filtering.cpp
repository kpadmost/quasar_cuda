#include <Rcpp.h>
#include <algorithm>

#include "includes/filtering.h"

#include "RcppExtended.hpp"


using namespace Rcpp;



// [[Rcpp::export]]
IntegerVector cppCountNInfSize(const NumericMatrix& inputMatrix) {
  const size_t inputWidth = inputMatrix.rows(); // number of spectrums
  const size_t inputHeight = inputMatrix.cols();
  std::vector<size_t> outputSizes(inputWidth);
  // copy onto kernel
  countNotInf(&inputMatrix[0], &outputSizes[0], inputWidth, inputHeight);
  std::vector<size_t>::iterator newHeight = std::max_element(outputSizes.begin(), outputSizes.end());
  return IntegerVector(outputSizes.begin(), outputSizes.end());
  // determine size
}

// [[Rcpp::export]]
NumericMatrix cppCopyNInf(const NumericMatrix& inputMatrix, const size_t newHeight) {
  const size_t inputWidth = inputMatrix.rows(); // number of spectrums
  const size_t inputHeight = inputMatrix.cols();
  NumericMatrix result(inputWidth, newHeight);
  copyIfNotInf(&inputMatrix[0], &result[0], inputWidth, inputHeight, newHeight);
  return result;
}

// [[Rcpp::export]]
List cppFilterWithValues(
    const NumericMatrix& wavelengthMatrix, 
    const NumericMatrix& spectrumMatrix, 
    const NumericMatrix& errorMatrix, 
    const NumericVector& sizes)
{
  const size_t width = spectrumMatrix.ncol();
  Rcout << " w " << width;
  const size_t height = spectrumMatrix.nrow();
  std::vector<size_t> real_sizes(sizes.size());
  std::copy(sizes.begin(), sizes.end(), real_sizes.begin());
  
  NumericMatrix wavelengthMatrix_c(wavelengthMatrix); 
  NumericMatrix spectrumMatrix_c(spectrumMatrix); 
  NumericMatrix errorMatrix_c(errorMatrix);
  try {
    filterNonpositive(
      &spectrumMatrix_c[0],
                       &wavelengthMatrix_c[0],
                       &wavelengthMatrix_c[0],
                       &real_sizes[0],
                       height
    );
  } catch(const std::runtime_error& e) {
    Rcpp::Rcout << e.what();
    exit(1);
  }
  return List::create(
    Named("wavelengthsMatrix") = wavelengthMatrix_c,
    Named("spectrumsMatrix") = spectrumMatrix_c,
    Named("errorsMatrix") = errorMatrix_c
  );
}

// [[Rcpp::export]]
List cppFilterWithWavelengthWindows(
    const NumericMatrix& wavelengthsMatrix,
    const NumericMatrix& spectrumsMatrix,
    const NumericMatrix& errorsMatrix,
    const IntegerVector& sizesVectorR,
    SEXP continuumWindowsVectorR
    ) {
  const size_t width = wavelengthsMatrix.rows();
  const size_t height = wavelengthsMatrix.cols();
  // work on copys
  NumericMatrix wavelengthsMatrixC(wavelengthsMatrix);
  NumericMatrix spectrumsMatrixC(spectrumsMatrix);
  NumericMatrix errorsMatrixC(errorsMatrix);
  
  std::vector<double2> continuumWindowsVector = as<std::vector<double2>>(continuumWindowsVectorR);
  for(auto i = continuumWindowsVector.begin(); i != continuumWindowsVector.end(); ++i)
    Rcout << "s " << i->x << " " << i->y << " ";
  std::vector<size_t> sizesVector(sizesVectorR.begin(), sizesVectorR.end());
  filterWithWavelengthsWindows(
    &wavelengthsMatrixC[0],
    &spectrumsMatrixC[0],
    &errorsMatrixC[0],
    &sizesVector[0],
    &continuumWindowsVector[0],
    continuumWindowsVector.size(),
    width,
    height
  );
  return List::create(
    Named("wavelengthsMatrix") = wavelengthsMatrixC,
    Named("spectrumsMatrix") = spectrumsMatrixC,
    Named("errorsMatrix") = errorsMatrixC
  );
  /*void filterWithWavelengthsWindows(
   double *h_wavelengths,
   double *h_spectrums,
   double *h_errors,
   const size_t *h_sizes,
   const double2 *h_windows,
   const size_t windows_number,
   const size_t width,
   const size_t height
  );*/
}