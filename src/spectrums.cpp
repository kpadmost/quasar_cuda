#include <Rcpp.h>
using namespace Rcpp;

#include "includes/quasar_spectrum.h"
#include "RcppExtended.hpp"
// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
NumericMatrix cppGenerateWavelenghtMatrix(SEXP params) {
  std::vector<double4> params_m = as< std::vector<double4>>(params);
  const size_t spectrum_number = params_m.size();
  NumericMatrix result(ASTRO_OBJ_SIZE, spectrum_number);
  try {
    generateWavelengths(
        &result[0],
        &params_m[0],
        spectrum_number
    );
  } catch(const std::runtime_error& e) {
    Rcpp::Rcout << e.what();
    exit(1);
  }
  return result;
}


//NOT TESTED!
// [[Rcpp:export]]
NumericMatrix cppSingleInterpolation(NumericMatrix xMatrix, NumericMatrix yMatrix, IntegerVector xSizes, 
                                     NumericMatrix sMatrix, NumericMatrix tMatrix, const size_t sSizes) {
  // copy to array
  const size_t spectrum_number = xMatrix.rows();
  const size_t spectrum_size = xMatrix.cols();
  std::vector<size_t> sizes(xSizes.size());
  std::copy(xSizes.begin(), xSizes.end(), sizes.begin());
  // call kernel 
  NumericMatrix output(spectrum_number, spectrum_size);
  singleInterpolation(
    &xMatrix[0], 
    &yMatrix[0], 
    &sizes[0], 
    spectrum_number, 
    &sMatrix[0],
    &tMatrix[0],
    spectrum_size,
    &output[0]
  );
  // enjoy
  return output;
}

//TODO:considerGrThan

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
      Named("wavelengthMatrix") = wavelengthMatrix_c,
      Named("spectrumsMatrix") = spectrumMatrix_c,
      Named("errorsMatrix") = errorMatrix_c
    );
}

