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
NumericMatrix cppGenerateWavelenghtMatrix(SEXP params, const unsigned int size) {
  std::vector<double4> params_m = as< std::vector<double4>>(params);
  const size_t spectrum_number = params_m.size();
  NumericMatrix result(ASTRO_OBJ_SIZE, size);
  generateWavelengths(
      &result[0],
      &params_m[0],
      spectrum_number
  );
  return NumericMatrix(1, 1);
}

