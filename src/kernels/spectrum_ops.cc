#include <Rcpp.h>
using namespace Rcpp;

#include "includes/spectrum_ops.h"
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
NumericMatrix cppGenerateWavelenghtMatrix(NumericMatrix params) {
  double4* params_m;
  memcpy(params_m, &params[0], sizeof(double4) * params.nrow());
  NumericMatrix result(4096, params.nrow());
  generateWavelengths(
      &result[0],
      params_m,
      params.nrow()
  );
  return NumericMatrix(1, 1);
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
timesTwo(42)
*/
