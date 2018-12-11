#include <Rcpp.h>
using namespace Rcpp;

#include "includes/quasar_spectrum.h"
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
  std::vector<double4> params_m;
  const size_t size = params.size();
  memcpy((void*)&params_m[0], (void*)&params[0], sizeof(double4) * params.size());
  
  for(int i = 0; i < size; ++ i) {
    double4 pr = params_m[0];
    Rcout << "a " << pr.x << " b " << pr.y << " z " << pr.z << std::endl;
  }
    
  NumericMatrix result(4096, size);
  // generateWavelengths(
  //     &result[0],
  //     params_m,
  //     params.nrow()
  // );
  return NumericMatrix(1, 1);
}

