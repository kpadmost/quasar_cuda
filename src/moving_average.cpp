#include <Rcpp.h>
using namespace Rcpp;
#include "includes/moving_average.h"
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
NumericMatrix cppMovingAverage(const NumericMatrix& inputMatrix, const IntegerVector& colsSizeVector, uint windowSize) {
  const size_t width = inputMatrix.rows();
  const size_t height = inputMatrix.cols();
  std::vector<uint> colsSizes(colsSizeVector.size());
  std::copy(colsSizeVector.begin(), colsSizeVector.end(), colsSizes.begin());
  NumericMatrix output(width, height);
  movingAverage(&inputMatrix[0], &output[0], &colsSizes[0], width, height, windowSize);
  return output;
}

