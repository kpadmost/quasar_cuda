#ifndef TOOLS_HPP
#define TOOLS_HPP
#include <Rcpp.h>

namespace quasarCuda {
  inline double* RcppNumericMatrixToCArray(Rcpp::NumericMatrix& input) {
    return &input[0];
  }

  inline double* RcppNumericVectorToCArray(Rcpp::NumericVector& input) {
    return &input[0];
  }

  inline Rcpp::NumericMatrix CArrayToRcppNumericMatrix(double* input, const size_t width, const size_t height) {
    //dataVec.insert(dataVec.end(), &dataArray[0], &dataArray[dataArraySize]);
    Rcpp::NumericMatrix result(width, height);
    for(int i = 0; i < width; ++i) 
      for(int j = 0; j < height; ++j) 
        result(j, i) = input[i * height + j];
    Rcpp::Rcout << result;
    return result;
  }
}

#endif
