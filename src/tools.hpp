#ifndef TOOLS_HPP
#define TOOLS_HPP
#include <Rcpp.h>

namespace quasarCuda {

  inline double* RcppNumericMatrixToCArray(Rcpp::NumericMatrix& input) {
    return &input[0];
  }

  inline double* RcppNumericMatrixToCArray(Rcpp::NumericMatrix& input, bool isN) {
    const uint width = input.cols();
    const uint height = input.cols();
    double *result = new double[width * height];
    for(int i = 0; i < width; ++i)
      for(int j = 0; j < height; ++j)
        result[i * height + j] = input(j, i);
    return result;
  }

  inline double* RcppNumericVectorToCArray(Rcpp::NumericVector& input) {
    return &input[0];
  }

  inline Rcpp::NumericMatrix CArrayToRcppNumericMatrix(double* input, const size_t width, const size_t height) {
    //dataVec.insert(dataVec.end(), &dataArray[0], &dataArray[dataArraySize]);
    Rcpp::NumericMatrix result(width, height);
    for(int i = 0; i < width; ++i) //cols
      for(int j = 0; j < height; ++j) //rows
        result(i, j) = input[i * height + j];
    return result;
  }
}

#endif
