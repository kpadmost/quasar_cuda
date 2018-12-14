#include <Rcpp.h>
#include <iterator>
#include <array>
#include <algorithm>
#include "tools.hpp"
#include "includes/matrix.h"

// [[Rcpp::export]]
SEXP cppMatrixLog10(Rcpp::NumericMatrix inputMatrix) {
  const size_t width = inputMatrix.cols();   //długość widma
  const size_t height = inputMatrix.rows();  //liczba kwazarów
  const size_t N = width * height;
  double* h_input = quasarCuda::RcppNumericMatrixToCArray(inputMatrix);
  matrixLog10(h_input, height, width);
  // todo: add conversion function
  // copy 
  // cl::Buffer bufferInput = cl::Buffer(context, CL_MEM_READ_WRITE, N * sizeof(cl_double));
  // cl::copy(queue, inputMatrix.begin(), inputMatrix.end(), bufferInput);
  // 
  // quasarcl::log10(*quasarclPtr, bufferInput, height, width);
  // 
  return inputMatrix;
}

// [[Rcpp::export]]
SEXP cppMatrixAddScalar(Rcpp::NumericMatrix inputMatrix, double scalar) {
  const size_t width = inputMatrix.cols();   //długość widma
  const size_t height = inputMatrix.rows();  //liczba kwazarów
  
  double* h_input = quasarCuda::RcppNumericMatrixToCArray(inputMatrix);
  matrixAddScalar(h_input, width, height, scalar);
  // todo: add conversion function
  return inputMatrix;
}

// [[Rcpp::export]]
SEXP cppMatrixMinusMatrix(
    Rcpp::NumericMatrix& inputMatrix,
    Rcpp::NumericMatrix& substrahendMatrix
) {
  const size_t width = inputMatrix.cols();   //długość widma
  const size_t height = inputMatrix.rows();  //liczba kwazarów
  const double* h_input = quasarCuda::RcppNumericMatrixToCArray(inputMatrix);
  const double* h_substrahend = quasarCuda::RcppNumericMatrixToCArray(substrahendMatrix);
  Rcpp::NumericMatrix result(inputMatrix.nrow(), inputMatrix.ncol());
  matrixSubstractMatrix(h_input, h_substrahend, &result[0], width, height);

  return result;
}

// [[Rcpp::export]]
SEXP cppMatrixDivideMatrix(
    Rcpp::NumericMatrix& inputMatrix,
    Rcpp::NumericMatrix& divisorMatrix
) {
  const size_t width = inputMatrix.cols();   //długość widma
  const size_t height = inputMatrix.rows();  //liczba kwazarów
  
  //TODO: check rows/cols?
  const double* h_input = quasarCuda::RcppNumericMatrixToCArray(inputMatrix);
  // for(int col = 0; col < width; ++col)
  //   for(int row = 0; row < height; ++row)
  //     Rcpp::Rcout << "m r " << row << " c " << col << " v " << h_input[col * height + row] << std::endl;
  const double* h_divisor = quasarCuda::RcppNumericMatrixToCArray(divisorMatrix);
  double* h_result = new double[width * height]; 
  
  matrixDivideMatrix(h_input, h_divisor, h_result, width, height);

  return quasarCuda::CArrayToRcppNumericMatrix(h_result, height, width);
}

//matrixMultiplyColVector
// [[Rcpp::export]]
SEXP cppMatrixMultiplyCol(
    Rcpp::NumericMatrix& inputMatrix,
    Rcpp::NumericVector& vector
) 
{
  const size_t width = inputMatrix.cols();   //długość widma
  const size_t height = inputMatrix.rows();  //liczba kwazarów
  const size_t length = vector.size();
  Rcpp::NumericMatrix result(width, height);
  double* output = quasarCuda::RcppNumericMatrixToCArray(result);
  double* input = quasarCuda::RcppNumericMatrixToCArray(inputMatrix);
  double* vectorA = quasarCuda::RcppNumericVectorToCArray(vector);
  matrixMultiplyColVector(input, output, vectorA, width, height, length);
  return result;
}

//TODO: fix!
// [[Rcpp::export]]
SEXP cppMatrixTranspose(
    Rcpp::NumericMatrix inputMatrix
) 
{
  const size_t width = inputMatrix.cols();   //długość widma
  const size_t height = inputMatrix.rows();  //liczba kwazarów;
//  double* h_idata = &inputMatrix[0];
  std::vector<double> h_input(width * height);
  // std::copy(inputMatrix.begin(), inputMatrix.end(), h_input.begin());
  std::vector<double> h_output(width * height);
  
  for(int j = 0; j < height; ++j) {
    for(int i = 0; i < width; ++i) 
      h_input[j * width + i] = inputMatrix(j, i);
  }
  
  for(int i = 0; i < 10; ++i)
    Rcpp::Rcout << inputMatrix[i] << std::endl;
  
  Rcpp::NumericMatrix result(width, height);
  matrixTranspose(&inputMatrix[0], &h_output[0], width, height);
   //result[i] = h_output[i];
 // std::copy(h_output.begin(), h_output.end(), result.begin());
 for(int j = 0; j < height; ++j) {
   for(int i = 0; i < width; ++i)
     result(j, i) = h_output[j * width + i];
     // h_input[j * width + i] = inputMatrix(j, i);
 }
 return result;
  
}


