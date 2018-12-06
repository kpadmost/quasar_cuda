#include <Rcpp.h>
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
  
  //TODO: check rows/cols?
  const double* h_input = quasarCuda::RcppNumericMatrixToCArray(inputMatrix);
  for(int col = 0; col < width; ++col)
    for(int row = 0; row < height; ++row)
      Rcpp::Rcout << "m r " << row << " c " << col << " v " << h_input[col * height + row] << std::endl;
  const double* h_substrahend = quasarCuda::RcppNumericMatrixToCArray(substrahendMatrix);
  double* h_result = new double[width * height]; 
  
  matrixSubstractMatrix(h_input, h_substrahend, h_result, width, height);
  // todo: add conversion function
  for(int i = 0; i < width; ++i)
    for(int j = 0; j < height; ++j)
      Rcpp::Rcout << "m r " << i << " c " << j << " v " << h_result[i * height + j] << std::endl;
  return quasarCuda::CArrayToRcppNumericMatrix(h_result, width, height);
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
  for(int col = 0; col < width; ++col)
    for(int row = 0; row < height; ++row)
      Rcpp::Rcout << "m r " << row << " c " << col << " v " << h_input[col * height + row] << std::endl;
  const double* h_divisor = quasarCuda::RcppNumericMatrixToCArray(divisorMatrix);
  double* h_result = new double[width * height]; 
  
  matrixDivideMatrix(h_input, h_divisor, h_result, width, height);

  return quasarCuda::CArrayToRcppNumericMatrix(h_result, width, height);
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