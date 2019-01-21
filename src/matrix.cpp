#include <Rcpp.h>
#include <iterator>
#include <array>
#include <algorithm>
#include "includes/matrix.h"

// [[Rcpp::export]]
Rcpp::NumericMatrix cppMatrixLog10(const Rcpp::NumericMatrix& inputMatrix) {
  const size_t width = inputMatrix.rows();   //number of quas
  const size_t height = inputMatrix.cols();  //length of spectrum
  Rcpp::NumericMatrix result = Rcpp::clone(inputMatrix);
  matrixLog10(&result[0], width, height);

  return result;
}

// [[Rcpp::export]]
SEXP cppMatrixAddScalar(Rcpp::NumericMatrix inputMatrix, double scalar) {
  const size_t width = inputMatrix.cols();   //długość widma
  const size_t height = inputMatrix.rows();  //liczba kwazarów
  
  matrixAddScalar(&inputMatrix[0], width, height, scalar);
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
  Rcpp::NumericMatrix result(inputMatrix.nrow(), inputMatrix.ncol());
  matrixSubstractMatrix(&inputMatrix[0], &substrahendMatrix[0], &result[0], width, height);

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
  double* h_result = new double[width * height]; 
  Rcpp::NumericMatrix output(height, width);
  matrixDivideMatrix(&inputMatrix[0], &divisorMatrix[0], &output[0], width, height);

  return output;
}

//matrixMultiplyColVector
// [[Rcpp::export]]
SEXP cppMatrixMultiplyCol(
    const Rcpp::NumericMatrix& inputMatrix,
    const Rcpp::NumericVector& vector
) 
{
  const size_t width = inputMatrix.rows();   //długość widma
  const size_t height = inputMatrix.cols();  //liczba kwazarów
  const size_t length = vector.size();
  Rcpp::NumericMatrix result(width, height);
  matrixMultiplyColVector(&inputMatrix[0], &result[0], &vector[0], width, height);
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


