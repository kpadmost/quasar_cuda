#include <Rcpp.h>

#include "includes/continuum.h"

#include "RcppExtended.hpp"
using namespace Rcpp;

// [[Rcpp::export]]
List cppCalculateCfunDcfun(const NumericMatrix& wavelengthsMatrix, SEXP rReglin, SEXP rCReglin) {
  std::vector<double8> reglin = as<std::vector<double8> >(rReglin);
  std::vector<double8> cReglin = as<std::vector<double8> >(rCReglin);
  const size_t width = wavelengthsMatrix.rows();
  const size_t height = wavelengthsMatrix.cols();
  NumericMatrix continuumMatrix(width, height);
  NumericMatrix dContinuumMatrix(width, height);
  calculateDcfun(
    &wavelengthsMatrix[0],
                      &dContinuumMatrix[0],
                      &continuumMatrix[0],
                      &reglin[0],
                      &cReglin[0],
                      width,
                      height
  );
  return List::create(
    Named("cfun") = continuumMatrix,
    Named("dcfun") = dContinuumMatrix
  );
}

/*extern "C"
 void fixContinuumReglin(
 double8 *h_continuum_reglin,
 double8 *h_reglin,
 const size_t size
 )*/

// [[Rcpp::export]]
List cppFixReglin(SEXP cReglinVector, SEXP reglinVector) {
  std::vector<double8> cReglinVectorOutput = Rcpp::as< std::vector<double8> > (cReglinVector);
  std::vector<double8> reglinVectorOutput = Rcpp::as< std::vector<double8> > (reglinVector);
  const size_t vectorSize = cReglinVectorOutput.size();
  fixContinuumReglin(&cReglinVectorOutput[0], &reglinVectorOutput[0], vectorSize);
  return List::create(
    Named("cReglinFixed") = wrap(cReglinVectorOutput),
    Named("reglinFixed") = wrap(reglinVectorOutput)
  );
}

/*extern "C"
 void calculateContinuumFunction(
 const double *h_wavelengths,
 double *h_cfun,
 const double8 *h_reglin_vector,
 const size_t width,
 const size_t height
 )*/
// [[Rcpp::export]]
NumericMatrix cppCalculateContinuumMatrix(const NumericMatrix& wavelengthsMatrix, SEXP reglinVector) {
  const size_t width = wavelengthsMatrix.rows();
  const size_t height = wavelengthsMatrix.cols();
  std::vector<double8> continuumReglinVector = as<std::vector<double8> >(reglinVector);
  NumericMatrix output(width, height);
  calculateContinuumFunction(
    &wavelengthsMatrix[0], 
                      &output[0],
                      &continuumReglinVector[0],
                      width,
                      height
  );
  return output;
}


// [[Rcpp::export]]
NumericVector cppReduceContinuumChisq(NumericVector& chisq, const IntegerVector sizes) {
  const size_t spectrumsNumber = chisq.size();
  std::vector<size_t> sizesVector(sizes.begin(), sizes.end());
  NumericVector chisqCopy(chisq.begin(), chisq.end());
  reduceChisqs(
    &chisqCopy[0],
              &sizesVector[0],
              spectrumsNumber
  );
  
  return chisqCopy;
  
}