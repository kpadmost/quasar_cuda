#include <Rcpp.h>
using namespace Rcpp;

#include "RcppExtended.hpp"

#include "includes/quasar_spectrum.h"
#include "includes/quasar_fefit.h"
#include "includes/rcppTools.h"
#include "includes/tools.h"
#include <cmath>
/*loadDefaultFeFitParams <- function() {
 list(
 fwhmn = 1600.0,
 fwhmt = 900.0,
 feFitRange = c(2200.0, 2650.0),
 isSubC = FALSE,
 fitType="fwin"
 )*/

NumericVector cppGaussPeak(const NumericVector& x, const double sigmaConv, const double miConv) {
  const double a = 1.0 / (sigmaConv * std::sqrt(2 * M_PI));
  NumericVector result = clone(x);
  std::transform(result.begin(), result.end(), result.begin(), 
    [&](const double& x) -> double {
      return a * std::exp(-0.5 * std::pow((x - miConv) / sigmaConv, 2));
    });
  return result;
}

//TODO:refactor
inline NumericVector cpuConvolve(const NumericVector& signal, const NumericVector& kernel, bool centered)
{
  Rcpp::NumericVector result(signal.size() + kernel.size() - 1, 0);
  
  for (size_t n = 0; n < result.size(); n++)
  {
    size_t kmin, kmax, k;
    
    kmin = (n >= kernel.size() - 1) ? n - (kernel.size() - 1) : 0;
    kmax = (n < signal.size() - 1) ? n : signal.size() - 1;
    
    for (k = kmin; k <= kmax; k++)
    {
      result[n] += signal[k] * kernel[n - k];
    }
  }
  
  if (centered)
  {
    size_t kernel_center = kernel.size() % 2 == 0 ? (kernel.size() - 1) / 2 : kernel.size() / 2;
    result.erase(result.begin(), result.begin() + kernel_center);
    result.erase(result.end() - (kernel.size() - 1 - kernel_center), result.end());
  }
  
  return result;
}

//fwhmConv <- (feFitParams$fwhmn ^ 2 - feFitParams$fwhmt ^ 2) ^ 0.5

//sigmaConv <- fwhmConv / (2 * ((2 * log(2)) ^ 0.5)) / C * 1e3

inline double calculateSigmaConvolution(const double fwhmn, const double fwhmt) {
  const double fwhmConv = std::sqrt(std::pow(fwhmn, 2) - std::pow(fwhmt, 2));
  return fwhmConv / (2 * (std::sqrt(2 * std::log(2)))) / C * 1e3;
}

/*loadDefaultFeFitParams <- function() {
 list(
 fwhmn = 1600.0,
 fwhmt = 900.0,
 feFitRange = c(2200.0, 2650.0),
 isSubC = FALSE,
 fitType="fwin"
 )*/
// [[Rcpp::export]]
NumericMatrix cppExpandTemplate(
    const NumericMatrix& wavelengths,
    const IntegerVector& sizes,
    const NumericVector& feWavelengths, 
    const NumericVector& feSpectrum,
    const List& feFitParams
) {
  const size_t width = wavelengths.rows();
  const size_t height = wavelengths.cols();
  std::vector<size_t> sizesVector = as<std::vector<size_t> > (sizes);
  NumericVector feWavelengthsLog10 = clone(feWavelengths);
  std::transform(feWavelengthsLog10.begin(), feWavelengthsLog10.end(), feWavelengthsLog10.begin(),
    [](const double x) -> double {
        return std::log10(x);
    });
  const double sigmaConv = calculateSigmaConvolution(
    as<double>(feFitParams["fwhmn"]),
    as<double>(feFitParams["fwhmt"])
  );
  const double miConv = feWavelengthsLog10[feWavelengthsLog10.size() / 2];
  NumericVector feWavelengthsSmooth = cppGaussPeak(feWavelengthsLog10, sigmaConv, miConv);
  NumericVector feSpectrumConvolved = cpuConvolve(feSpectrum, feWavelengthsSmooth, true);
  NumericMatrix spectrumZeros(width, height);
  NumericMatrix results(width, height);
  singleInterpolation(
    &wavelengths[0],
    &spectrumZeros[0],
    &sizesVector[0],
    width,
    &feWavelengths[0],
    &feSpectrumConvolved[0],
    feWavelengths.size(),
    &results[0]
  );
  return results;
}

// [[Rcpp::export]]
NumericVector cppCalculateFeScaleRates(
  const NumericMatrix& spectrumsMatrix,
  const NumericMatrix& templateMatrix,
  const IntegerVector& sizes,
  SEXP feFitParams
) {
  const size_t width = spectrumsMatrix.rows();
  const size_t height = spectrumsMatrix.cols();
  std::vector<size_t> sizesVector = as<std::vector<size_t> >(sizes);
  NumericVector result(width);
  calculateReglinSimplified(
    &templateMatrix[0],
    &spectrumsMatrix[0],
    width,
    height,
    &sizesVector[0],
    &result[0]
  );
  return result;
}

// [[Rcpp::export]]
NumericVector cppCalculateFeReducesChisq(
  const NumericMatrix& spectrumsMatrix,
  const NumericMatrix& templateMatrix,
  const NumericMatrix& errorsMatrix,
  const IntegerVector& sizes
) {
  const size_t width = spectrumsMatrix.rows();
  const size_t height = spectrumsMatrix.cols();
  std::vector<size_t> sizesVector = as<std::vector<size_t> >(sizes);
  NumericVector chisqVector = cppChisq(spectrumsMatrix, templateMatrix, errorsMatrix, sizes);
  reduceFeChisq(&chisqVector[0], &sizesVector[0], width);
  return chisqVector;
}

