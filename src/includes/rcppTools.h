#ifndef RCPP_TOOLS_H
#define RCPP_TOOLS_H

#include <Rcpp.h>

Rcpp::NumericVector cppChisq(
    const Rcpp::NumericMatrix& x,
    const Rcpp::NumericMatrix& y, 
    const Rcpp::NumericMatrix& e,
    const Rcpp::IntegerVector& sizes
);

#endif