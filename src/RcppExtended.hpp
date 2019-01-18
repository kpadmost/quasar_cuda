#ifndef RCPPEXTENDED_HPP
#define RCPPEXTENDED_HPP

#include <RcppCommon.h>
#include "kernels/cuda_rcpp_common.h"
#include <builtin_types.h>
const double C = 299792458.0;

namespace Rcpp {
template <> double8 as (SEXP);
  template <> SEXP wrap(const double8& obj);

  template <> double4 as (SEXP);
  template <> SEXP wrap(const double4& obj);
  
  // use builtin cudas
  template <> double2 as (SEXP);
  template <> SEXP wrap(const double2& obj);
}
#include <Rcpp.h>
#endif

