#ifndef RCPPEXTENDED_HPP
#define RCPPEXTENDED_HPP

#include <RcppCommon.h>

#include <builtin_types.h>

namespace Rcpp {
  template <> double4 as (SEXP);
  template <> SEXP wrap(const double4& obj);
}
#include <Rcpp.h>
#endif

