#ifndef RCPPEXTENDED_HPP
#define RCPPEXTENDED_HPP

#include <RcppCommon.h>

#include <builtin_types.h>

namespace Rcpp {
  template <> double4 as (SEXP);
  template <> SEXP wrap(const double4& obj);
  
  // use builtin cudas
  template <> double2 as (SEXP);
  template <> SEXP wrap(const double2& obj);
}
#include <Rcpp.h>
#endif

