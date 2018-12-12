
#include "RcppExtended.hpp"


namespace Rcpp {
  template <> double4 as (SEXP vec_) {
    NumericVector vector(vec_);
    double4 res;
    res.x = vector[0];
    res.y = vector[1];
    res.z = vector[2];
    res.w = vector[3];
    return res;
  };
    
  template <> SEXP wrap(const double4& obj) {
    NumericVector result(4);
    result[0] = obj.x;
    result[1] = obj.y;
    result[2] = obj.z;
    result[3] = obj.w;
    return result;
  }
}
