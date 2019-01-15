
#include "RcppExtended.hpp"


namespace Rcpp {
  template <> double8 as (SEXP vec_) {
    NumericVector vector(vec_);
    double8 res;
    res.x1 = vector[0];
    res.x2 = vector[1];
    res.x3 = vector[2];
    res.x4 = vector[3];
    res.x5 = vector[4];
    res.x6 = vector[5];
    res.x7 = vector[6];
    res.x8 = vector[7];
    
    return res;
  };
  
  template <> SEXP wrap(const double8& obj) {
    NumericVector result(8);
    result[0] = obj.x1;
    result[1] = obj.x2;
    result[2] = obj.x3;
    result[3] = obj.x4;
    result[4] = obj.x5;
    result[5] = obj.x6;
    result[6] = obj.x7;
    result[7] = obj.x8;
    return result;
  }
  
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
  
  template <> double2 as (SEXP vec_) {
    NumericVector vector(vec_);
    double2 res;
    res.x = vector[0];
    res.y = vector[1];
    return res;
  };
  
  template <> SEXP wrap(const double2& obj) {
    NumericVector result(2);
    result[0] = obj.x;
    result[1] = obj.y;
    return result;
  }
}
