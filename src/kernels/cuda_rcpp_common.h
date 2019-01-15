#ifndef CUDA_RCPP_COMMON_H
#define CUDA_RCPP_COMMON_H
#include <cuda_runtime.h>


#ifdef __CUDACC__
typedef struct __align__(16)
#else 
typedef struct
#endif
  double8{
  double x1;
  double x2;
  double x3;
  double x4;
  double x5;
  double x6;
  double x7;
  double x8;
} double8;

#endif