#ifndef CUDA_TOOLS_CC
#define CUDA_TOOLS_CC


#include "cuda_rcpp_common.h"

#include <device_launch_parameters.h>
#include <stdio.h>
#include <iostream>
#include <math_constants.h>



#define BLOCK_DIM 32
#define ASTRO_OBJ_SIZE 4096
#define C 299792458.0


template<typename T> 
struct cuda_plus {
  __device__
  inline T operator()(const T a, const T b) {
    return a + b;
  } 
};

template<typename T> 
struct cuda_minus{
  __device__
  inline T operator()(const T a, const T b) {
    return a - b;
  } 
};

template <typename F>
__device__
double4 double4Op(const double4 a, const double4 b, F f) {
  double4 result;
  result.x = f(a.x, b.x);
  result.y = f(a.y, b.y);
  result.z = f(a.z, b.z);
  result.w = f(a.w, b.w);
  return result;
}


#define checkCudaErrors(val) check( (val), #val, __FILE__, __LINE__)

template<typename T>
void check(T err, const char* const func, const char* const file, const int line) {
  if (err != cudaSuccess) {
    std::cerr << "CUDA error at: " << file << ":" << line << std::endl;
    std::cerr << cudaGetErrorString(err) << " enum: " << err << " " << func << std::endl;
    throw std::runtime_error(cudaGetErrorString(err));
  }
}

inline uint calculateBlockNumber(const size_t size, const uint thrPerBlock) {
  const uint rem = (size % thrPerBlock);
  uint global = size;
  if(rem > 0 )
    global += thrPerBlock - rem ;
  return global / thrPerBlock;
}

namespace nv_tools {

inline cudaError_t init_cuda() {
  cudaError_t cudaStatus;
  
  // Choose which GPU to run on, change this on a multi-GPU system.
  cudaStatus = cudaSetDevice(0);
  if (cudaStatus != cudaSuccess) {
    fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
    
  }
  return cudaStatus;
}





inline cudaError_t cuda_malloc_wrap(void **devPtr, size_t size) {
  cudaError_t cudaStatus = cudaMalloc(devPtr, size);
  if (cudaStatus != cudaSuccess) {
    fprintf(stderr, "cudaMalloc failed!");
  }
  return cudaStatus;
}

inline cudaError_t cuda_memcpy_wrap( void* dst, const void* src, size_t count, cudaMemcpyKind kind ) {
  cudaError_t cudaStatus = cudaMemcpy(dst, src, count, kind);
  if(cudaStatus != cudaSuccess) {
    fprintf(stderr, "cudaMalloc failed!");
  }
  return cudaStatus;
}
}



#endif