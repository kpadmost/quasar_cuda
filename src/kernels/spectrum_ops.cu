#include "cuda_tools.cc"

// generate code



extern "C"
void generateWavelengths(
  double* h_output,
  double4* params,
  const size_t spectrum_size
) 
{
  checkCudaErrors(cudaSetDevice(0));
}