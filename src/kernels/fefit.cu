#include "cuda_tools.h"

__global__
void reduce_fe_chisqs
	(
		double *chisqs, 		// Bufor z chisq
		const size_t * filtered_sizes, 	// Ilość znaczących elementów po filtracji
		const size_t size				// Ilość chisq i filtered_sizes
	)
{
	// gid0 - numer elementu z chisqs (czyli jednego chisq)
	const uint gid0 = blockDim.x * blockIdx.x + threadIdx.x;

	if(gid0 >= size)
		return;
	
	chisqs[gid0] = chisqs[gid0] / ((double)filtered_sizes[gid0] - 1.0);	
}

extern "C"
void reduceFeChisq(
  double *h_chisq,
  const size_t *h_sizes,
  const size_t width
)
{
  // allocating kernel data
  double *d_chisq;
  size_t *d_sizes;
  const size_t double_vector_size = width * sizeof(double);
  const size_t size_t_vector_size = width * sizeof(size_t);
  
  checkCudaErrors(cudaMalloc((void**)&d_chisq, double_vector_size));
  checkCudaErrors(cudaMalloc((void**)&d_sizes, size_t_vector_size));
  // copying data
  checkCudaErrors(cudaMemcpy(d_chisq, h_chisq, double_vector_size, cudaMemcpyHostToDevice));
  checkCudaErrors(cudaMemcpy(d_sizes, h_sizes, size_t_vector_size, cudaMemcpyHostToDevice));
  // instatiating kernel
  const size_t threadsPerBlock = BLOCK_DIM;
  const size_t blocksPerGrid = width / threadsPerBlock;
  // calling kernel
  reduce_fe_chisqs<<<blocksPerGrid, threadsPerBlock>>>(d_chisq, d_sizes, width);
  cudaDeviceSynchronize();
  checkCudaErrors(cudaGetLastError());
  //copying memory back
  checkCudaErrors(cudaMemcpy(h_chisq, d_chisq, double_vector_size, cudaMemcpyDeviceToHost));
  // free memory
  cudaFree(d_chisq);
  cudaFree(d_sizes);
}