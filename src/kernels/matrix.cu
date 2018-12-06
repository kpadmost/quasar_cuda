#include <cuda_runtime.h>
#include <stdio.h>

#include "cuda_tools.cc"

#define THREADS_NUM 512

__global__ void matrix_log10
	(
		double *matrix,	// Macierz
		uint row_size,			// Rozmiar wiersza macierzy.
		uint col_size
	)
{
	// gid0 - numer wiersza macierzy input
	uint gid0 = blockIdx.x * blockDim.x + threadIdx.x;
	// gid1 - numer elementu w wierszu.
	uint gid1 = blockIdx.y * blockDim.y + threadIdx.y;
	if(gid1 >= row_size || gid0 >= col_size)
	{
		return;
	}
	
	uint idx = gid0 * row_size + gid1;
	double m = matrix[idx];
	m = log10(m);
	matrix[idx] = m;
}

__global__
void matrix_add_scalar
	(
		double * matrix,
		uint row_size, //TODO: change to width height
		uint col_size,
		const double scalar
	)
{
	// gid0 - numer wiersza macierzy input
	const uint gid0 = blockIdx.x * blockDim.x + threadIdx.x;
	// gid1 - numer elementu w wierszu.
	const uint gid1 = blockIdx.y * blockDim.y + threadIdx.y;

	if(gid1 >= row_size || gid0 >= col_size)
		return;
	
	const uint idx = gid0 * row_size + gid1;
	matrix[idx] = matrix[idx] + scalar;
}

__global__
void matrix_minus_matrix
	(
		const double* minuend,	// Macierz
					// Rozmiar wiersza macierzy.
		double* subtrahend,		// Macierz do odjęcia
		double* output,	// Wynik
		const uint row_size,
		const uint col_size
	)
{
	// gid0 - numer wiersza macierzy input
	uint gid0 = blockIdx.x * blockDim.x + threadIdx.x;
	// gid1 - numer elementu w wierszu.
	uint gid1 = blockIdx.y * blockDim.y + threadIdx.y;

	if(gid1 >= row_size || gid0 >= col_size)
	{
		return;
	}
	
	const uint idx = gid0 * row_size + gid1;
	output[idx] = minuend[idx] - subtrahend[idx];
}

__global__
void matrix_divide_matrix
	(
		const double* divident,	// Macierz
					// Rozmiar wiersza macierzy.
		const double* divisor,		// Macierz do odjęcia
		double* output,	// Wynik
		const uint row_size,
		const uint col_size
	)
{
	// gid0 - numer wiersza macierzy input
	uint gid0 = blockIdx.x * blockDim.x + threadIdx.x;
	// gid1 - numer elementu w wierszu.
	uint gid1 = blockIdx.y * blockDim.y + threadIdx.y;

	if(gid1 >= row_size || gid0 >= col_size)
	{
		return;
	}
	
	const uint idx = gid0 * row_size + gid1;
	output[idx] = divident[idx] / divisor[idx];
}

__global__
void matrix_multiply_vector
	(
		const double *matrix,
		const double *vector,	// Wektor, których zawiera co najmniej
		double *output,	// Wynik
		const uint row_size,
		const uint col_size
	)
{
	// gid0 - numer wiersza macierzy input
	uint gid0 = blockIdx.x * blockDim.x + threadIdx.x;
	// gid1 - numer elementu w wierszu (numer kolumny).
	uint gid1 = blockIdx.y * blockDim.y + threadIdx.y;

	if(gid1 >= row_size || gid0 >= col_size)
		return;
	
	uint idx = gid0 * row_size + gid1;
	uint col_idx = gid1;
	
	output[idx] = matrix[idx] * vector[col_idx];
}

__global__
void matrix_transpose(
  double *input,
  double *output,
  const uint row_size,
  const uint col_size
)
{

}


extern "C"
void matrixAddScalar(double* h_input, 
               const size_t width, 
               const size_t height,
               const double scalar
) {
  // initialize variables
  double *d_input;
  size_t size = width * height * sizeof(double);
  cudaError_t cudaStatus;
  // initialize device
  cudaStatus = cudaSetDevice(0);
  if (cudaStatus != cudaSuccess) {
      fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
      return;
  }
  // kernel memory allocation
  cudaStatus = cudaMalloc((void**)&d_input, size);
  if (cudaStatus != cudaSuccess) {
      fprintf(stderr, "cudaMalloc failed!");
      return;
  }
  
  cudaStatus = cudaMemcpy(d_input, h_input, size, cudaMemcpyHostToDevice);
  if (cudaStatus != cudaSuccess) {
      fprintf(stderr, "cudaMemcpy failed!");
      return;
  }

  // run kernel
  dim3 threadsPerBlock(16, 16);
  dim3 blocksPerGrid(1, 1);
  blocksPerGrid.x = ceil(double(width)/double(threadsPerBlock.x));
  blocksPerGrid.y = ceil(double(height)/double(threadsPerBlock.y));
  matrix_add_scalar<<<blocksPerGrid, threadsPerBlock>>>(
    d_input,
    width,
    height,
    scalar
  );
  
  // device to host memory copy
  cudaStatus = cudaMemcpy(h_input, d_input, size, cudaMemcpyDeviceToHost);
  if (cudaStatus != cudaSuccess) {
      fprintf(stderr, "cudaMemcpy failed!");
  }
  // free memory
  cudaFree(d_input);
}

extern "C"
void matrixLog10(double* h_input, const size_t width, size_t height)
{
    double *d_input = 0;
    size_t size = width * height * sizeof(double);
    //initialize device 
    checkCudaErrors(cudaSetDevice(0));
  
    // device memory allocation
    checkCudaErrors(cudaMalloc((void**)&d_input, size));
    // device memory copying
    checkCudaErrors(cudaMemcpy(d_input, h_input, size, cudaMemcpyHostToDevice));
    //kernel invocation
    dim3 threadsPerBlock(16, 16);
    dim3 blocksPerGrid(1, 1);
    blocksPerGrid.x = ceil(double(width)/double(threadsPerBlock.x));
    blocksPerGrid.y = ceil(double(height)/double(threadsPerBlock.y));
    matrix_log10<<<blocksPerGrid, threadsPerBlock>>>(d_input, width, height);
    checkCudaErrors(cudaGetLastError());
    //device to host memory copy
    checkCudaErrors(cudaMemcpy(h_input, d_input, size, cudaMemcpyDeviceToHost));
    
    cudaFree(d_input);
}

extern "C"
void matrixMultiplyColVector(
    const double* h_input,
    double* h_output,
    const double* h_vector,
    const size_t width,
    const size_t height,
    const size_t length
) 
{
  double *d_input = 0, *d_output = 0, *d_vector = 0;
  const size_t size = width * height * sizeof(double);
  const size_t vSize = length * sizeof(double);
  
  checkCudaErrors(cudaSetDevice(0));
  checkCudaErrors(cudaMalloc((void**)&d_input, size));
  checkCudaErrors(cudaMalloc((void**)&d_output, size));
  checkCudaErrors(cudaMalloc((void**)&d_vector, vSize)); //TODO: check!
  checkCudaErrors(cudaMemcpy(d_input, h_input, size, cudaMemcpyHostToDevice));
  checkCudaErrors(cudaMemcpy(d_vector, h_vector, vSize, cudaMemcpyHostToDevice));
  //kernel invocation
  dim3 threadsPerBlock(16, 16);
  dim3 blocksPerGrid(1, 1);
  blocksPerGrid.x = ceil(double(width)/double(threadsPerBlock.x));
  blocksPerGrid.y = ceil(double(height)/double(threadsPerBlock.y));
  matrix_multiply_vector<<<blocksPerGrid, threadsPerBlock>>>(
    d_input,
    d_vector, 
    d_output, 
    width, 
    height
  );
  checkCudaErrors(cudaGetLastError());
  checkCudaErrors(cudaMemcpy(h_output, d_output, size, cudaMemcpyDeviceToHost));
  cudaFree(d_input);
  cudaFree(d_output);
  cudaFree(d_vector);
}

extern "C"
void matrixSubstractMatrix(
  const double* h_input,
  const double* h_subtrahend,
  double* h_output,
  const size_t width,
  const size_t height
) 
{
   double* d_input = 0, *d_subtrahend = 0, *d_output = 0;
    size_t size = width * height * sizeof(double);
    //initialize device 
    checkCudaErrors(cudaSetDevice(0));
    
    // device memory allocation
    checkCudaErrors(cudaMalloc((void**)&d_input, size));
    
    checkCudaErrors(cudaMalloc((void**)&d_subtrahend, size));
    checkCudaErrors(cudaMalloc((void**)&d_output, size));
    
    // device memory copying
    checkCudaErrors(cudaMemcpy(d_input, h_input, size, cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(d_subtrahend, h_subtrahend, size, cudaMemcpyHostToDevice));
    //kernel invocation
    dim3 threadsPerBlock(512, 2);
    dim3 blocksPerGrid(1, 1);
    blocksPerGrid.x = ceil(double(width)/double(threadsPerBlock.x));
    blocksPerGrid.y = ceil(double(height)/double(threadsPerBlock.y));
    matrix_minus_matrix<<<blocksPerGrid, threadsPerBlock>>>(
      d_input,
      d_subtrahend,
      d_output,
      width, 
      height
    );
    //device to host memory copy
    checkCudaErrors(cudaMemcpy(h_output, d_output, size, cudaMemcpyDeviceToHost));
    cudaFree(d_input);
    cudaFree(d_subtrahend);
    cudaFree(d_output);
}


extern "C"
void matrixDivideMatrix(
    double* h_divided,
    double* h_divisor,
    double* h_output,
    const size_t width,
    const size_t height
)
{
   double* d_divided = 0, *d_divisor = 0, *d_output = 0;
    size_t size = width * height * sizeof(double);
    //initialize device 
    checkCudaErrors(cudaSetDevice(0));
    
    // device memory allocation
    checkCudaErrors(cudaMalloc((void**)&d_divided, size));
    
    checkCudaErrors(cudaMalloc((void**)&d_divisor, size));
    checkCudaErrors(cudaMalloc((void**)&d_output, size));
    
    // device memory copying
    checkCudaErrors(cudaMemcpy(d_divided, h_divided, size, cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(d_divisor, h_divisor, size, cudaMemcpyHostToDevice));
    //kernel invocation
    dim3 threadsPerBlock(512, 2);
    dim3 blocksPerGrid(1, 1);
    blocksPerGrid.x = ceil(double(width)/double(threadsPerBlock.x));
    blocksPerGrid.y = ceil(double(height)/double(threadsPerBlock.y));
    matrix_divide_matrix<<<blocksPerGrid, threadsPerBlock>>>(
      d_divided,
      d_divisor,
      d_output,
      width, 
      height
    );
    //device to host memory copy
    checkCudaErrors(cudaMemcpy(h_output, d_output, size, cudaMemcpyDeviceToHost));
    cudaFree(d_divided);
    cudaFree(d_divisor);
    cudaFree(d_output);
}



extern "C"
void matrixTranspose(
    double *h_input,
    double *h_output,
    const size_t realWidth,
    const size_t realHeight
) 
{
    double *d_input = 0, *d_output = 0;
    // for efficient kernel computations, we should expand matrices
    
    size_t size = realWidth * realHeight * sizeof(double);
    //initialize device 
    checkCudaErrors(cudaSetDevice(0));
  
    // device memory allocation
    checkCudaErrors(cudaMalloc((void**)&d_input, size));
    checkCudaErrors(cudaMalloc((void**)&d_output, size));
    // device memory copying
    checkCudaErrors(cudaMemcpy(d_input, h_input, size, cudaMemcpyHostToDevice));
    //kernel invocation
    dim3 threadsPerBlock(512, 2);
    dim3 blocksPerGrid(1, 1);
    blocksPerGrid.x = ceil(double(realWidth)/double(threadsPerBlock.x));
    blocksPerGrid.y = ceil(double(realHeight)/double(threadsPerBlock.y));
    matrix_transpose<<<blocksPerGrid, threadsPerBlock>>>(d_input, d_output, realWidth, realHeight);
    checkCudaErrors(cudaGetLastError());
    
    //device to host memory copy
    checkCudaErrors(cudaMemcpy(h_output, d_output, size, cudaMemcpyDeviceToHost));
    
    cudaFree(d_input);
    cudaFree(d_output);
}