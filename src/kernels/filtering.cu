#include "cuda_tools.cc"

__global__
void count_if_not_inf
	(
		const double *input,	// Macierz wejściowa.
		const size_t width,
		const size_t height,
    size_t *result  // Wyniki, liczby elementów != INF dla każdegow wiersza.
	) 
{
	const uint gid = blockDim.x * blockIdx.x + threadIdx.x;		
	if(gid >= width)
	{
		return;
	}

	uint idx = gid;
	const uint idx_max = idx + height * width;
	
	uint count = 0;
	while(idx < idx_max)
	{
		if(input[idx] != CUDART_INF) ++count;
		idx += width;
	}
	result[gid] = count;
}

__global__ void copy_if_not_inf
	(
		const double * input,	// Macierz wejściowa.
		double * output,// Macierz wynikowa.
		const uint width,
		const uint height,
		const uint output_height	// Rozmiar wiersza macierzy output.
	)
{
	const uint gid = blockIdx.x * blockDim.x + threadIdx.x;		
	if(gid >= width)
	{
		return;
	}

	// indeksy 
	uint input_idx = gid;
	const uint input_end = input_idx + height * width;

	uint output_idx = gid;
	uint output_end = output_idx + output_height * width;

	for(; input_idx < input_end && output_idx < output_end; input_idx+=width)
	{
		double in = input[input_idx];
		if(in != CUDART_INF)
		{			
			output[output_idx] = in;
			output_idx+=width;
		}
		
	}	
}

extern "C"
void countNotInf(
  double *h_input,
  size_t *h_output,
  const size_t width,
  const size_t height
)
{
  // allocating kernel data
  double *d_input;
  size_t *d_output;
  const size_t matrix_size = width * height * sizeof(double);
  const size_t output_size = width * sizeof(size_t);
  checkCudaErrors(cudaMalloc((void**)&d_output, output_size));
  checkCudaErrors(cudaMalloc((void**)&d_input, matrix_size));
  // copying data
  checkCudaErrors(cudaMemcpy(d_input, h_input, matrix_size, cudaMemcpyHostToDevice));
  // instatiating kernel
  const dim3 threadsPerBlock(BLOCK_DIM, BLOCK_DIM, 1);
  const dim3 blocksPerGrid(width / threadsPerBlock.x, height / threadsPerBlock.y, 1);
  // calling kernel
  count_if_not_inf<<<blocksPerGrid, threadsPerBlock>>>(d_input, width, height, d_output);
  cudaDeviceSynchronize();
  checkCudaErrors(cudaGetLastError());
  //copying memory back
  checkCudaErrors(cudaMemcpy(h_output, d_output, output_size, cudaMemcpyDeviceToHost));
  // free memory
  cudaFree(d_input);
  cudaFree(d_output);
}

//consider CPU
extern "C"
void copyIfNotInf(
  const double *h_input,
  double *h_output,
  const size_t width,
  const size_t height,
  const size_t newHeight
)
{
  // allocating kernel data
  double *d_input, *d_output;
  const size_t input_matrix_size = width * height * sizeof(double);
  const size_t output_matrix_size = width * newHeight * sizeof(double);
  checkCudaErrors(cudaMalloc((void**)&d_output, output_matrix_size));
  checkCudaErrors(cudaMalloc((void**)&d_input, input_matrix_size));
  // copying data
  checkCudaErrors(cudaMemcpy(d_input, h_input, input_matrix_size, cudaMemcpyHostToDevice));
  // instatiating kernel
  const size_t threadsPerBlock = BLOCK_DIM;
  const size_t blocksPerGrid = calculateBlockNumber(width, threadsPerBlock);
  // calling kernel
  copy_if_not_inf<<<blocksPerGrid, threadsPerBlock>>>(d_input, d_output, width, height, newHeight);
  cudaDeviceSynchronize();
  checkCudaErrors(cudaGetLastError());
  //copying memory back
  checkCudaErrors(cudaMemcpy(h_output, d_output, output_matrix_size, cudaMemcpyDeviceToHost));
  // free memory
  cudaFree(d_input);
  cudaFree(d_output);
}

/* example func, copypaste as template
extern "C"
void template_f(
  double *h_input,
  double *h_output,
  const size_t width,
  const size_t height
)
{
  // allocating kernel data
  double *d_input, d_output;
  const size_t matrix_size = width * height * sizeof(double);
  checkCudaErrors(cudaMalloc((void**)&d_output, matrix_size));
  checkCudaErrors(cudaMalloc((void**)&d_input, matrix_size));
  // copying data
  checkCudaErrors(cudaMemcpy(d_input, h_input, matrix_size, cudaMemcpyHostToDevice));
  // instatiating kernel
  const dim3 threadsPerBlock(BLOCK_DIM, BLOCK_DIM, 1);
  const dim3 blocksPerGrid(width / threadsPerBlock.x, height / threadsPerBlock.y, 1);
  // calling kernel
  kernel_f<<<blocksPerGrid, threadsPerBlock>>>(d_input, d_output);
  cudaDeviceSynchronize();
  checkCudaErrors(cudaGetLastError());
  //copying memory back
  checkCudaErrors(cudaMemcpy(h_output, d_output, matrix_size, cudaMemcpyDeviceToHost));
  // free memory
  cudaFree(d_input);
  cudaFree(d_output);
}*/
