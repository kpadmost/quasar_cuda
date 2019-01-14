#include "cuda_tools.cc"

__global__
void moving_average
	(
		const double *input,
		double *output,
		const size_t width,
		const uint *cols_heights,
		const uint window_width
	)
{
	// column
	uint gid = blockIdx.x * blockDim.x + threadIdx.x;
	uint idx = gid;

	if(gid >= width)
	{
		return;
	}

	uint col_height = cols_heights[gid];
		
	unsigned long end = idx + col_height * width;
	uint i = 0;
	double lastSum = 0;
	double result = 0;

	//
	// Dla pierwszych window_width - 1 elementów [0; window_width - 1)
	//

	// Wartość średniej kroczącej dla pierwsztch window_width elementów jest
	// wyznaczana tak, że dla i-tego elementu jest równe śrędniej arytmetycznej
	// od elementu o indeksie 0 do elementu o indeksie i+window_width
	while(i < window_width && idx < end)
	{
		lastSum += input[idx];			
		idx += width;	
		i++;
	}

	i = 0;
	while(i < window_width && idx < end)
	{	
		result = lastSum / ((double)(window_width + i));
		output[idx - (window_width * width)] = result;
		
		double new_v = input[idx];
		lastSum = lastSum + new_v;

		idx += width;	
		i++;
	}

	//
	// Dla elementów o indeksie z przedziału [window_width; ilość_elementów - window_width]
	//
	
	double fwindow_width = (double)(2 * window_width);
	while(idx < end)
	{			
		result = lastSum / fwindow_width;
		output[idx - (window_width * width)] = result;

		double new_v = input[idx];
		double old = input[idx - (2 * window_width * width)];
		lastSum = lastSum - old + new_v;	

		idx += width;	
	}
	
	//
	// Dla elementów o indeksie z przedziału (ilość_elementów - window_width; ilość_elementów]
	//

	// Wartość średniej kroczącej dla i-tego elementu jest średnią arytmetyczną elementów
	// o indeksach w przedziale (i; ilość_elementów)

	lastSum = 0.0f;
	idx -= 2 * window_width * width;
	while(idx < end)
	{		
		lastSum += input[idx];
		idx += width;	
	}

	idx -= window_width * width;	
	i = 2 * window_width;
	while(idx < end)
	{			
		result = lastSum / ((double)(i));
		output[idx] = result;
		
		double old = input[idx - window_width * width];
		lastSum = lastSum - old;	

		idx += width;	
		i--;
	}
}


extern "C"
void movingAverage( // centered moving average
  const double *h_input,
  double *h_output,
  const uint *h_cols,
  const size_t width,
  const size_t height,
  const uint window
) 
{
    // malloc
    double *d_input, *d_output;
    uint *d_cols;
    const size_t matrix_size = width * height * sizeof(double);
    const size_t cols_vec_size = width * sizeof(size_t);
    checkCudaErrors(cudaMalloc((void**)&d_input, matrix_size));
    checkCudaErrors(cudaMalloc((void**)&d_output, matrix_size));
    checkCudaErrors(cudaMalloc((void**)&d_cols, cols_vec_size));
    
    checkCudaErrors(cudaMemcpy(d_input, h_input, matrix_size, cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(d_cols, h_cols, cols_vec_size, cudaMemcpyHostToDevice));
    
    const uint threadsPerBlock = BLOCK_DIM;
    const uint blocksPerGrid =  calculateBlockNumber(width, BLOCK_DIM);
    moving_average<<<blocksPerGrid, threadsPerBlock>>>(d_input, d_output, width, d_cols, window);
    cudaDeviceSynchronize();
    checkCudaErrors(cudaGetLastError());
    checkCudaErrors(cudaMemcpy(h_output, d_output, matrix_size, cudaMemcpyDeviceToHost));
    cudaFree(d_input);
    cudaFree(d_output);
    cudaFree(d_cols);
}