#include "cuda_tools.cc"

#include <cfloat>

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


// Filtruje tylko te dane, których odpowiadająca wartość
// długości fali z wavelengths_matrix mieści się, w jakimś oknie widmowym.
//
__global__ void filter_with_wavelength_windows
	(
		double *wavelengths_matrix,	// Długości fal widm
		double *spectrums_matrix,	// Widma
		double *errors_matrix,		// Błędy pomiaru widm
		const size_t  *sizes,		 	// Rozmiary widm w spectrums_matrix
		const double2 *windows,	// Okna	długości fal
		const size_t windows_size 		// Ilość okien
	)	
{
	// gid0 - numer widma kwazaru
	const uint gid0 = blockIdx.x * blockDim.x + threadIdx.x;
	// gid1 - numer elementu w widmie
	const uint gid1 = blockIdx.y * blockDim.y + threadIdx.y;
	// indeks 
	uint idx = (gid0 * ASTRO_OBJ_SIZE) + gid1;

	uint size = sizes[idx % (gridDim.x *  blockDim.x)];

	// Zero wskazuje na brak dopasownia do któregokolwiek okna

	double wl_result = wavelengths_matrix[idx];
	

	bool global_flag = false;
	for(int window_idx = 0; window_idx < windows_size; window_idx++)
	{
    const double2 window = windows[window_idx];
		const int window_flag = wl_result >= window.x && wl_result <= window.y;
		  
		// Jeżeli oba warunki spełnione to mamy dopasowania
		// W przeciwnym przypadku zostawiamy to co było.

		//global_flag = select(global_flag, 1, window_flag == 1);
		if(window_flag == 1) {
		  global_flag = true;
		  break; // test performance
		}
		  
	}
	const uint row_idx = idx / (gridDim.x * blockDim.x);
	if(!global_flag || row_idx >= size) {
	  wavelengths_matrix[idx] = CUDART_INF;
	  spectrums_matrix[idx] = CUDART_INF;
	  errors_matrix[idx] = CUDART_INF;
	}
}


// set all as infinity
__global__ void filter_nonpositive
	(		
		double * spectrums_matrix,	// Widma
		double * a_matrix,	// Dowolna macierz wielkości spectrums_matrix
		double * b_matrix,	// Dowolna macierz wielkości spectrums_matrix
		const size_t *sizes 	// Rozmiary widm w spectrums_matrix
	)
{
	const uint gid0 = blockIdx.x * blockDim.x + threadIdx.x;
	const uint gid1 = blockIdx.y * blockDim.y + threadIdx.y;
	// Indeks elementu
	const uint idx = (gid0 * ASTRO_OBJ_SIZE) + gid1;
	
	const uint size = sizes[idx % gridDim.x * blockDim.x ];
	const uint row_idx = idx / (gridDim.x * blockDim.x);

	const double spectrum = spectrums_matrix[idx];
  
  if(row_idx >= size || spectrum < DBL_MIN) {
    spectrums_matrix[idx] = CUDART_INF;
    a_matrix[idx] = CUDART_INF;
    b_matrix[idx] = CUDART_INF;
  }
}


extern "C"
void countNotInf(
  const double *h_input,
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
  const size_t threadsPerBlock = BLOCK_DIM;
  const size_t blocksPerGrid = calculateBlockNumber(width, threadsPerBlock);
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

extern "C"
void filterNonpositive(
  double *h_spectrums,
  double *h_wavelengths,
  double *h_errors,
  const size_t *h_real_sizes,
  const size_t width
)
{
  double *d_spectrums, *d_wavelengths, *d_errors;
  size_t *d_real_sizes;
  const size_t matrix_size = width * ASTRO_OBJ_SIZE * sizeof(double);
  const size_t sizes_v_size = width * sizeof(size_t);
  checkCudaErrors(cudaMalloc((void**)&d_spectrums, matrix_size));
  checkCudaErrors(cudaMalloc((void**)&d_wavelengths, matrix_size));
  checkCudaErrors(cudaMalloc((void**)&d_errors, matrix_size));
  checkCudaErrors(cudaMalloc((void**)&d_real_sizes, sizes_v_size));
  // copying data
  checkCudaErrors(cudaMemcpy(d_spectrums, h_spectrums, matrix_size, cudaMemcpyHostToDevice));
  checkCudaErrors(cudaMemcpy(d_wavelengths, h_wavelengths, matrix_size, cudaMemcpyHostToDevice));
  checkCudaErrors(cudaMemcpy(d_errors, h_errors, matrix_size, cudaMemcpyHostToDevice));
  checkCudaErrors(cudaMemcpy(d_real_sizes, h_real_sizes, sizes_v_size, cudaMemcpyHostToDevice));
  // instatiating kernel
  const dim3 threadsPerBlock(1, BLOCK_DIM * 4, 1);
  const dim3 blocksPerGrid(width / threadsPerBlock.x, ASTRO_OBJ_SIZE / threadsPerBlock.y, 1);
  // calling kernel
  filter_nonpositive<<<blocksPerGrid, threadsPerBlock>>>(
    d_spectrums,
    d_wavelengths, 
    d_errors, 
    d_real_sizes
  );
  checkCudaErrors(cudaGetLastError());
  cudaDeviceSynchronize();
  //copying memory back
  checkCudaErrors(cudaMemcpy(h_spectrums, d_spectrums, matrix_size, cudaMemcpyDeviceToHost));
  checkCudaErrors(cudaMemcpy(h_wavelengths, d_wavelengths, matrix_size, cudaMemcpyDeviceToHost));
  checkCudaErrors(cudaMemcpy(h_errors, d_errors, matrix_size, cudaMemcpyDeviceToHost));
  // free memory
  cudaFree(d_spectrums);
  cudaFree(d_wavelengths);
  cudaFree(d_errors);
  cudaFree(d_real_sizes);
}


// filtering based on marking all of vals which not fit in windows as Inf
extern "C"
void filterWithWavelengthsWindows(
  double *h_wavelengths,
  double *h_spectrums,
  double *h_errors,
  const size_t *h_sizes,
  const double2 *h_windows,
  const size_t windows_number,
  const size_t width,
  const size_t height
)
{
  // allocating kernel data
  double *d_wavelengths, *d_spectrums, *d_errors;
  size_t *d_sizes;
  double2 *d_windows;
  const size_t windows_size = windows_number * sizeof(double2);
  const size_t matrix_size = width * height * sizeof(double);
  const size_t vector_size = width  * sizeof(size_t);
  checkCudaErrors(cudaMalloc((void**)&d_wavelengths, matrix_size));
  checkCudaErrors(cudaMalloc((void**)&d_spectrums, matrix_size));
  checkCudaErrors(cudaMalloc((void**)&d_errors, matrix_size));
  checkCudaErrors(cudaMalloc((void**)&d_sizes, vector_size));
  checkCudaErrors(cudaMalloc((void**)&d_windows, windows_size));
  // copying data
  checkCudaErrors(cudaMemcpy(d_wavelengths, h_wavelengths, matrix_size, cudaMemcpyHostToDevice));
  checkCudaErrors(cudaMemcpy(d_spectrums, h_spectrums, matrix_size, cudaMemcpyHostToDevice));
  checkCudaErrors(cudaMemcpy(d_errors, h_errors, matrix_size, cudaMemcpyHostToDevice));
  checkCudaErrors(cudaMemcpy(d_sizes, h_sizes, vector_size, cudaMemcpyHostToDevice));
  checkCudaErrors(cudaMemcpy(d_windows, h_windows, windows_size, cudaMemcpyHostToDevice));
  // instatiating kernel
  const dim3 threadsPerBlock(1, BLOCK_DIM, 1);
  const dim3 blocksPerGrid(width / threadsPerBlock.x, height / threadsPerBlock.y, 1);

  filter_with_wavelength_windows<<<blocksPerGrid, threadsPerBlock>>>(
    d_wavelengths,
    d_spectrums,
    d_errors,
    d_sizes,
    d_windows,
    windows_number
  );
  cudaDeviceSynchronize();
  checkCudaErrors(cudaGetLastError());
  //copying memory back
  checkCudaErrors(cudaMemcpy(h_wavelengths, d_wavelengths, matrix_size, cudaMemcpyDeviceToHost));
  checkCudaErrors(cudaMemcpy(h_spectrums, d_spectrums, matrix_size, cudaMemcpyDeviceToHost));
  checkCudaErrors(cudaMemcpy(h_errors, d_errors, matrix_size, cudaMemcpyDeviceToHost));
  // free memory
  cudaFree(d_wavelengths);
  cudaFree(d_spectrums);
  cudaFree(d_errors);
  cudaFree(d_sizes);
  cudaFree(d_windows);
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
