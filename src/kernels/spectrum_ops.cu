#include "cuda_tools.h"

#include <math_constants.h>

// generate code

__global__
void generateWavelengthsMatrix
	(
		double  *output_matrix,
		double4 *params_vector
	)
{
	// gid0 - numer widma kwazaru
	uint gid0 = blockIdx.x * blockDim.x + threadIdx.x;
	
	// parametry a, b i z kwazaru
  double4 abz_ = params_vector[gid0];
	
	// gid1 - numer elementu widma (indeks od 0 do 4095)
	uint gid1 = blockIdx.y * blockDim.y + threadIdx.y;
	uint idx = gid0 * gridDim.y * blockDim.y  + gid1;
	
	// Wyliczenie lambdy dla tego gid1
	double wl = abz_.x  + (abz_.y * (double)(gid1));
	wl = pow(10.0,wl);
	
	// Uwzglednienie przesuniecia kosmologicznego - przejscie do ukladu emitujacego z wykorzystaniem redshiftu
	output_matrix[idx] = wl / (abz_.z + 1.0);
	
	return;
}


__global__
void single_interpolation
	(
		const double *x_matrix,	// Długości dla widm
		const double *y_matrix,	// Widma po kolumnach
		const size_t  *sizes,			// Rozmiary widm
		const uint size,				// Ilość widm (ilość kolumn)
		const double *s_matrix,	// Długości fal widma, które dodajemy
		const double *t_matrix,	// Widmo, które dodajemy
		const uint to_add_spectrum_size, 	
		double *output_matrix		// Macierz zapisu wyniku sumy
	)
{
	// gid0 - numer widma kwazaru
	uint gid0 = blockIdx.x * blockDim.x + threadIdx.x;
	if(gid0 >= size)
	{
		return;
	}	
	
	uint spectrum_size = sizes[gid0];

	// Indek elementu ze wavelengths_matrix/spectrums_matrix
	uint idx = gid0;
	uint idx_max = idx + spectrum_size * size;	
	
	// Indeks elementu z to_add_wavelengths/to_add_spectrum
	uint to_add_idx = 0;
	uint to_add_idx_max = to_add_spectrum_size - 1;
	
	double wavelength = x_matrix[idx];
	double value = y_matrix[idx];
	
	double to_add_wl 	= s_matrix[to_add_idx];
	double to_add_wl_next 	= s_matrix[to_add_idx+1];		
	double to_add_value 	= t_matrix[to_add_idx];
	double to_add_value_next = t_matrix[to_add_idx+1];	
	

	int idx_max_flag = 0;
	int to_add_idx_max_flag = 0;

	while(1)
	{				
		if(wavelength >= to_add_wl)
		{			
			if(wavelength <= to_add_wl_next)
			{	
				double a = (to_add_value_next - to_add_value)/(to_add_wl_next - to_add_wl);
				double b = to_add_value - (a * to_add_wl);
				value = value + (a * wavelength + b);		
				output_matrix[idx] = value;		
				
				idx += size;
				// Sprawdzanie czy idx nie przekroczył idx_max
				// Zanim odczytamy dane.
				idx_max_flag = idx < idx_max ? 0 : 1;
				if(idx_max_flag)
				{
					break;
				}
				wavelength = x_matrix[idx];
				value = y_matrix[idx];
			}
			else
			{
				to_add_idx++;
				// Sprawdzanie czy to_add_idx nie przekroczył to_add_idx_max
				// Zanim odczytamy dane.
				to_add_idx_max_flag = to_add_idx < to_add_idx_max ? 0 : 1;
				if(to_add_idx_max_flag)
				{
					break;
				}
				to_add_wl = to_add_wl_next;				
				to_add_wl_next = s_matrix[to_add_idx+1];
				
				to_add_value = to_add_value_next;
				to_add_value_next = t_matrix[to_add_idx+1];
			}
		}
		else
		{	
			output_matrix[idx] = value;
		
			idx += size;
			// Sprawdzanie czy idx nie przekroczył idx_max
			// Zanim odczytamy dane.
			idx_max_flag = idx < idx_max ? 0 : 1;
			if(idx_max_flag)
			{
				break;
			}
			wavelength = x_matrix[idx];
			value = y_matrix[idx];
		}	
	}
	
	while(idx < idx_max)
	{
		value = y_matrix[idx];
		output_matrix[idx] = value; 
		idx += size;
	}	
}







extern "C"
void generateWavelengths(
  double* h_output,
  const double4* h_params,
  const size_t spectrum_number
) 
{
  double *d_output = 0;
  double4 *d_params = 0;
  checkCudaErrors(cudaSetDevice(0));
  const size_t matrix_size = spectrum_number * ASTRO_OBJ_SIZE * sizeof(double);
  checkCudaErrors(cudaMalloc((void**)&d_output, matrix_size));
  checkCudaErrors(cudaMalloc((void**)&d_params, spectrum_number * sizeof(double4)));
  checkCudaErrors(cudaMemcpy(d_params, h_params, spectrum_number * sizeof(double4), cudaMemcpyHostToDevice));
  dim3 threadsPerBlock(1, 512);
  dim3 blocksPerGrid(1, 1);
  blocksPerGrid.x = spectrum_number / threadsPerBlock.x;
  blocksPerGrid.y = ASTRO_OBJ_SIZE / threadsPerBlock.y;
  generateWavelengthsMatrix<<<blocksPerGrid, threadsPerBlock>>>(d_output, d_params);
  checkCudaErrors(cudaGetLastError());
  cudaDeviceSynchronize();
  checkCudaErrors(cudaMemcpy(h_output, d_output, matrix_size, cudaMemcpyDeviceToHost));
  cudaFree(d_output);
  cudaFree(d_params);
}

extern "C"
void singleInterpolation(
    const double *h_matrix_x,
    const double *h_matrix_y,
    const size_t *h_sizes_x,
    const size_t size, // number of quasars(coln)
    const double *h_matrix_s,
    const double *h_matrix_t,
    const size_t size_s,
    double *h_output
) {
  double *d_matrix_x, *d_matrix_y, *d_matrix_s, *d_matrix_t, *d_output;
  size_t *d_sizes_x;
  const size_t x_matrix_size = size * ASTRO_OBJ_SIZE * sizeof(double);
  const size_t s_matrix_size = size * size_s * sizeof(double);
  const size_t sizes_vector_size = size * sizeof(size_t);
  checkCudaErrors(cudaMalloc((void**)&d_output, x_matrix_size));
  checkCudaErrors(cudaMalloc((void**)&d_matrix_x, x_matrix_size));
  checkCudaErrors(cudaMalloc((void**)&d_matrix_y, x_matrix_size));
  checkCudaErrors(cudaMalloc((void**)&d_matrix_s, s_matrix_size));
  checkCudaErrors(cudaMalloc((void**)&d_matrix_t, s_matrix_size));
  
  checkCudaErrors(cudaMalloc((void**)&d_sizes_x, sizes_vector_size));
  
  
  
  checkCudaErrors(cudaMemcpy(d_matrix_x, h_matrix_x, x_matrix_size, cudaMemcpyHostToDevice));
  checkCudaErrors(cudaMemcpy(d_matrix_y, h_matrix_y, x_matrix_size, cudaMemcpyHostToDevice));
  checkCudaErrors(cudaMemcpy(d_matrix_s, h_matrix_s, s_matrix_size, cudaMemcpyHostToDevice));
  checkCudaErrors(cudaMemcpy(d_matrix_t, h_matrix_t, s_matrix_size, cudaMemcpyHostToDevice));
  checkCudaErrors(cudaMemcpy(d_sizes_x, h_sizes_x, sizes_vector_size, cudaMemcpyHostToDevice));
  
  const uint threadsPerBlock = BLOCK_DIM;
  const uint blocksPerGrid = calculateBlockNumber(size, threadsPerBlock);
  /*void single_interpolation
	(
		double *x_matrix,	// Długości dla widm
		double *y_matrix,	// Widma po kolumnach
		size_t  *sizes,			// Rozmiary widm
		const uint size,				// Ilość widm (ilość kolumn)
		const double *s_matrix,	// Długości fal widma, które dodajemy
		const double *t_matrix,	// Widmo, które dodajemy
		const uint to_add_spectrum_size, 	
		double *output_matrix		// Macierz zapisu wyniku sumy
	)*/
  single_interpolation<<<blocksPerGrid, threadsPerBlock>>>(
    d_matrix_x,
    d_matrix_y,
    d_sizes_x,
    size,
    d_matrix_s,
    d_matrix_t,
    size_s,
    d_output
    );
  checkCudaErrors(cudaGetLastError());
  checkCudaErrors(cudaMemcpy(h_output, d_output, x_matrix_size, cudaMemcpyDeviceToHost));
  cudaFree(d_matrix_x);
  cudaFree(d_matrix_y);
  cudaFree(d_matrix_s);
  cudaFree(d_matrix_t);
  cudaFree(d_output);
  cudaFree(d_sizes_x);
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