#include "cuda_tools.cc"

#include <math_constants.h>
#include <cfloat>
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
		double *x_matrix,	// Długości dla widm
		double *y_matrix,	// Widma po kolumnach
		size_t  *sizes,			// Rozmiary widm
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
		const size_t &windows_size 		// Ilość okien
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
	for(uint window_idx = 0; window_idx < windows_size; window_idx++)
	{
    const double2 window = windows[window_idx];
		const int window_flag = window.x >= wl_result && window.y <= wl_result;
		  
		// Jeżeli oba warunki spełnione to mamy dopasowania
		// W przeciwnym przypadku zostawiamy to co było.

		//global_flag = select(global_flag, 1, window_flag == 1);
		if(window_flag == 1) {
		  global_flag = true;
		  break; // test performance
		}
		  
	}
	const uint row_idx = idx / gridDim.x * blockDim.x;
	if(!global_flag || row_idx >= size) {
	  wavelengths_matrix[idx] = CUDART_INF;
	  spectrums_matrix[idx] = CUDART_INF;
	  errors_matrix[idx] = CUDART_INF;;
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
  printf("x %d y %d", blocksPerGrid.x,  blocksPerGrid.y);
  generateWavelengthsMatrix<<<blocksPerGrid, threadsPerBlock>>>(d_output, d_params);
  checkCudaErrors(cudaGetLastError());
  cudaDeviceSynchronize();
  checkCudaErrors(cudaMemcpy(h_output, d_output, matrix_size, cudaMemcpyDeviceToHost));
  cudaFree(d_output);
  cudaFree(d_params);
}

extern "C"
void singleInterpolation(
    double *h_matrix_x,
    double *h_matrix_y,
    size_t *h_sizes_x,
    const size_t size, // number of quasars(coln)
    double *h_matrix_s,
    double *h_matrix_t,
    const size_t size_s,
    double *h_output
) {
  double *d_matrix_x, *d_matrix_y, *d_matrix_s, *d_matrix_t, *d_output;
  size_t *d_sizes_y;
  const size_t x_matrix_size = size * ASTRO_OBJ_SIZE * sizeof(double);
  
  checkCudaErrors(cudaMalloc((void**)&d_output, x_matrix_size));
  checkCudaErrors(cudaMalloc((void**)&d_matrix_x, x_matrix_size));
  checkCudaErrors(cudaMalloc((void**)&d_matrix_y, x_matrix_size));
  
  checkCudaErrors(cudaMalloc((void**)&d_sizes_y, size * sizeof(double)));
  
  checkCudaErrors(cudaMalloc((void**)&d_matrix_s, x_matrix_size));
  checkCudaErrors(cudaMalloc((void**)&d_matrix_t, x_matrix_size));
  
  checkCudaErrors(cudaMemcpy(d_matrix_x, h_matrix_x, x_matrix_size, cudaMemcpyHostToDevice));
  checkCudaErrors(cudaMemcpy(d_matrix_y, h_matrix_y, x_matrix_size, cudaMemcpyHostToDevice));
  checkCudaErrors(cudaMemcpy(d_matrix_s, h_matrix_s, x_matrix_size, cudaMemcpyHostToDevice));
  checkCudaErrors(cudaMemcpy(d_matrix_t, h_matrix_t, x_matrix_size, cudaMemcpyHostToDevice));
  
  const uint threadsPerBlock = 16;
  const uint blocksPerGrid = calculateBlockNumber(size, threadsPerBlock);
  
  single_interpolation<<<blocksPerGrid, threadsPerBlock>>>(d_matrix_x,d_matrix_y,
    d_sizes_y,
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
  cudaFree(d_sizes_y);
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
  // calling kernel
  /*__global__ void filterWithWavelengthWindows
	(
		double *wavelengths_matrix,	// Długości fal widm
		double *spectrums_matrix,	// Widma
		double *errors_matrix,		// Błędy pomiaru widm
		uint  *sizes,		 	// Rozmiary widm w spectrums_matrix
		const double2 *windows,	// Okna	długości fal
		uint windows_size 		// Ilość okien
	)	*/
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