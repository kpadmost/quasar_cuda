#include "cuda_tools.cc"
// done
__global__ 
void fix_reglin_results
	(
		double8 *c_reglin_results, //  Parametry prostej regresji liniowej cont
		double8 *reglin_results, // Parametry prostej regresji liniowej
		const size_t size	// Ilość parametrów
	)
{
	// gid0
	const uint gid0 =  blockIdx.x * blockDim.x + threadIdx.x ;

	if(gid0 >= size)
	{
		return;
	}
	
	// r =  (a, b, sia2, sib2, siy2, x_sum, x^2_sum, y_sum)
	double8 r = reglin_results[gid0];	

	// b = 10^b
	r.x2 = pow(10.0, r.x2);
	// sia2 = sia2^0.5
	r.x3 = pow(r.x3, 0.5);
	// sib2 = 10^b * log(10) * sib2^0.5	
	r.x4 = r.x2 * log(10.0) * pow(r.x4, 0.5);
	// siy2 = siy2^0.5
	r.x5 = pow(r.x5, 0.5);

	reglin_results[gid0] = r;

	r = c_reglin_results[gid0];
	r.x2 = pow(10.0, r.x2);
	c_reglin_results[gid0] = r;
}

// done
__global__
void calc_cfun_dcfun
	(
		const double *wavelengths_matrix, 	// Długości fali dla widm [orygniał]
		double *dcfuns_matrix, 	// 
		double *cfuns_matrix, // Kontinuum
		const double8 *reglin_results, 	// Parametry prostej regresji liniowej
		const double8 *c_reglin_results 	// Parametry prostej regresji liniowej
	)
{
	const uint gid0 =  blockIdx.x * blockDim.x + threadIdx.x;
	const uint gid1 =  blockIdx.y * blockDim.y + threadIdx.y;	

	// Obliczenie indeksu elementu
	const uint idx = gid0 * ASTRO_OBJ_SIZE + gid1;

	const uint col_idx = idx % (gridDim.x * blockDim.x);  

	//         (0,    1,    2,    3,    4, ....)
	// cpar =  (a, 10^b, sia2, sib2, siy2, ...)
	double8 cpar = reglin_results[col_idx];
	
	const double wavelength = wavelengths_matrix[idx];

	const double cfun = cpar.x2 * pow(wavelength, cpar.x1);
	cfuns_matrix[idx] = cfun;
		
	//        (0, 1,           2,                         3,        4, ...)
	// par =  (a, 10^b, sia2^0.5, 10^b * log(10) * sib2^0.5, siy2^0.5, ...)
	const double8 par = c_reglin_results[col_idx];

	// dcfun= sib * la^a + sia * b * a * log(lambda)
	double dcfun = par.x4 * pow(wavelength, par.x1);
	dcfun += par.x3 * par.x2 * par.x1 * log(wavelength);
	dcfuns_matrix[idx] = dcfun;
}

//done
__global__ 
void calculate_cfun // calculate continuum baseline
	(
		const double *wavelengths_matrix_filtered, 	// Długości fali dla widm po filtracji
		double *cfuns_filtered,	// Kontinuum dla elementów po filtracji (w oknach)
		const uint filtered_size,			// Maksymalną ilość znaczących elementów po filtracji
		const double8 *c_reglin_results 	// Parametry prostej regresji liniowej
	)
{
	const uint gid0 = blockIdx.x * blockDim.x + threadIdx.x;
	const uint gid1 = blockIdx.y * blockDim.y + threadIdx.y;	

	// Obliczenie indeksu elementu
	const uint idx = gid0 * filtered_size + gid1;
	const uint col_idx = idx % (gridDim.x * blockDim.x);  

	//         (0,    1,    2,    3, ....)
	// cpar =  (a, 10^b, sia2, sib2, ...)
	const double8 cpar = c_reglin_results[col_idx];
	
	const double wavelength = wavelengths_matrix_filtered[idx];

	double cfun_filtered = cpar.x2 * pow(wavelength, cpar.x1);
	cfuns_filtered[idx] = cfun_filtered;			
}

// done
__global__
void reduce_continuum_chisqs
	(
		double *chisqs, 		// Bufor z chisq
		const size_t * filtered_sizes, 	// Ilość znaczących elementów po filtracji
		const size_t size				// Ilość chisq i filtered_sizes
	)
{
	// gid0 - numer elementu z chisqs (czyli jednego chisq)
	uint gid0 = blockIdx.x * blockDim.x + threadIdx.x;
	if(gid0 >= size)
	{
		return;
	}
	chisqs[gid0] = chisqs[gid0] / (double)filtered_sizes[gid0];	
}


extern "C"
void calculateDcfun(
  const double *h_wavelengths,
  double *h_dcfuns,
  double *h_cfuns,
  const double8 *h_reglin, 	// Parametry prostej regresji liniowej
	const double8 *h_creglin 	// Parametry prostej regresji liniowej
  const size_t width,
  const size_t height
)
{
  // allocating kernel data
  double *d_wavelengths, *d_dcfuns, *d_cfuns;
  double8 *d_reglin, d_creglin;
  const size_t matrix_size = width * height * sizeof(double);
  const size_t vector_size = width * sizeof(double8);
  checkCudaErrors(cudaMalloc((void**)&d_wavelengths, matrix_size));
  checkCudaErrors(cudaMalloc((void**)&d_dcfuns, matrix_size));
  checkCudaErrors(cudaMalloc((void**)&d_cfuns, matrix_size));
  checkCudaErrors(cudaMalloc((void**)&d_reglin, vector_size));
  checkCudaErrors(cudaMalloc((void**)&d_input, vector_size));
  // copying data
  checkCudaErrors(cudaMemcpy(d_wavelengths, h_wavelengths, matrix_size, cudaMemcpyHostToDevice));
  checkCudaErrors(cudaMemcpy(d_reglin, h_reglin, vector_size, cudaMemcpyHostToDevice));
  checkCudaErrors(cudaMemcpy(d_creglin, h_creglin, vector_size, cudaMemcpyHostToDevice));
  // instatiating kernel
  const dim3 threadsPerBlock(1, BLOCK_DIM, 1);
  const dim3 blocksPerGrid(width / threadsPerBlock.x, height / threadsPerBlock.y, 1);
  // calling kernel
  calc_cfun_dcfun<<<blocksPerGrid, threadsPerBlock>>>(
    d_wavelengths,
    d_dcfuns,
    d_cfuns,
    d_reglin,
    d_creglin
  );
	
  cudaDeviceSynchronize();
  checkCudaErrors(cudaGetLastError());
  //copying memory back
  checkCudaErrors(cudaMemcpy(h_dcfuns, d_dcfuns, matrix_size, cudaMemcpyDeviceToHost));
  checkCudaErrors(cudaMemcpy(h_cfuns, d_cfuns, matrix_size, cudaMemcpyDeviceToHost));
  // free memory
  cudaFree(d_wavelengths);
  cudaFree(d_dcfuns);
  cudaFree(d_cfuns);
  cudaFree(d_reglin);
  cudaFree(d_creglin);
}


extern "C"
void reduceChisqs(
  double *h_chisq,
  const size_t *h_sizes,
  const size_t width
)
{
  // allocating kernel data
  double *d_chisq;
  size_t *d_sizes;
  const size_t vector_double_size = width * sizeof(double);
  const size_t vector_size_t_size = width * sizeof(size_t);
  checkCudaErrors(cudaMalloc((void**)&d_chisq, vector_double_size));
  checkCudaErrors(cudaMalloc((void**)&d_sizes, vector_size_t_size));
  // copying data
  checkCudaErrors(cudaMemcpy(d_chisq, h_chisq, vector_double_size, cudaMemcpyHostToDevice));
  checkCudaErrors(cudaMemcpy(d_sizes, h_sizes, vector_size_t_size, cudaMemcpyHostToDevice));
  // instatiating kernel
  const size_t threadsPerBlock = BLOCK_DIM;
  const size_t blocksPerGrid = calculateBlockNumber(width, threadsPerBlock);
  // calling kernel
  reduce_continuum_chisqs<<<blocksPerGrid, threadsPerBlock>>>(d_chisq, d_sizes, width);
  cudaDeviceSynchronize();
  checkCudaErrors(cudaGetLastError());
  //copying memory back
  checkCudaErrors(cudaMemcpy(h_chisq, d_chisq, vector_double_size, cudaMemcpyDeviceToHost));
  // free memory
  cudaFree(d_chisq);
  cudaFree(d_sizes);
}


extern "C"
void fixContinuumReglin(
  double8 *h_continuum_reglin,
  double8 *h_reglin,
  const size_t size
)
{
  // allocating kernel data
  double8 *d_reglin, *d_continuum_reglin;
  const size_t vector_size = size * sizeof(double8);
  checkCudaErrors(cudaMalloc((void**)&d_reglin, vector_size));
  checkCudaErrors(cudaMalloc((void**)&d_continuum_reglin, vector_size));
  // copying data
  checkCudaErrors(cudaMemcpy(d_reglin, h_reglin, vector_size, cudaMemcpyHostToDevice));
  checkCudaErrors(cudaMemcpy(d_continuum_reglin, h_continuum_reglin, vector_size, cudaMemcpyHostToDevice));
  // instatiating kernel
  const int threadsPerBlock = BLOCK_DIM;
  const int blocksPerGrid = calculateBlockNumber(size, BLOCK_DIM);
  // calling kernel
  fix_reglin_results<<<blocksPerGrid, threadsPerBlock>>>(d_continuum_reglin, d_reglin, size);
  cudaDeviceSynchronize();
  checkCudaErrors(cudaGetLastError());
  //copying memory back
  checkCudaErrors(cudaMemcpy(d_reglin, h_reglin, vector_size, cudaMemcpyDeviceToHost));
  checkCudaErrors(cudaMemcpy(d_continuum_reglin, h_continuum_reglin, vector_size, cudaMemcpyDeviceToHost));
  // free memory
  cudaFree(d_reglin);
  cudaFree(d_continuum_reglin);
}

extern "C"
void calculateContinuumFunction(
  const double *h_wavelengths,
  double *h_cfun,
  const double8 *h_reglin_vector,
  const size_t width,
  const size_t height
)
{
  // allocating kernel data
  double *d_wavelengths, *d_cfun;
  double8 *d_reglin_vector;
  const size_t matrix_size = width * height * sizeof(double);
  const size_t vector_size = width * sizeof(double8);
  checkCudaErrors(cudaMalloc((void**)&d_wavelengths, matrix_size));
  checkCudaErrors(cudaMalloc((void**)&d_cfun, matrix_size));
  checkCudaErrors(cudaMalloc((void**)&d_reglin_vector, vector_size));
  // copying data
  checkCudaErrors(cudaMemcpy(d_wavelengths, h_wavelengths, matrix_size, cudaMemcpyHostToDevice));
  checkCudaErrors(cudaMemcpy(d_cfun, h_cfun, matrix_size, cudaMemcpyHostToDevice));
  
  checkCudaErrors(cudaMemcpy(d_reglin_vector, h_reglin_vector, vector_size, cudaMemcpyHostToDevice));
  // instatiating kernel
  const dim3 threadsPerBlock(1, BLOCK_DIM, 1);
  const dim3 blocksPerGrid(width / threadsPerBlock.x, calculateBlockNumber(height, threadsPerBlock.y) / threadsPerBlock.y, 1);
  // calling kernel
  calculate_cfun<<<blocksPerGrid, threadsPerBlock>>>(
    d_wavelengths,
    d_cfun,
    height,
    d_reglin_vector
  );
  cudaDeviceSynchronize();
  checkCudaErrors(cudaGetLastError());
  //copying memory back
  checkCudaErrors(cudaMemcpy(h_cfun, d_cfun, matrix_size, cudaMemcpyDeviceToHost));
  // free memory
  cudaFree(d_wavelengths);
  cudaFree(d_cfun);
  cudaFree(d_reglin_vector);
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