#include "cuda_tools.cc"

__global__ 
void fix_reglin_results
	(
		double8 *c_reglin_results, //  Parametry prostej regresji liniowej cont
		double8 *reglin_results, // Parametry prostej regresji liniowej
		uint size	// Ilość parametrów
	)
{
	// gid0
	uint gid0 =  blockIdx.x * blockDim.x + threadIdx.x ;

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

__global__
void calc_cfun_dcfun
	(
		const double *wavelengths_matrix, 	// Długości fali dla widm [orygniał]
		const double *dcfuns_matrix, 	// 
		const double *cfuns_matrix, // Kontinuum
		double8 *reglin_results, 	// Parametry prostej regresji liniowej
		double8 *c_reglin_results 	// Parametry prostej regresji liniowej
	)
{
	const uint gid0 =  blockIdx.x * blockDim.x + threadIdx.x;
	const uint gid1 =  blockIdx.y * blockDim.y + threadIdx.y;	

	// Obliczenie indeksu elementu
	const uint idx = gid0 * ASTRO_OBJ_SPEC_SIZE + gid1;

	const uint col_idx = idx % (gridDim.x * blockDim.x);  

	//         (0,    1,    2,    3,    4, ....)
	// cpar =  (a, 10^b, sia2, sib2, siy2, ...)
	double8 cpar = reglin_results[col_idx];
	
	double wavelength = wavelengths_matrix[idx];

	double cfun = cpar.s1 * pow(wavelength, cpar.s0);
	cfuns_matrix[idx] = cfun;
		
	//        (0, 1,           2,                         3,        4, ...)
	// par =  (a, 10^b, sia2^0.5, 10^b * log(10) * sib2^0.5, siy2^0.5, ...)
	double8 par = c_reglin_results[col_idx];

	// dcfun= sib * la^a + sia * b * a * log(lambda)
	double dcfun = par.s3 * pow(wavelength, par.s0);
	dcfun += par.s2 * par.s1 * par.s0 * log(wavelength);
	dcfuns_matrix[idx] = dcfun;
}

__kernel void calc_cw
	(
		__global double * wavelengths_matrix_filtered, 	// Długości fali dla widm po filtracji
		__global double * cfuns_filtered,	// Kontinuum dla elementów po filtracji (w oknach)
		uint filtered_size,			// Maksymalną ilość znaczących elementów po filtracji
		__global double8 * c_reglin_results 	// Parametry prostej regresji liniowej
	)
{
	uint gid0 = get_global_id(0);
	uint gid1 = get_global_id(1);	

	// Obliczenie indeksu elementu
	uint idx = gid0 * filtered_size + gid1;
	uint col_idx = idx % get_global_size(0);  

	//         (0,    1,    2,    3, ....)
	// cpar =  (a, 10^b, sia2, sib2, ...)
	double4 cpar = c_reglin_results[col_idx].lo;
	
	double wavelength = wavelengths_matrix_filtered[idx];

	double cfun_filtered = cpar.s1 * pow(wavelength, cpar.s0);
	cfuns_filtered[idx] = cfun_filtered;			
}

__kernel void reduce_continuum_chisqs
	(
		__global double * chisqs, 		// Bufor z chisq
		__constant uint * filtered_sizes, 	// Ilość znaczących elementów po filtracji
		uint size				// Ilość chisq i filtered_sizes
	)
{
	// gid0 - numer elementu z chisqs (czyli jednego chisq)
	uint gid0 = get_global_id(0);

	if(gid0 >= size)
	{
		return;
	}
		
	double filtered_size, chisq;

	filtered_size = (double)filtered_sizes[gid0];
	chisq = chisqs[gid0];

	chisq /= filtered_size;

	chisqs[gid0] = chisq;	
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