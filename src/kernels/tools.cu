#include "cuda_tools.cc"

__device__
double4 cDouble4(
  const double x,
  const double y,
  const double z,
  const double w
) {
  double4 result;
  result.x = x;
  result.y = y;
  result.z = z;
  result.w = w;
  return result;
}

__device__
double2 cDouble2(
  const double x,
  const double y
) {
  double2 result;
  result.x = x;
  result.y = y;
  return result;
}



template<typename T> 
struct cuda_plus {
  __device__
  T operator()(const T a, const T b) {
    return a + b;
  } 
};

template<typename T> 
struct cuda_minus{
  __device__
  T operator()(const T a, const T b) {
    return a - b;
  } 
};

template <typename F>
__device__
double4 double4Op(const double4 a, const double4 b, F f) {
  double4 result;
  result.x = f(a.x, b.x);
  result.y = f(a.y, b.y);
  result.z = f(a.z, b.z);
  result.w = f(a.w, b.w);
  return result;
}

// Oblicza współczynniki prostej regresjii.
//
//
// Struktura wyniku: (a, b, sia2, sib2, siy2, x_sum, x^2_sum, y_sum)
// done
__global__
void reglin(
  const double * xs,	// Macierz x'ów
  const double * ys,	// Macierz y'ów
  const size_t *cols_sizes, 	// Tablica ilośći znaczących elememntów 
					   	// z kolejnych kolumn macierzy input
  double8 *results,	// Tablica ze wszystkim współczynnikami dla każdego wiersza.
  uint width,
  uint height
	) 
{
	const uint gid = blockIdx.x * blockDim.x + threadIdx.x;		
	if(gid >= width)
	{
		return;
	}

	uint idx = gid;
	const uint col_height = cols_sizes[gid];
	const uint idx_max = idx + col_height * width;

	double x_sum = 0.0;
	double x2_sum = 0.0;
	double y_sum = 0.0;
	double xy_sum = 0.0;

	double n = (double)col_height;
	
	
	double4 c1 = cDouble4(0.0, 0.0, 0.0, 0.0);
	while (idx < idx_max) 
	{
		double4 t, u;
		const double x = xs[idx];
		const double y = ys[idx];

		u = double4Op(cDouble4(x, pow(x, 2), y, x*y), c1, cuda_minus<double>());
		t = double4Op(cDouble4(x_sum, x2_sum, y_sum, xy_sum), u, cuda_plus<double>());
		const double4 tmp = double4Op(t, cDouble4(x_sum, x2_sum, y_sum, xy_sum), cuda_minus<double>());
		c1 = double4Op(tmp, u, cuda_minus<double>());
		
		x_sum 	= t.x;
		x2_sum 	= t.y;
		y_sum	= t.z;
		xy_sum	= t.w;
		
	
		idx += width;
	}

	// x - x_sum
	// y - x2_sum
	// z - y_sum
	// w - xy_sum

	// Obliczenie współczynników a i b
	const double dd = n * x2_sum - pow(x_sum, 2);
	const double a = ((n * xy_sum) - (x_sum * y_sum)) / dd;
	const double b = ((x2_sum * y_sum) - (x_sum * xy_sum)) / dd;

	// OBLICZENIE SUMY POTRZEBNEJ DO OBLICZENI SIY2
	idx = gid;
	double siy2 = 0;

	double c = 0;
	while(idx < idx_max) 
	{
		double x, y, t, u;
		x = xs[idx];
		y = ys[idx];
		
		u = pow(y - b - (a * x), 2) - c;
		t = siy2 + u;
		c = (t - siy2) - u;
		
		siy2 = t;

		idx += width;
	}

	// ZAPIS WYNIKÓW DO PAMIĘCI GLOBALNEJ			
	siy2 = (1.0 / (n - 2.0) ) * siy2;
	
	double8 result;
	result.x1 = a;
	result.x2 = b;
	result.x3 = sqrt(n * siy2 / dd);
	result.x4 = sqrt(siy2 * x2_sum / dd);
	result.x5 = sqrt(siy2);
	result.x6 = x_sum;
	result.x7 = x2_sum;
	result.x8 = y_sum;
	
	results[gid] = result;
}


// Oblicza chisq między odpowiadającymi sobie kolumnami
// z macierzy fs i ys.
//
__global__
void chisq(
  const double *xs,	
	const double *ys,	
	const double *errs,	
  const size_t *cols_sizes, 	
  double *results,
  const size_t width
) 
{
	const uint gid = blockIdx.x * blockDim.x + threadIdx.x;		
	if(gid >= width)
	{
		return;
	}

	uint idx = gid;
	const uint col_height = cols_sizes[gid];
	const uint idx_max = idx + col_height * width;

	double chi2 = 0.0;

	while (idx < idx_max) 
	{
    chi2 += pow((ys[idx] - xs[idx]) / errs[idx], 2);
		idx += width;
	}		
	results[gid] = chi2;
}


__global__
void reglin_simplified(
  const double * xs,	// Macierz x'ów
  const double * ys,	// Macierz y'ów
  const size_t width,
  const size_t height,
  const size_t *cols_sizes, 	// Tablica ilośći znaczących elememntów 
					   	// z kolejnych kolumn macierzy input
  double *results	// Tablica ze wszystkim współczynnikami dla każdego wiersza.
	) 
{
	const uint gid = blockIdx.x * blockDim.x + threadIdx.x;		
	if(gid >= width)
	{
		return;
	}

	uint idx = gid;
	uint col_height = cols_sizes[gid];
	uint idx_max = idx + col_height * width;

	double x2_sum = 0;
	double xy_sum = 0;	
	{
		while (idx < idx_max) 
		{
			double2 t, u;
			const double x = xs[idx];
			const double y = ys[idx];

			u = cDouble2(x * x,  x*y);			
			t = cDouble2(x2_sum + u.x, xy_sum + u.y);
			
			x2_sum 	= t.x;
			xy_sum 	= t.y;

			idx += width;
		}
	}
	double r = (xy_sum / x2_sum);
	if(isnan(r))
	{
		r = 0.0;
	}	
	results[gid] = (double)(r);
}


extern "C"
void calculateReglinSimplified(
  const double *h_x,
  const double *h_y,
  const size_t width,
  const size_t height,
  const size_t *h_cols_sizes,
  double *h_results
) {
  // allocating kernel data
  double *d_x, *d_y, *d_results;
  size_t *d_cols_sizes;
  const size_t matrix_size = width * height * sizeof(double);
  const size_t vector_size_t_size = width *  sizeof(size_t);
  const size_t vector_double_size = width *  sizeof(double);
  checkCudaErrors(cudaMalloc((void**)&d_x, matrix_size));
  checkCudaErrors(cudaMalloc((void**)&d_y, matrix_size));
  checkCudaErrors(cudaMalloc((void**)&d_cols_sizes, vector_size_t_size));
  checkCudaErrors(cudaMalloc((void**)&d_results, vector_double_size));
  // copying data
  checkCudaErrors(cudaMemcpy(d_x, h_x, matrix_size, cudaMemcpyHostToDevice));
  checkCudaErrors(cudaMemcpy(d_y, h_y, matrix_size, cudaMemcpyHostToDevice));
  checkCudaErrors(cudaMemcpy(d_cols_sizes, h_cols_sizes, vector_size_t_size, cudaMemcpyHostToDevice));
  // instatiating kernel
  const size_t threadsPerBlock = BLOCK_DIM;
  const size_t blocksPerGrid = calculateBlockNumber(width, threadsPerBlock);
  // calling kernel
  reglin_simplified<<<blocksPerGrid, threadsPerBlock>>>(d_x, d_y, width, height, d_cols_sizes, d_results);
  cudaDeviceSynchronize();
  checkCudaErrors(cudaGetLastError());
  //copying memory back
  checkCudaErrors(cudaMemcpy(h_results, d_results, vector_double_size, cudaMemcpyDeviceToHost));
  // free memory
  cudaFree(d_x);
  cudaFree(d_y);
  cudaFree(d_cols_sizes);
  cudaFree(d_results);
}


extern "C"
void calculateChisq(
  const double *h_x,
  const double *h_y,
  const double *h_e,
  const size_t *h_col_sizes,
  double *h_results,
  const size_t width,
  const size_t height
) {
  // allocating kernel data
  double *d_x, *d_y, *d_e, *d_results;
  size_t *d_col_sizes;
  const size_t matrix_size = width * height * sizeof(double);
  const size_t vector_size_t_size = width * sizeof(size_t);
  const size_t vector_double_size = width * sizeof(double);
  checkCudaErrors(cudaMalloc((void**)&d_x, matrix_size));
  checkCudaErrors(cudaMalloc((void**)&d_y, matrix_size));
  checkCudaErrors(cudaMalloc((void**)&d_e, matrix_size));
  checkCudaErrors(cudaMalloc((void**)&d_col_sizes, vector_size_t_size));
  checkCudaErrors(cudaMalloc((void**)&d_results, vector_double_size));
  // copying data
  checkCudaErrors(cudaMemcpy(d_x, h_x, matrix_size, cudaMemcpyHostToDevice));
  checkCudaErrors(cudaMemcpy(d_y, h_y, matrix_size, cudaMemcpyHostToDevice));
  checkCudaErrors(cudaMemcpy(d_e, h_e, matrix_size, cudaMemcpyHostToDevice));
  checkCudaErrors(cudaMemcpy(d_col_sizes, h_col_sizes, vector_size_t_size, cudaMemcpyHostToDevice));
  // instatiating kernel
  const size_t threadsPerBlock = BLOCK_DIM;
  const size_t blocksPerGrid = calculateBlockNumber(width, threadsPerBlock);
  // calling kernel
  chisq<<<blocksPerGrid, threadsPerBlock>>>(
    d_x,
    d_y,
    d_e,
    d_col_sizes,
    d_results,
    width
  );
  cudaDeviceSynchronize();
  checkCudaErrors(cudaGetLastError());
  //copying memory back
  checkCudaErrors(cudaMemcpy(h_results, d_results, vector_double_size, cudaMemcpyDeviceToHost));
  // free memory
  cudaFree(d_x);
  cudaFree(d_y);
  cudaFree(d_e);
  cudaFree(d_col_sizes);
  cudaFree(d_results);
}


extern "C"
void calculateReglin(
  const double *h_x,
  const double *h_y,
  const size_t *h_col_sizes,
  double8 *h_reglin_results,
  const size_t width,
  const size_t height
)
{
  // allocating kernel data
  double *d_x, *d_y;
  size_t *d_col_sizes;
  double8 *d_reglin_results;
  const size_t matrix_size = width * height * sizeof(double);
  const size_t result_size = width * sizeof(double8);
  const size_t vector_size = width * sizeof(double);
  checkCudaErrors(cudaMalloc((void**)&d_x, matrix_size));
  checkCudaErrors(cudaMalloc((void**)&d_y, matrix_size));
  checkCudaErrors(cudaMalloc((void**)&d_col_sizes, vector_size));
  checkCudaErrors(cudaMalloc((void**)&d_reglin_results, result_size));
  // copying data
  checkCudaErrors(cudaMemcpy(d_x, h_x, matrix_size, cudaMemcpyHostToDevice));
  checkCudaErrors(cudaMemcpy(d_y, h_y, matrix_size, cudaMemcpyHostToDevice));
  checkCudaErrors(cudaMemcpy(d_col_sizes, h_col_sizes, vector_size, cudaMemcpyHostToDevice));
  // instatiating kernel
  const size_t threadsPerBlock = BLOCK_DIM;
  const size_t blocksPerGrid = calculateBlockNumber(width, threadsPerBlock);
  // calling kernel
  reglin<<<blocksPerGrid, threadsPerBlock>>>(
    d_x,
    d_y,
    d_col_sizes,
    d_reglin_results,
    width,
    height
  );
  cudaDeviceSynchronize();
  checkCudaErrors(cudaGetLastError());
  //copying memory back
  checkCudaErrors(cudaMemcpy(h_reglin_results, d_reglin_results, result_size, cudaMemcpyDeviceToHost));
  // free memory
  cudaFree(d_x);
  cudaFree(d_y);
  cudaFree(d_col_sizes);
  cudaFree(d_reglin_results);
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