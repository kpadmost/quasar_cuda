#include "cuda_tools.cc"

// generate code

__global__
void generateWavelengthsMatrix
	(
		double * output_matrix
		const double4* params_vector,
	)
{
	// gid0 - numer widma kwazaru
	uint gid0 = get_global_id(0);	
	
	// parametry a, b i z kwazaru
	double4 abz_;
	__local double4 local_abz_;

	if (get_local_id(0) == 0)
	{
  		local_abz_ = abz_buffer[gid0];
	}
	barrier(CLK_LOCAL_MEM_FENCE);
	abz_ = local_abz_;
	
	// gid1 - numer elementu widma (indeks od 0 do 4095)
	uint gid1 = get_global_id(1);
	uint idx = gid0 * get_global_size(1) + gid1;
	
	// Wyliczenie lambdy dla tego gid1
	double wl = (abz_.x * (double)(gid1)) + abz_.y;
	wl = pow((double)(10),wl);
	
	// Uwzglednienie przesuniecia kosmologicznego - przejscie do ukladu emitujacego z wykorzystaniem redshiftu
	wavelengths_matrix[idx] = wl / (abz_.z + (double)(1));
	
	return;
}


extern "C"
void generateWavelengths(
  double* h_output,
  double4* params,
  const size_t spectrum_number
) 
{
  double* d_output = 0, d_params;
  checkCudaErrors(cudaSetDevice(0));
  checkCudaErrors(cudaMalloc(d_output), spectrum_number * ASTRO_OBJ_SIZE * sizeof(double));
  checkCudaErrors(cudaMalloc(d_params), spectrum_size * sizeof(double4));
  
}