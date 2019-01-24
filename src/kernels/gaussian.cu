#include "cuda_tools.h"
#include <cfloat>
#include <vector>
__device__
inline double4 createDouble4(
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
inline double max_in_double4(const double4 val)
{
	return max(max(val.x, val.y), max(val.z, val.w));
}

// Gaussian
__device__
inline double f(const double x, const double a, const double b, const double c)
{
	return a * exp( -0.5 * pow(x - b, 2) / pow(c, 2));
}

// dfda
__device__
inline double dfda(double x, double b, const double c)
{
	return exp( -0.5 * pow(x - b, 2) / pow(c, 2));
}

// dfdb
__device__
double dfdb(const double x, const double a, const double b, const double c)
{
	return ((a * (x - b)) / pow(c, 2)) * exp( -0.5 * pow(x - b, 2) / pow(c, 2));
}

// dfdc
__device__
double dfdc(const double x, const double a, const double b, const double c)
{
	return ((a * pow(x - b, 2)) / pow(c, 3)) * exp( -0.5 * pow(x - b, 2) / pow(c, 2));
}

//
// Dopasowuje do podanych wartości (y'ów i x'ów) gaussiana (funkcje Gaussa),
// tj. znajduję współczynniki a, b i c funkcji 
//
// f(x) = a * exp( - 0.5 * (x - b)^2 / c^2)
//
// , dla których funkcja najlepiej przybliża zadane wartości.
//
// Dane są ułożone kolumnami w macierzach ys i xs.
//
// Kernel korzysta z metody Levenberg–Marquardt dla nieliniowego problemu
// najmniejszych kwadratów. 
// Wikipedia: http://en.wikipedia.org/wiki/Levenberg%E2%80%93Marquardt_algorithm
//
// Potrzebne definicje
#define MAX_ALFA 1.0e+10
#define MIN_ALFA 1.0e-7
#define ALFA_CHANGE 10
//
__global__
void fit_gaussian(
	const double * xs,	// Argumenty (x'y) dla których dopasowujemy
	const double * ys,	// Wartości, które dopasowujemy
	const size_t *cols_sizes, 	// Tablica ilośći znaczących elememntów 
	double4 * results,	// Wyniki dopasowania,
	const size_t width,
				   	// z kolejnych kolumn macierzy ys
	const int max_iterations
					// na początku znajdują się tam początkowe
					// wartości dla współczynników a, b i c.
	)
{
	const uint gid = blockIdx.x * blockDim.x + threadIdx.x;	
	if(gid >= width)
	{
		return;
	}

	// W metodzie GN iteracyjnie dfrom vectorochodzimy do "najlepszych"
	// współczynników rozwiązując za każdym razem równanie
	// (J^T * J + alfa * diag(JT * J)) * dB = J^T * R  dla dB
	//
	// (Można rozwiązać równienie dB = (J^T * J + alfa * diag(JT * J))^-1 * J^T * R
	// ale odwracanie macierzy jest problematyczne, a tu są tylko
	// 3 parametry więc łatwiej użyć wzorów Cramera.)
	//
	// B - współczynniki f(x)
	// dB - zmiana wartości współczynników.
	// R - różnice y - f(x, B) 
	// J - jakobian dla funkcji f(x)

	// J^T * J
	double JTJ[3][3];
	// Współczynniki
 	
	// J^T * R	
	double JTR[3]; 
	double4 B = results[gid];
	// Flaga, że dopasowanie nie doszło do skutku
	B.w = 0.0;

	// Indeks
	uint idx = gid;
	const size_t col_height = cols_sizes[gid];
	if(col_height == 0)
	{
		return;
	}
	const uint idx_max = idx + col_height * width;

	//
	double alfa = 100;

	// Dopasowanie dla aktualnych współczynników
	double chisq_B = 0.0;
	// Dopasowanie dla nowych współczynników
	double chisq_newB = 0.0;
	
	for(uint i = 0; i < max_iterations; i++)
	{		
		// Zerowanie danych
		JTJ[0][0] = 0;
		JTJ[0][1] = 0;
		JTJ[0][2] = 0;
		JTJ[1][0] = 0;
		JTJ[1][1] = 0;
		JTJ[1][2] = 0;
		JTJ[2][0] = 0;
		JTJ[2][1] = 0;
		JTJ[2][2] = 0;

		JTR[0] = 0;
		JTR[1] = 0;
		JTR[2] = 0;

		chisq_B = 0.0;

		idx = gid;
		while(idx < idx_max)
		{
			double x = xs[idx];

			// dfda
			double dfda_x = dfda(x, B.y, B.z);
			//jacobian[idx] = dfd_;
			// dfdb
			double dfdb_x = dfdb(x, B.x, B.y, B.z);
			//jacobian[idx + width] = dfd_;
			// dfdc
			double dfdc_x = dfdc(x, B.x, B.y, B.z);
			//jacobian[idx + 2 * width] = dfd_

			double y = ys[idx];
		
  		JTJ[0][0] += dfda_x * dfda_x;
  		JTJ[0][1] += dfda_x * dfdb_x;
  		JTJ[0][2] += dfda_x * dfdc_x;

  		JTJ[1][0] += dfdb_x * dfda_x;
  		JTJ[1][1] += dfdb_x * dfdb_x;
  		JTJ[1][2] += dfdb_x * dfdc_x;

  		JTJ[2][0] += dfdc_x * dfda_x;
  		JTJ[2][1] += dfdc_x * dfdb_x;
  		JTJ[2][2] += dfdc_x * dfdc_x;

       			// R[idx]
			double r = y - f(x, B.x, B.y, B.z);
			// chisq_B
			chisq_B += pow(r, 2);	
			//JT * R
  		JTR[0] += dfda_x * r;
  		JTR[1] += dfdb_x * r;
  		JTR[2] += dfdc_x * r;   

			idx += width;  
		}
		chisq_B /= 2.0;

		double diagJTJ[3];
		diagJTJ[0] = JTJ[0][0];
		diagJTJ[1] = JTJ[1][1];
		diagJTJ[2] = JTJ[2][2];

		// Metoda największego spadku w LM
		// (modyfikacja względem algorytmu Gaussa Newtona)
		JTJ[0][0] += alfa * diagJTJ[0];
		JTJ[1][1] += alfa * diagJTJ[1];
		JTJ[2][2] += alfa * diagJTJ[2];

    		// (JT * J + alfa * diag(JT * J) ) * dB = JT * R jest równaniem typu Ax = b
    		// A = (JT * J), dB = x, JT * R = b
    		// Rozwiązanie za pomocą wzorów Cramera, wikipedia: http://en.wikipedia.org/wiki/Cramer%27s_rule
    		// x_i = det(A_i)/det(A)

		double detA = 
   		JTJ[0][0] * (JTJ[1][1] * JTJ[2][2] - JTJ[1][2] * JTJ[2][1]) -
    		JTJ[0][1] * (JTJ[1][0] * JTJ[2][2] - JTJ[1][2] * JTJ[2][0]) + 
    		JTJ[0][2] * (JTJ[1][0] * JTJ[2][1] - JTJ[1][1] * JTJ[2][0]) ;

		double detA1 =
    		JTR[0]	  * (JTJ[1][1] * JTJ[2][2] - JTJ[1][2] * JTJ[2][1]) -
    		JTJ[0][1] * (  JTR[1]  * JTJ[2][2] - JTJ[1][2] * JTR[2]   ) + 
    		JTJ[0][2] * (  JTR[1]  * JTJ[2][1] - JTJ[1][1] * JTR[2]   ) ;

		double detA2 = 
    		JTJ[0][0] * (JTR[1]    * JTJ[2][2] - JTJ[1][2] * JTR[2]   ) -
    		JTR[0]	  * (JTJ[1][0] * JTJ[2][2] - JTJ[1][2] * JTJ[2][0]) + 
    		JTJ[0][2] * (JTJ[1][0] * JTR[2]    - JTR[1]    * JTJ[2][0]) ;

		double detA3 = 
    		JTJ[0][0] * (JTJ[1][1] * JTR[2]    - JTR[1]    * JTJ[2][1]) -
    		JTJ[0][1] * (JTJ[1][0] * JTR[2]    - JTR[1]    * JTJ[2][0]) + 
    		JTR[0]	  * (JTJ[1][0] * JTJ[2][1] - JTJ[1][1] * JTJ[2][0]) ;

		if(fabs(detA) < DBL_MIN)
		{			
			break;
		}				

		// Zmiana i sprawdzenie warunków stopu
		{    		
			double4 dB = createDouble4(detA1/detA, detA2/detA, detA3/detA, 0.0);
			
			// B(k+1) = B(k) + dB	
			double4 newB = double4Op(B, dB, cuda_plus<double>());

					
			// Obliczenie dopasowanie dla nowych współczynników
			// jeżeli następuje pierwsza zmiana.
			chisq_newB = 0.0;						
			idx = gid;			
			while(idx < idx_max)
			{
				const double x = xs[idx];
				const double fx = f(x, newB.x, newB.y, newB.z);
				const double y = ys[idx];
			
				// 
				chisq_newB += pow(y - fx, 2);
				idx += width;  
			}
			chisq_newB /= 2.0;
			
			// Sprawdzenie, czy nowe współczynniki są lepsze					
			if(chisq_newB < chisq_B)
			{
				// B(k+1) = B(k)+ dB	
    				B = newB;

				// Modyfikacja w stronę metody Gaussa-Newtonwa
				alfa = max(alfa/ALFA_CHANGE, MIN_ALFA);	
			}
			else
			{
				// Zwiększamy udział metody największego spadku
				// aż dojdziemy do maksymalnego wpływu.
				while(alfa != MAX_ALFA && i < max_iterations)
				{
					i++;

					// Modyfikacja w stronę metody największego spadku
					alfa = min(alfa*ALFA_CHANGE, MAX_ALFA);	

					// Metoda największego spadku w LM
					// (modyfikacja względem algorytmu Gaussa Newtona)
					JTJ[0][0] += (alfa - alfa/ALFA_CHANGE) * diagJTJ[0];
					JTJ[1][1] += (alfa - alfa/ALFA_CHANGE) * diagJTJ[1];
					JTJ[2][2] += (alfa - alfa/ALFA_CHANGE) * diagJTJ[2];

					detA = 
   					JTJ[0][0] * (JTJ[1][1] * JTJ[2][2] - JTJ[1][2] * JTJ[2][1]) -
    					JTJ[0][1] * (JTJ[1][0] * JTJ[2][2] - JTJ[1][2] * JTJ[2][0]) + 
    					JTJ[0][2] * (JTJ[1][0] * JTJ[2][1] - JTJ[1][1] * JTJ[2][0]) ;

					detA1 =
    					JTR[0]	  * (JTJ[1][1] * JTJ[2][2] - JTJ[1][2] * JTJ[2][1]) -
    					JTJ[0][1] * (  JTR[1]  * JTJ[2][2] - JTJ[1][2] * JTR[2]   ) + 
    					JTJ[0][2] * (  JTR[1]  * JTJ[2][1] - JTJ[1][1] * JTR[2]   ) ;

					detA2 = 
    					JTJ[0][0] * (JTR[1]    * JTJ[2][2] - JTJ[1][2] * JTR[2]   ) -
    					JTR[0]	  * (JTJ[1][0] * JTJ[2][2] - JTJ[1][2] * JTJ[2][0]) + 
    					JTJ[0][2] * (JTJ[1][0] * JTR[2]    - JTR[1]    * JTJ[2][0]) ;

					detA3 = 
    					JTJ[0][0] * (JTJ[1][1] * JTR[2]    - JTR[1]    * JTJ[2][1]) -
    					JTJ[0][1] * (JTJ[1][0] * JTR[2]    - JTR[1]    * JTJ[2][0]) + 
    					JTR[0]	  * (JTJ[1][0] * JTJ[2][1] - JTJ[1][1] * JTJ[2][0]) ;

					if(fabs(detA) < DBL_MIN)
					{			
						break;
					}	

					dB = createDouble4(detA1/detA, detA2/detA, detA3/detA, 0.0);
			
					// B(k+1) = B(k) + dB	
					newB = double4Op(B, dB, cuda_plus<double>());

					
					// Obliczenie dopasowanie dla nowych współczynników
					// jeżeli następuje pierwsza zmiana.
					chisq_newB = 0.0;						
					idx = gid;			
					while(idx < idx_max)
					{
						const double x = xs[idx];
						const double fx = f(x, newB.x, newB.y, newB.z);
						const double y = ys[idx];
			
						// 
						chisq_newB += pow(y - fx, 2);
						idx += width;  
					}
					chisq_newB /= 2.0;	

					if(chisq_newB < chisq_B)
					{
						// B(k+1) = B(k)+ dB	
						B = newB;

						// Modyfikacja w stronę metody Gaussa-Newtonwa
						alfa = max(alfa/ALFA_CHANGE, MIN_ALFA);	
						break;
					}					
				}
				
				// Nie udało się osiągnąć lepszego wyniku dla
				// największego alfa, więc kończymy cały algorytm.
				if(alfa == MAX_ALFA)
				{
					break;
				}									
			}
			
		}
	};

	// Zapisanie flagi, ze dopasowanie doszło do skutku.
	B.w = 1.0;
	// Zapisanie wyników	
	//results[gid] = convert_double4(B);
	results[gid] = B;
}

__global__
void calc_fwhm
	(
		const double4 *fitParams,	// Współczynniki funkcji Gaussa
		double *fwhms,	// 
		const uint width
	)
{
	const uint idx = blockDim.x * blockIdx.x + threadIdx.x;		
	if(idx >= width)
	{
		return;
	}
	double4 gaussParams = fitParams[idx];

	double b = gaussParams.y;
	double c = gaussParams.z;

	// przeliczenie jednoski współczynnika c (sigma): A -> km/s
	const double tmp = pow(1.0 + c/b, 2);
 	const double c_kms = C * (tmp  - 1.0) / (tmp + 1.0) * 1.0e-3;
	fwhms[idx] = c_kms * 2.0 * sqrt(2.0 * log(2.0));
}

//
// Oblicza funkcje Gaussa, gdzie xs jest macierzą z argumentami funkcji.
//


__global__
void calc_gaussian_chisq
	(
		const double *xs,	// Argumenty dla których obliczamy gaussiana
		const double *ys,	// Argumenty dla których obliczamy gaussiana
		const double *errors,
		const double4 *gparams,	// Współczynniki funkcji Gaussa
		const uint *cols_sizes,
		double *chisqs,
		const uint width
	)
{
	const uint gid = blockDim.x * blockIdx.x + threadIdx.x;
	if(gid >= width)
	{
		return;
	}

	// Ilość elementów w kolumnie
	const uint col_height = cols_sizes[gid];

	// Indeks
	uint idx = gid;	
	const uint idx_max = gid + col_height * width;

	const double4 abc = gparams[gid];
	
	double chisq = 0.0;
	double c = 0.0;
	while (idx < idx_max) 
	{
		double f, y, e, t, u, x;
		y = ys[idx];
		x = xs[idx];
		f = abc.x * exp(-0.5 * pow((x - abc.y) / abc.z, 2));
		e = errors[idx];

		u = pow((y-f) / e, 2) - c;
		t = chisq + u;
		c = (t - chisq) - u;

		chisq = t;

		idx += width;
	}

	// Zapis
	chisqs[gid] = chisq;
}


extern "C"
void calculateGaussianChisq(
  const double *h_x,
  const double *h_y,
  const double *h_e,
  const double4 *h_gparams,
  const uint *h_col_sizes,
  double *h_chisq,
  const uint width,
  const uint height
) {
  // allocating kernel data
  double *d_x, *d_y, *d_e, *d_chisq;
  double4 *d_gparams;
  uint *d_col_sizes;
  const size_t matrix_size = width * height * sizeof(double);
  const size_t double_vector_size = width * sizeof(double);
  const size_t double4_vector_size = width  * sizeof(double4);
  const size_t uint_vector_size = width * sizeof(uint);
  checkCudaErrors(cudaMalloc((void**)&d_x, matrix_size));
  checkCudaErrors(cudaMalloc((void**)&d_y, matrix_size));
  checkCudaErrors(cudaMalloc((void**)&d_e, matrix_size));
  checkCudaErrors(cudaMalloc((void**)&d_chisq, double_vector_size));
  checkCudaErrors(cudaMalloc((void**)&d_gparams, double4_vector_size));
  checkCudaErrors(cudaMalloc((void**)&d_col_sizes, uint_vector_size));
  // copying data
  checkCudaErrors(cudaMemcpy(d_x, h_x, matrix_size, cudaMemcpyHostToDevice));
  checkCudaErrors(cudaMemcpy(d_y, h_y, matrix_size, cudaMemcpyHostToDevice));
  checkCudaErrors(cudaMemcpy(d_e, h_e, matrix_size, cudaMemcpyHostToDevice));
  checkCudaErrors(cudaMemcpy(d_gparams, h_gparams, double4_vector_size, cudaMemcpyHostToDevice));
  checkCudaErrors(cudaMemcpy(d_col_sizes, h_col_sizes, uint_vector_size, cudaMemcpyHostToDevice));
  // instatiating kernel
  const size_t threadsPerBlock = BLOCK_DIM;
  const size_t blocksPerGrid = calculateBlockNumber(width, threadsPerBlock);
  // calling kernel
  calc_gaussian_chisq<<<blocksPerGrid, threadsPerBlock>>>(
    d_x,
    d_y,
    d_e,
    d_gparams,
    d_col_sizes,
    d_chisq,
    width
  );
  cudaDeviceSynchronize();
  checkCudaErrors(cudaGetLastError());
  //copying memory back
  checkCudaErrors(cudaMemcpy(h_chisq, d_chisq, double_vector_size, cudaMemcpyDeviceToHost));
  // free memory
  cudaFree(d_x);
  cudaFree(d_y);
  cudaFree(d_e);
  cudaFree(d_gparams);
  cudaFree(d_col_sizes);
  cudaFree(d_chisq);
}

__global__
void calc_gaussian(
		const double *xs,	// Argumenty dla których obliczamy gaussiana
		double *ys,	// Wyniki funkcji
		const double4 *gparams,	// Współczynniki funkcji Gaussa
		const uint *cols_sizes,
		const uint width		
	)
{
	const uint gid0 = blockDim.x * blockIdx.x + threadIdx.x;	
	const uint gid1 = blockDim.y * blockIdx.y + threadIdx.y;	
	const uint col_height = cols_sizes[gid0];
	if(gid0 >= width || gid1 >= col_height)
	{
	  
		return;
	}
	// Indeks
	const uint idx = gid0 + width * gid1;	
	
	// Pobranie x
	const double x = xs[idx];

	// Pobranie parametrów dla tej kolumny x'ów

	
	const double4 abc = gparams[gid0];
	// Zapis
	ys[idx] = abc.x * exp(-0.5 * pow((x - abc.y) / abc.z, 2));
}

extern "C"
void calculateGaussian(
  const double *h_x,
  double *h_y,
  const double4 *h_gaussian,
  const uint *h_sizes,
  const uint width,
  const uint height
) {
  // allocating kernel data
  double *d_x, *d_y;
  double4 *d_gaussian;
  uint *d_sizes;
  
  const size_t matrix_size = width * height * sizeof(double);
  const size_t double4_vector_size = width * sizeof(double4);
  const size_t uint_vector_size = width * sizeof(uint);
  checkCudaErrors(cudaMalloc((void**)&d_x, matrix_size));
  checkCudaErrors(cudaMalloc((void**)&d_y, matrix_size));
  checkCudaErrors(cudaMalloc((void**)&d_gaussian, double4_vector_size));
  checkCudaErrors(cudaMalloc((void**)&d_sizes, uint_vector_size));
  // copying data
  checkCudaErrors(cudaMemcpy(d_x, h_x, matrix_size, cudaMemcpyHostToDevice));
  checkCudaErrors(cudaMemcpy(d_gaussian, h_gaussian, double4_vector_size, cudaMemcpyHostToDevice));
  checkCudaErrors(cudaMemcpy(d_sizes, h_sizes, uint_vector_size, cudaMemcpyHostToDevice));
  // instatiating kernel
  const dim3 threadsPerBlock(1, BLOCK_DIM, 1);
  const dim3 blocksPerGrid(width, calculateBlockNumber(height, threadsPerBlock.y), 1);
  // calling kernel
  calc_gaussian<<<blocksPerGrid, threadsPerBlock>>>(d_x, d_y, d_gaussian, d_sizes, width);
  cudaDeviceSynchronize();
  checkCudaErrors(cudaGetLastError());
  //copying memory back
  checkCudaErrors(cudaMemcpy(h_y, d_y, matrix_size, cudaMemcpyDeviceToHost));
  // free memory
  cudaFree(d_x);
  cudaFree(d_y);
  cudaFree(d_gaussian);
  cudaFree(d_sizes);
}
extern "C"
void calculateFwhm(
  const double4 *h_fit_params,
  double *h_fwhm,
  const uint width
) {
  // allocating kernel data
  double4 *d_fit_params;
  double *d_fwhm;
  const size_t double_vector_size = width * sizeof(double);
  const size_t double4_vector_size = width * sizeof(double4);
  checkCudaErrors(cudaMalloc((void**)&d_fit_params, double4_vector_size));
  checkCudaErrors(cudaMalloc((void**)&d_fwhm, double_vector_size));
  // copying data
  checkCudaErrors(cudaMemcpy(d_fit_params, h_fit_params, double4_vector_size, cudaMemcpyHostToDevice));
  // instatiating kernel
  const uint threadsPerBlock = BLOCK_DIM;
  const uint blocksPerGrid = calculateBlockNumber(width, threadsPerBlock);
  // calling kernel
  calc_fwhm<<<blocksPerGrid, threadsPerBlock>>>(d_fit_params, d_fwhm, width);
  cudaDeviceSynchronize();
  checkCudaErrors(cudaGetLastError());
  //copying memory back
  checkCudaErrors(cudaMemcpy(h_fwhm, d_fwhm, double_vector_size, cudaMemcpyDeviceToHost));
  // free memory
  cudaFree(d_fit_params);
  cudaFree(d_fwhm);
}

extern "C"
void fitGaussian(
  const double *h_x,
  const double *h_y,
  const size_t *h_col_sizes,
  double4 *h_results,
  const size_t width,
  const size_t height,
  const int max_iterations
) {
  // allocating kernel data
  double *d_x, *d_y;
  double4 *d_results;
  size_t *d_col_sizes;
  const size_t matrix_size = width * height * sizeof(double);
  const size_t double4_vector_size = width * sizeof(double4);
  const size_t size_t_vector_size = width * sizeof(size_t);
  checkCudaErrors(cudaMalloc((void**)&d_x, matrix_size));
  checkCudaErrors(cudaMalloc((void**)&d_y, matrix_size));
  checkCudaErrors(cudaMalloc((void**)&d_col_sizes, size_t_vector_size));
  checkCudaErrors(cudaMalloc((void**)&d_results, double4_vector_size));
  // copying data
  checkCudaErrors(cudaMemcpy(d_x, h_x, matrix_size, cudaMemcpyHostToDevice));
  checkCudaErrors(cudaMemcpy(d_y, h_y, matrix_size, cudaMemcpyHostToDevice));
  checkCudaErrors(cudaMemcpy(d_col_sizes, h_col_sizes, size_t_vector_size, cudaMemcpyHostToDevice));
  checkCudaErrors(cudaMemcpy(d_results, h_results, double4_vector_size, cudaMemcpyHostToDevice));
  // instatiating kernel
  const size_t threadsPerBlock = BLOCK_DIM;
  const size_t blocksPerGrid = calculateBlockNumber(width, threadsPerBlock);
  // calling kernel
  fit_gaussian<<<blocksPerGrid, threadsPerBlock>>>(
    d_x,
    d_y,
    d_col_sizes,
    d_results,
    width,
    max_iterations
  );
  cudaDeviceSynchronize();
  checkCudaErrors(cudaGetLastError());
  //copying memory back
  checkCudaErrors(cudaMemcpy(h_results, d_results, double4_vector_size, cudaMemcpyDeviceToHost));
  // free memory
  cudaFree(d_x);
  cudaFree(d_y);
  cudaFree(d_col_sizes);
  cudaFree(d_results);
}
