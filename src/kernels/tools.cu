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
		const double x = xs[idx], y = ys[idx];

		u = double4Op(cDouble4(x, pow(x, 2), y, x*y), c1, cuda_minus<double>());
		t = double4Op(cDouble4(x_sum, x2_sum, y_sum, xy_sum), u, cuda_plus<double>());
		x_sum 	= t.x;
		x2_sum 	= t.y;
		y_sum	= t.z;
		xy_sum	= t.w;
		
		const double4 tmp = double4Op(t, cDouble4(x_sum, x2_sum, y_sum, xy_sum), cuda_minus<double>());
		c1 = double4Op(tmp, u, cuda_minus<double>());

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