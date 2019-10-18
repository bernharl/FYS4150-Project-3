# include <omp.h>
# include "weights.h"
# include <stdexcept>
# define PI 3.1415926535897932384626433
using std::cout;
using std::cos;
using std::sin;
using std::sqrt;
using std::exp;

double int_func_cart(double alpha, double x1, double y1, double z1, double x2, double y2, double z2)
  /*
  Computes the integrand in the integral expression for the quantum mechanical
  expectation value of the correlation energy between two electrons which repel
  each other via the classical Coulomb interaction in cartesian coordinates.
  ------------
  alpha: double
    Constant in the exponential term of the integrand.
  x1, y1, z1: double
  Cartesian coordinates for the vector r1
  x2, y2, z2: double
    Cartesian coordinates for the vector r2
  */
{
    double r1 = sqrt(x1 * x1 + y1 * y1 + z1 * z1);
    double r2 = sqrt(x2 * x2 + y2 * y2 + z2 * z2);
    double r12 = sqrt((x1 - x2) * (x1 - x2) +
                      (y1 - y2) * (y1 - y2) +
                      (z1 - z2) * (z1 - z2));
    if (r12 <= ZERO)
    {
      return 0;
    }
    return exp(-2 * alpha * (r1 + r2)) / r12;

}
double int_func_spherical(double u1, double u2,
                          double theta1, double theta2,
                          double phi1, double phi2)
  /*
  Computes the integrand in the integral expression for the quantum mechanical
  expectation value of the correlation energy between two electrons which repel
  each other via the classical Coulomb interaction in spherical coordinates.
  ------------
  u1, u2: double
    Dimensionless radial distance
  x1, y1, z1: double
  Cartesian coordinates for the vector r1
  x2, y2, z2: double
    Cartesian coordinates for the vector r2
  */
{
    double cos_b = cos(theta1) * cos(theta2) + sin(theta1) * sin(theta2) * cos(phi1 - phi2);
    double r12 = u1 * u1 + u2 * u2 - 2 * u1 * u2 * cos_b;
    if (r12 <= ZERO)
    {
      return 0;
    }
    return 1.0 / sqrt(r12);

}

double debug_integrand(double u1, double u2,
                       double theta1, double theta2,
                       double phi1, double phi2)
  /*
  Integrand function with known analytical integral. Used in unit test.
  ------------
  u1: double
    VENTE MED DETTE!!!!!
  x1, y1, z1: double
  Cartesian coordinates for the vector r1
  x2, y2, z2: double
    Cartesian coordinates for the vector r2
  */
{
    double integrand = 1;
    return integrand;

}

double gauleg_quad(double a, double b, int N, double alpha)
{
  /*
  Calculates the integral in cartesian coordinates
  using gaussian quadrature with legendre plynomials.
  ------------
  a: double
    Lower limit of the integral
  b: double
    Upper limit of the integral
  N: int
    number of grid points between a and b
  alpha: double
    Constant in the exponential term of the integrand.
  */
    double *x = new double[N];
    double *w = new double[N];

    double I = 0;
    gauleg(a, b, x, w, N);
    for (int i = 0; i < N; i++){
    for (int j = 0; j < N; j++){
    for (int k = 0; k < N; k++){
    for (int l = 0; l < N; l++){
    for (int m = 0; m < N; m++){
    for (int n = 0; n < N; n++){
        I += w[i] * w[j] * w[k] * w[l] * w[m] * w[n]
             * int_func_cart( alpha, x[i], x[j], x[k], x[l], x[m], x[n]);
    }}}}}}
    delete [] x;
    delete [] w;
    return I;
}

double gauss_quad_improved(int N, double alpha,
                           double (*int_func)(double, double,
                                              double, double,
                                              double, double))
  /*
  Calculates the integral in spherical coordinates using
  gaussian quadrature with legendre plynomials and
  laguerre polynomials on the radial axis.
  ------------
  N: int
    number of grid points between limits a and b
  alpha: double
    Constant in the exponential term of the integrand.
  int_func: double
    Integrand function.
  */
{   double *u = new double[N + 1];
    double *theta = new double[N];
    double *phi = new double[N];

    double *w_u = new double[N + 1];
    double *w_theta = new double[N];
    double *w_phi = new double[N];

    double I = 0;
    gauss_laguerre(u, w_u, N, 2);
    gauleg(0, PI, theta, w_theta, N);
    gauleg(0, 2 * PI, phi, w_phi, N);
    for (int i = 1; i <= N; i++){
    for (int j = 1; j <= N; j++){
    for (int k = 0; k < N; k++){
    for (int l = 0; l < N; l++){
    for (int m = 0; m < N; m++){
    for (int n = 0; n < N; n++){
        I += w_u[i] * w_u[j] * w_theta[k] * w_theta[l] * w_phi[m] * w_phi[n]
             * int_func(u[i], u[j], theta[k], theta[l], phi[m], phi[n]) * sin(theta[k]) * sin(theta[l]);
    }}}}}}

    delete [] u;
    delete [] theta;
    delete [] phi;

    delete [] w_u;
    delete [] w_theta;
    delete [] w_phi;
    return I / (32 * pow(alpha, 5));
}

std::pair<double, double> monte_carlo(double a, double b, int N, double lambda, double alpha, int number_of_threads)
  /*
  Calculates the integral in cartesian coordinates
  using monte-carlo integration.
  ------------
  a: double
    Lower limit of the integral
  b: double
    Upper limit of the integral
  N: int
    number of grid points between a and b
  alpha: double
    Constant in the exponential term of the integrand.
  number_of_threads: int
    Number of threads on  the
    processor the program is run
  */
{
    double I;
    double var;
    double func_val;
    double f = 0;
    double f_2 = 0;
    uniform_real_distribution<double> uniform(-lambda, lambda);
    mt19937 generator;
    double x1;
    double y1;
    double z1;
    double x2;
    double y2;
    double z2;
    #pragma omp parallel reduction (+:f, f_2) num_threads(number_of_threads) private(x1, x2, y1, y2, z1, z2, func_val, generator)
    generator.seed(omp_get_wtime() + omp_get_thread_num());
    #pragma omp for
    for (int i = 0; i < N; i++)
    {
        x1 = uniform(generator);
        x2 = uniform(generator);
        y1 = uniform(generator);
        y2 = uniform(generator);
        z1 = uniform(generator);
        z2 = uniform(generator);
        func_val = int_func_cart(alpha, x1, y1, z1, x2, y2, z2);
        f += func_val;
        f_2 += func_val * func_val;
    }
    double common_factor = pow(b - a, 6);
    I = f * common_factor  / ((double) N);
    f_2 *= pow(common_factor, 2) / ( (double) N);
    var = f_2 - I * I;
    std::pair<double, double> results = make_pair(I, var);
    return results;
}

std::pair<double, double> monte_carlo_improved(int N, double alpha, int number_of_threads)
  /*
  Calculates the integral in spherical coordinates using monte-carlo
  integration with importance and exponential distribution for the
  radial axis.
  ------------
  N: int
    number of grid points between a and b
  alpha: double
    Constant in the exponential term of the integrand.
  number_of_threads: int
    Number of threads on  the
    processor the program is run
  */
{
  double I;
  double func_val;
  double var;
  double f = 0;
  double f_2 = 0;

  exponential_distribution<double> exponential(1);
  uniform_real_distribution<double> uniform_theta(0, PI);
  uniform_real_distribution<double> uniform_phi(0, 2 * PI);

  mt19937 generator;

  double u1;
  double u2;
  double theta1;
  double theta2;
  double phi1;
  double phi2;
  #pragma omp parallel reduction (+:f, f_2) num_threads(number_of_threads) private(u1, u2, theta1, theta2, phi1, phi2, func_val, generator)
  generator.seed(omp_get_wtime() + omp_get_thread_num());
  #pragma omp for
  for (int i = 0; i < N; i++)
  {
      u1 = exponential(generator);
      u2 = exponential(generator);
      theta1 = uniform_theta(generator);
      theta2 = uniform_theta(generator);
      phi1 = uniform_phi(generator);
      phi2 = uniform_phi(generator);
      func_val = int_func_spherical(u1, u2, theta1, theta2, phi1, phi2)  * u1 * u1 * u2 * u2 * sin(theta1) * sin(theta2);
      f += func_val;
      f_2 += func_val * func_val;
  }

  double common_factor = 4 * pow(PI, 4) / pow(2 * alpha, 5);
  I = f * common_factor / ((double) N);

  f_2 *= common_factor * common_factor / ((double) N);
  var = f_2 - I * I;

  std::pair<double, double> results = make_pair(I, var);
  return results;
}
