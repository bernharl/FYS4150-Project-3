# include <omp.h>
# include "weights.h"
# include <stdexcept>
# define PI 3.1415926535897932384626433
using std::cout;
using std::cos;
using std::sin;
using std::sqrt;
using std::exp;

#include "integration.cpp"


int main(int argc, char *argv[])
{
  int number_of_threads;

  if (argc != 2) {
    throw std::invalid_argument("Please provide exactly one command line argument deciding the number of threads, use -1 to use all available");
  }
  else {
    number_of_threads = atoi(argv[1]);
  }

  if (number_of_threads == -1) {
      number_of_threads = omp_get_max_threads();
  }

  else if (number_of_threads > omp_get_max_threads()) {
      throw std::invalid_argument("Number of threads chosen is higher than the max available!");
  }

  else if (number_of_threads <= 0) {
    throw std::invalid_argument("Please choose a positive, non-zero amount of threads. To use all threads, choose -1");
  }

  double analy =  5 * PI * PI / (16 * 16); // Analytical solution to the integral.
  double alpha = 2; // Constant factor in the exponential function of the integrand.
  int N = 1e7; // Grid resolution for gaussian quadrature 
               // and sample size for Monte Carlo methods.




  double lambda = 1.5; //Factor determining upper and lower limits.
  double error_gauleg; // Error for Gaussian quadrature using laguerre polynomials.
  double error_gauss_improved; // Error for Gaussian quadrature using legendre plynomials.

  // Exercise a) and b) output
  ofstream outfileab;
  outfileab.open("Exercise_a_b.txt");
  outfileab << " N: " << " Error gualeg_quad: " << " Error gauss_quad_improved: " << endl;
  for (int i = 1; i <= 6; i++)
  { 
    error_gauleg = std::fabs(gauleg_quad( -lambda, lambda, i, 2.0) - analy);
    error_gauss_improved = std::fabs(gauss_quad_improved(i, 2.0, &int_func_spherical) - analy);
    outfileab << setw(20) << setprecision(10) << i
            << setw(20) << setprecision(10) << error_gauleg
            << setw(20) << setprecision(10) << error_gauss_improved
            << endl;
  cout << i << endl;
  }
  outfileab.close();

 
  lambda = 1.5;
  double a = - lambda; // Upper limit of the integral.
  double b = lambda;   // Lower limit of the integral.

  double t_start;
  std::pair<double, double> results_MC;
  std::pair<double, double> results_MC_par;
  double t_end;
  double CPU_time;
  double CPU_time_par;
  double integral_MC = results_MC.first;
  double confidence_MC = results_MC.second;
  

  //Output for exercise c) and d)
  ofstream outfilecd;
  outfilecd.open("montecarlo.txt");
  outfilecd << " N: " << " Integral: " << " Variance: " << " CPU_time " << " CPU_time_par " << endl;
  for (int i=1; i<=8; i++)
  {
    N = std::pow(10, i);
    t_start = omp_get_wtime();
    std::pair<double, double> results_MC_par = monte_carlo(a, b, N, lambda, alpha, 2);
    t_end = omp_get_wtime();
    CPU_time_par = 1000.0 * (t_end - t_start);

    t_start = omp_get_wtime();
    std::pair<double, double> results_MC = monte_carlo(a, b, N, lambda, alpha, 1);
    t_end = omp_get_wtime();
    CPU_time = 1000.0 * (t_end - t_start);

    outfilecd << setw(20) << setprecision(10) << N 
            << setw(20) << setprecision(10) << results_MC.first
            << setw(20) << setprecision(10) << results_MC.second 
            << setw(20) << setprecision(10) << CPU_time
            << setw(20) << setprecision(10) << CPU_time_par
            << endl;
  }
  outfilecd.close();

  ofstream outfileimp;
  outfileimp.open("montecarlo_improved.txt");
  outfileimp << " N: " << " Integral: " << " Variance: " << "CPU_time" << endl;


  // Output for exercise c) and d)
  for (int i=1; i<=8; i++)
  {
    N = std::pow(10, i);

    t_start = omp_get_wtime();
    std::pair<double, double> results_MC_par = monte_carlo_improved(N, alpha, 2);
    t_end = omp_get_wtime();
    CPU_time_par = 1000.0 * (t_end - t_start);

    t_start = omp_get_wtime();
    std::pair<double, double> results_MC = monte_carlo_improved(N, alpha, 1);
    t_end = omp_get_wtime();
    CPU_time = 1000.0 * (t_end - t_start);

    outfileimp << setw(20) << setprecision(10) << N 
               << setw(20) << setprecision(10) << results_MC.first
               << setw(20) << setprecision(10) << results_MC.second 
               << setw(20) << setprecision(10) << CPU_time
               << setw(20) << setprecision(10) << CPU_time_par
               << endl;
  }
  outfileimp.close();
  
  
  /*
  //Output for exercise e)
  //This part of the code has to be run separately for three different
  //compiler flags to produce three separate text files with the data
  //needed for exercise e).

  ofstream outfilepar;
  outfilepar.open("montecarlo_paro3.txt");
  outfilepar << " N: " << " Integral: " << " Variance: " << "CPU_time" << endl;
  for (int i=1; i<=8; i++)
  {
    N = std::pow(10, i);
    t_start = omp_get_wtime();
    cout << "hei" << endl;
    std::pair<double, double> results_MC = monte_carlo_improved(N, alpha, number_of_threads);
    t_end = omp_get_wtime();
    CPU_time = 1000.0 * (t_end - t_start);
    outfilepar << setw(20) << setprecision(10) << N 
               << setw(20) << setprecision(10) << results_MC.first
               << setw(20) << setprecision(10) << results_MC.second 
               << setw(20) << setprecision(10) << CPU_time
               << endl;
  }
  outfilepar.close();
  */
  return 0;
}