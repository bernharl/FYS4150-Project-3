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




int main()
{

  double analy =  5 * PI * PI / (16 * 16); // Analytical solution to the integral.
  double alpha = 2; // Constant factor in the exponential function of the integrand.
  int N = 1e7; // Grid resolution for gaussian quadrature
               // and sample size for Monte Carlo methods.

  double lambda; //Factor determining upper and lower limits.
  double error_gauleg; // Error for Gaussian quadrature using laguerre polynomials.
  double error_gauss_improved; // Error for Gaussian quadrature using legendre plynomials.

  lambda = 1.5;
  double a = - lambda; // Upper limit of the integral.
  double b = lambda;   // Lower limit of the integral.

  double t_start;
  std::pair<double, double> results_MC;
  std::pair<double, double> results_MC_par;
  double t_end;
  double CPU_time;
  double CPU_time_par;

  // Exercise a) and b) output
  cout << "Initializing gaussian quadrature brute force and improved" << endl;

  ofstream outfileab;
  outfileab.open("Exercise_a_b.txt");
  outfileab << " N: " << " Error gualeg_quad: " << " Error gauss_quad_improved: " << endl;
  for (int i = 1; i <= 30; i++)
  {

    t_start = omp_get_wtime();
    error_gauleg = std::fabs(gauleg_quad( -lambda, lambda, i, 2.0) - analy);
    t_end = omp_get_wtime();
    CPU_time = 1000.0 * (t_end - t_start);
    if (i == 30)
    {
      cout << "CPU time gaussian quadrature  for N = 30:  " << CPU_time << " ms" <<endl;
    }
    t_start = omp_get_wtime();
    error_gauss_improved = std::fabs(gauss_quad_improved(i, 2.0, &int_func_spherical) - analy);
    t_end = omp_get_wtime();
    CPU_time = 1000.0 * (t_end - t_start);
        if (i == 30)
    {
      cout << "CPU time gaussian quadrature improved for N = 30: " << CPU_time << " ms" << endl;
    }

    outfileab << setw(20) << setprecision(10) << i
              << setw(20) << setprecision(10) << error_gauleg
              << setw(20) << setprecision(10) << error_gauss_improved
              << endl;
  }
  outfileab.close();



  
  //Output for exercise c) and d)

  lambda = 1.5;
  cout << "Computing regular Monte Carlo method with lambda = 2" << endl;
  ofstream outfilecd_2;
  outfilecd_2.open("montecarlo_lambda_2.txt");
  outfilecd_2 << " N: " << " Integral: " << " Variance: " << " CPU_time " << " CPU_time_par " << endl;
  for (int i=1; i<=8; i++)
  {
    N = std::pow(10, i);
    t_start = omp_get_wtime();
    results_MC_par = monte_carlo(-lambda, lambda, N / 2.0 , alpha, 2);
    t_end = omp_get_wtime();
    CPU_time_par = 1000.0 * (t_end - t_start);

    if (i == 8)
    {
      cout << "CPU time Monte Carlo two threads with N = 1e8: " << CPU_time_par << " ms" << endl;
    }

    t_start = omp_get_wtime();
    results_MC = monte_carlo(-lambda, lambda, N, alpha, 1);
    t_end = omp_get_wtime();
    CPU_time = 1000.0 * (t_end - t_start);

    if (i == 8)
    {
      cout << "CPU time Monte Carlo one thread with N = 1e8: " << CPU_time << " ms" << endl;
    }


    outfilecd_2 << setw(20) << setprecision(10) << N
            << setw(20) << setprecision(10) << results_MC.first
            << setw(20) << setprecision(10) << results_MC.second
            << setw(20) << setprecision(10) << CPU_time
            << setw(20) << setprecision(10) << CPU_time_par
            << endl;
  }
  outfilecd_2.close();


  lambda = 3;
  cout << "Computing regular Monte Carlo method lambda = 3" << endl;
  ofstream outfilecd_3;
  outfilecd_3.open("montecarlo_lambda_3.txt");
  outfilecd_3 << " N: " << " Integral: " << " Variance: " << " CPU_time " << " CPU_time_par " << endl;
  for (int i=1; i<=8; i++)
  {
    N = std::pow(10, i);
    t_start = omp_get_wtime();
    results_MC_par = monte_carlo(-lambda, lambda, N / 2.0 , alpha, 2);
    t_end = omp_get_wtime();
    CPU_time_par = 1000.0 * (t_end - t_start);

    if (i == 8)
    {
      cout << "CPU time Monte Carlo two threads with N = 1e8: " << CPU_time_par << " ms" << endl;
    }

    t_start = omp_get_wtime();
    results_MC = monte_carlo(-lambda, lambda, N, alpha, 1);
    t_end = omp_get_wtime();
    CPU_time = 1000.0 * (t_end - t_start);

    if (i == 8)
    {
      cout << "CPU time Monte Carlo one thread with N = 1e8: " << CPU_time << " ms" << endl;
    }


    outfilecd_3 << setw(20) << setprecision(10) << N
            << setw(20) << setprecision(10) << results_MC.first
            << setw(20) << setprecision(10) << results_MC.second
            << setw(20) << setprecision(10) << CPU_time
            << setw(20) << setprecision(10) << CPU_time_par
            << endl;
  }
  outfilecd_3.close();


  ofstream outfileimp;
  outfileimp.open("montecarlo_improved.txt");
  outfileimp << " N: " << " Integral: " << " Variance: " << "CPU_time" << endl;


  // Output for exercise c) and d)
  cout << "Computing improved Monte Carlo method" << endl;
  for (int i=1; i<=8; i++)
  {
    N = std::pow(10, i);

    t_start = omp_get_wtime();
    results_MC_par = monte_carlo_improved(N / 2.0, alpha, 2);
    t_end = omp_get_wtime();
    CPU_time_par = 1000.0 * (t_end - t_start);

    if (i == 8)
    {
      cout << "CPU time Monte Carlo improved two threads with N = 1e8: " << CPU_time_par << " ms" << endl;
    }

    t_start = omp_get_wtime();
    results_MC = monte_carlo_improved(N, alpha, 1);
    t_end = omp_get_wtime();
    CPU_time = 1000.0 * (t_end - t_start);
    if (i == 8)
    {
      cout << "CPU time Monte Carlo improved one thread with N = 1e8: " << CPU_time << " ms" << endl;
    }
    outfileimp << setw(20) << setprecision(10) << N
               << setw(20) << setprecision(10) << results_MC.first
               << setw(20) << setprecision(10) << results_MC.second
               << setw(20) << setprecision(10) << CPU_time
               << setw(20) << setprecision(10) << CPU_time_par
               << endl;
  }
  outfileimp.close();

  return 0;
}
