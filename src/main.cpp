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

    if (argc != 2){
      throw std::invalid_argument("Please provide exactly one command line argument deciding the number of threads, use -1 to use all available");
    }
    else {
      number_of_threads = atoi(argv[1]);
    }

    if (number_of_threads == -1){
      number_of_threads = omp_get_max_threads();
    }
    else if (number_of_threads > omp_get_max_threads()){
      throw std::invalid_argument("Number of threads chosen is higher than the max available!");
    }
    else if (number_of_threads <= 0){
      throw std::invalid_argument("Please choose a positive, non-zero amount of threads. To use all threads, choose -1");
    }
    double analy =  5 * PI * PI / (16 * 16);
    double alpha = 2;
    int N = 1e7;
    /*
    double t_start = omp_get_wtime();
    std::pair<double, double> results_MC_improved = monte_carlo_improved(N, alpha, number_of_threads);
    double t_end = omp_get_wtime();
    double CPU_time = 1000.0 * (t_end - t_start);

    double integral_MC_improved = results_MC_improved.first;
    double confidence_MC_improved = results_MC_improved.second;
    cout << "Improved MC: Time: "<< CPU_time << " ms | Numerical integral: " << integral_MC_improved << "| Analytical: " << analytical << "| Variance: " << confidence_MC_improved << endl;


    double lambda = 1.5;
    double a = - lambda;
    double b = lambda;

    t_start = omp_get_wtime();
    std::pair<double, double> results_MC = monte_carlo(a, b, N, lambda, alpha, number_of_threads);
    t_end = omp_get_wtime();
    CPU_time = 1000.0 * (t_end - t_start);

    double integral_MC = results_MC.first;
    double confidence_MC = results_MC.second;
    cout << "Brute MC: Time: "<< CPU_time << " ms | Numerical integral: " << integral_MC << "| Analytical: " << analytical << "| Variance: " << confidence_MC << endl;
    */

    // Exercise a) output
    /*
    arma::mat gauleg_res = arma::zeros <arma::mat> (21, 21);
    int iterator1 = 1;
    int iterator2;
    
    for (int i=1; i<=20; i++){  
        cout << iterator1 << endl;
        gauleg_res(iterator1, 0) = i;
        iterator1++;
        iterator2 = 1;
    for (double j=1; j<=3; j+=0.1){
        N = std::round(std::pow(10,i));
        gauleg_res(0, iterator2) = N;
        gauleg_res(i, j) = gauleg_quad(-j, j, N, 2.0);
        iterator2++;
    }}
    
    cout << gauleg_res << endl;
    gauleg_res.save("Exercise_a.txt", arma::arma_ascii);
    */
  double lambda = 1.5;
  double error_gauleg;
  double error_gauss_improved;
  /*
  cout << gauss_quad_improved(30, 2.0) << "  " << analy << endl;
  
  ofstream outfile;
  outfile.open("Exercise_a_b.txt");
  outfile << " N: " << " Error gualeg_quad: " << " Error gauss_quad_improved: " << endl;
  for (int i = 1; i <= 6; i++)
  { 
    error_gauleg = std::fabs(gauleg_quad( -lambda, lambda, i, 2.0) - analy);
    error_gauss_improved = std::fabs(gauss_quad_improved( i, 2.0) - analy);
    outfile << setw(20) << setprecision(10) << i
            << setw(20) << setprecision(10) << error_gauleg
            << setw(20) << setprecision(10) << error_gauss_improved
            << endl;
            
  }
  outfile.close();
  */

  lambda = 1.5;
  double a = - lambda;
  double b = lambda;

  double t_start;
  std::pair<double, double> results_MC;
  double t_end;
  double CPU_time;

  double integral_MC = results_MC.first;
  double confidence_MC = results_MC.second;
  /*
  ofstream outfile;
  outfile.open("montecarlo.txt");
  outfile << " N: " << " Integral: " << " Variance: " << "CPU_time" << endl;
  for (int i=1; i<=8; i++)
  {
    N = std::pow(10, i);
    t_start = omp_get_wtime();
    std::pair<double, double> results_MC = monte_carlo(a, b, N, lambda, alpha, 1);
    t_end = omp_get_wtime();
    CPU_time = 1000.0 * (t_end - t_start);

    outfile << setw(20) << setprecision(10) << N 
            << setw(20) << setprecision(10) << results_MC.first
            << setw(20) << setprecision(10) << results_MC.second 
            << setw(20) << setprecision(10) << CPU_time
            << endl;
  }
  outfile.close();

  ofstream outfileimp;
  outfileimp.open("montecarlo_improved.txt");
  outfileimp << " N: " << " Integral: " << " Variance: " << "CPU_time" << endl;
  for (int i=1; i<=8; i++)
  {
    N = std::pow(10, i);
    t_start = omp_get_wtime();
    std::pair<double, double> results_MC = monte_carlo_improved(N, alpha, 1);
    t_end = omp_get_wtime();
    CPU_time = 1000.0 * (t_end - t_start);

    outfileimp << setw(20) << setprecision(10) << N 
               << setw(20) << setprecision(10) << results_MC.first
               << setw(20) << setprecision(10) << results_MC.second 
               << setw(20) << setprecision(10) << CPU_time
               << endl;
  }
  outfileimp.close();
  */
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
  
  /*
  // Exercise b) output
  ofstream outfile;
  outfile.open("Exercise_b.txt");
  outfile << " N " << " Integral " << endl;
  for (int i = 1; i <= 5; i++)
  { 
    N = std::round(std::pow(10, i));
    outfile << setw(20) << setprecision(10) << N 
            << setw(20) << setprecision(10) << gauss_quad_improved( N, 2) 
            << endl;
  }
  outfile.close();
  */
  return 0;
}