# include <omp.h>
# include "weights.h"
# include <stdexcept>
# define PI 3.1415926535897932384626433
using std::cout;
using std::cos;
using std::sin;
using std::sqrt;
using std::exp;
using std::string;
#include "integration.cpp"


int main(int argc, char *argv[])
{
  int number_of_threads;
  string filename_O;
  if (argc != 3) {
    throw std::invalid_argument("Please provide exactly two command line arguments. The first deciding the number of threads, use -1 to use all available. The should be the same as the compiler flag used!");
  }
  else {
    number_of_threads = atoi(argv[1]);
    filename_O = argv[2];
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


  /*
  Output for exercise e)
  This part of the code has to be run separately for three different
  compiler flags to produce three separate text files with the data
  needed for exercise e).
  */
  int alpha = 2;
  int N;
  double t_start;
  double t_end;
  double CPU_time;
  string filename = "montecarlo_paro";
  filename += filename_O;
  string extension = ".txt";
  filename += extension;
  ofstream outfilepar;
  outfilepar.open(filename);
  outfilepar << " N: " << " Integral: " << " Variance: " << "CPU_time" << endl;
  for (int i=1; i<=8; i++)
  {
    N = std::pow(10, i);
    t_start = omp_get_wtime();
    std::pair<double, double> results_MC = monte_carlo_improved(N, alpha, 1);
    t_end = omp_get_wtime();
    CPU_time = 1000.0 * (t_end - t_start);
    outfilepar << setw(20) << setprecision(10) << N
               << setw(20) << setprecision(10) << results_MC.first
               << setw(20) << setprecision(10) << results_MC.second
               << setw(20) << setprecision(10) << CPU_time
               << endl;
  }
  outfilepar.close();

}
