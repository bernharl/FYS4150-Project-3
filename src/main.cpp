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
    double analytical =  5 * PI * PI / (16 * 16);
    double alpha = 2;
    int N = 1e7;

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


    return 0;
}