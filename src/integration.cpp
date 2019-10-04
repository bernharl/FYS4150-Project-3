# include "weights.h"
# define PI 3.1415926535897932384626433
using std::cout;

inline double sin_mapping(double y)
{
    return acos(1 - 2 * y);
}

double int_func_cart(double alpha, double x1, double y1, double z1, double x2, double y2, double z2)
{
    double r1 = sqrt(x1 * x1 + y1 * y1 + z1 * z1);
    double r2 = sqrt(x2 * x2 + y2 * y2 + z2 * z2);
    double r12 = sqrt((x1 - x2) * (x1 - x2) +
                      (y1 - y2) * (y1 - y2) +
                      (z1 - z2) * (z1 - z2));
    //cout << r12 << endl;
    //exit(0);
    //cout << exp(-2 * alpha * (r1 + r2)) / r12 << endl;
    //exit(0);
    if (r12 <= EPS)
    {
      //cout << "jeff" << endl;
      return 0;
    }
    //cout << "jeff" << endl;
    //cout << exp(-2 * alpha * (r1 + r2)) / r12 << endl;
    return exp(-2 * alpha * (r1 + r2)) / r12;

}
double int_func_spherical(double alpha, double u1, double u2,
                          double theta1, double theta2,
                          double phi1, double phi2)
{

    double cos_b = cos(theta1) * cos(theta2) + sin(theta1) * sin(theta2) * cos(phi1 - phi2);
    double r12 = sqrt(u1 * u1 + u2 * u2 - 2 * u1 * u2 * cos_b) / (2 * alpha);
    if (r12 <= EPS)
    {
      return 0;
    }
    return sin(theta1) * sin(theta2) / r12;

}

double gauleg_quad(double a, double b, int N, double alpha)
{
    double *x = new double[N];
    double *w = new double[N];

    double I = 0;
    gauleg(a, b, x, w, N);
    for (int i = 0; i < N; i++){
            //cout << x[i] << endl;
    for (int j = 0; j < N; j++){
    for (int k = 0; k < N; k++){
    for (int l = 0; l < N; l++){
    for (int m = 0; m < N; m++){
    for (int n = 0; n < N; n++){
        //cout << I << endl;
        //cout << int_func( alpha, x[i], x[j], x[k], x[l], x[m], x[n]) << endl;
        I += w[i] * w[j] * w[k] * w[l] * w[m] * w[n]
             * int_func_cart( alpha, x[i], x[j], x[k], x[l], x[m], x[n]);
    }}}}}}
    delete [] x;
    delete [] w;
    return I;
}

double gauss_quad_improved(int N, double alpha)
{
    double *u = new double[N];
    double *theta = new double[N];
    double *phi = new double[N];

    double *w_u = new double[N];
    double *w_theta = new double[N];
    double *w_phi = new double[N];

    double I = 0;
    gauss_laguerre(u, w_u, N, 2);
    gauleg(0, PI, theta, w_theta, N);
    gauleg(0, 2 * PI, phi, w_phi, N);
    for (int i = 0; i < N; i++){
    for (int j = 0; j < N; j++){
    for (int k = 0; k < N; k++){
    for (int l = 0; l < N; l++){
    for (int m = 0; m < N; m++){
    for (int n = 0; n < N; n++){
        //cout << I << endl;
        //cout << int_func( alpha, x[i], x[j], x[k], x[l], x[m], x[n]) << endl;
        I += w_r[i] * w_r[j] * w_theta[k] * w_theta[l] * w_phi[m] * w_phi[n]
             * int_func_spherical(alpha, u[i], u[j], theta[k], theta[l], phi[m], phi[n]);
    }}}}}}

    delete [] u;
    delete [] theta;
    delete [] phi;

    delete [] w_u;
    delete [] w_theta;
    delete [] w_phi;
    return I;
}

std::pair<double, double> monte_carlo(double a, double b, double N, double lambda, double alpha)
{
    double I;
    double var;
    double f = 0;
    double f_2 = 0;
    std::mt19937 generator (clock());
    std::uniform_real_distribution<double> uniform(-lambda, lambda);
    double x1;
    double y1;
    double z1;
    double x2;
    double y2;
    double z2;
    for (int i = 0; i < N; i++)
    {
        x1 = uniform(generator);
        y1 = uniform(generator);
        z1 = uniform(generator);
        x2 = uniform(generator);
        y2 = uniform(generator);
        z2 = uniform(generator);
        f += int_func_cart(alpha, x1, y1, z1, x2, y2, z2);
        x1 = uniform(generator);
        y1 = uniform(generator);
        z1 = uniform(generator);
        x2 = uniform(generator);
        y2 = uniform(generator);
        z2 = uniform(generator);
        f_2 += int_func_cart(alpha, x1, y1, z1, x2, y2, z2)
             * int_func_cart(alpha, x1, y1, z1, x2, y2, z2);
    }
    I = f * pow(b - a, 6)  / N;
    var = f_2 / N - f * f / (N * N);
    std::pair<double, double> results = make_pair(I, var);
    //results[0] = 12; // Calculated integral
    //results[1] = 13;  // Confidence interval
    return results;
}

std::pair<double, double> improved_monte_carlo(double N, double alpha)
{
    double I;
    double var;
    double f = 0;
    double f_2 = 0;
    std::mt19937 generator (clock());
    std::exponential_distribution<double> distribution(1);
    std::uniform_real_distribution<double> uniform(0, 1);
    std::uniform_real_distribution<double> uniform_phi(0, 2 * PI);

    double u1;
    double u2;
    double theta1;
    double theta2;
    double phi1;
    double phi2;
    for (int i = 0; i < N; i++)
    {
        u1 = distribution(generator);
        u2 = distribution(generator);
        theta1 = sin_mapping(uniform(generator));
        theta2 = sin_mapping(uniform(generator));
        phi1 = uniform_phi(generator);
        phi2 = uniform_phi(generator);
        f += int_func_spherical(alpha, u1, u2, theta1, theta2, phi1, phi2)
           / (exp(-u1 - u2) * sin(theta1) * sin(theta2));
        u1 = distribution(generator);
        u2 = distribution(generator);
        theta1 = sin_mapping(uniform(generator));
        theta2 = sin_mapping(uniform(generator));
        phi1 = uniform_phi(generator);
        phi2 = uniform_phi(generator);
        f_2 += int_func_spherical(alpha, u1, u2, theta1, theta2, phi1, phi2)
             * int_func_spherical(alpha, u1, u2, theta1, theta2, phi1, phi2);
    }
    I = f * pow(2 * PI, 2)  / (N * pow(2 * alpha, 6));
    var = f_2 / N - f * f / (N * N);
    std::pair<double, double> results = make_pair(I, var);
    //results[0] = 12; // Calculated integral
    //results[1] = 13;  // Confidence interval
    return results;
}

int main()
{
    /*
    double lambda = 1;
    double a = - lambda;
    double b = lambda;
    double alpha = 2;
    int N = 20;
    //double integral = gauleg_quad(a, b, N, alpha);
    double integral = gauss_quad_improved(N, alpha);
    cout << integral << " " << 5 * PI * PI / (16 * 16) << endl;
    */
    
    /*
    double lambda = 1.5;
    double a = - lambda;
    double b = lambda;
    double alpha = 2;
    int N = 1e6;
    //double integral = gauleg_quad(a, b, N, alpha);
    std::pair<double, double> results_MC = monte_carlo(a, b, N, lambda, alpha);
    double integral_MC = results_MC.first;
    double confidence_MC = results_MC.second;
    cout << integral_MC << " " << 5 * PI * PI / (16 * 16) << " " << confidence_MC << endl;
    */
    double alpha = 2;
    int N = 1e6;
    //double integral = gauleg_quad(a, b, N, alpha);
    std::pair<double, double> results_MC = improved_monte_carlo(N, alpha);
    double integral_MC = results_MC.first;
    double confidence_MC = results_MC.second;
    cout << integral_MC << " " << 5 * PI * PI / (16 * 16) << " " << confidence_MC << endl;
    
    return 0;
}