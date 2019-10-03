# include "weights.h"
# define PI 3.1415926535897932384626433
using namespace std;

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
double int_func_spherical(double alpha, double r1, double r2,
                          double theta1, double theta2,
                          double phi1, double phi2)
{
    double cos_b = cos(theta1) * cos(theta2) + sin(theta1) * sin(theta2) * cos(phi1 - phi2);
    double r12 = sqrt(r1 * r1 + r2 * r2 - 2 * r1 * r2 * cos_b);
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
    return r1 * r1 * r2 * r2 * exp(-2 * alpha * (r1 + r2)) / r12;

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
    double *r1 = new double[N+1];
    double *r2 = new double[N+1];
    double *theta1 = new double[N];
    double *theta2 = new double[N];
    double *phi1 = new double[N];
    double *phi2 = new double[N];

    double *w_r1 = new double[N+1];
    double *w_r2 = new double[N+1];
    double *w_theta1 = new double[N];
    double *w_theta2 = new double[N];
    double *w_phi1 = new double[N];
    double *w_phi2 = new double[N];


    double I = 0;
    gauss_laguerre(r1, w_r1, N+1, alpha);
    gauss_laguerre(r2, w_r2, N+1, alpha);
    gauleg(0, PI, theta1, w_theta1, N);
    gauleg(0, PI, theta2, w_theta2, N);
    gauleg(0, 2 * PI, phi1, w_phi1, N);
    gauleg(0, 2 * PI, phi2, w_phi2, N);
    for (int i = 0; i < N; i++){
            //cout << x[i] << endl;
    for (int j = 0; j < N; j++){
    for (int k = 0; k < N; k++){
    for (int l = 0; l < N; l++){
    for (int m = 0; m < N; m++){
    for (int n = 0; n < N; n++){
        //cout << I << endl;
        //cout << int_func( alpha, x[i], x[j], x[k], x[l], x[m], x[n]) << endl;
        I += w_r1[i] * w_r2[j] * w_theta1[k] * w_theta2[l] * w_phi1[m] * w_phi2[n]
             * int_func_spherical(alpha, r1[i], r2[j], theta1[k], theta2[l], phi1[m], phi2[n]);
    }}}}}}
    delete [] r1;
    delete [] r2;
    delete [] theta1;
    delete [] theta2;
    delete [] phi1;
    delete [] phi2;

    delete [] w_r1;
    delete [] w_r2;
    delete [] w_theta1;
    delete [] w_theta2;
    delete [] w_phi1;
    delete [] w_phi2;
    return I;
}

int main()
{
    double lambda = 1;
    double a = - lambda;
    double b = lambda;
    double alpha = 1.9;
    int N = 20;

    //double integral = gauleg_quad(a, b, N, alpha);
    double integral = gauss_quad_improved(N, alpha);
    cout << integral << " " << 5 * PI * PI / (16 * 16) << endl;
    return 0;
}
