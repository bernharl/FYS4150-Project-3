# include "weights.h"
using namespace std;

double int_func(double alpha, double x1, double y1, double z1, double x2, double y2, double z2)
{
    double r1 = sqrt(x1 * x1 + y1 * y1 + z1 * z1);
    double r2 = sqrt(x2 * x2 + y2 * y2 + z2 * z2);
    double r12 = sqrt((x1 - x2) * (x1 * x2) +
                      (y1 - y2) * (y1 * y2) +
                      (y1 - y2) * (y1 * y2));
    return exp(-2 * alpha * (r1 + r2)) / r12;

}

double gauleg_quad(double a, double b, int N, double alpha)
{
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
        cout << I << endl;
        I += w[i] * w[j] * w[k] * w[l] * w[m] * w[n]
             * int_func( alpha, x[i], x[j], x[k], x[l], x[m], x[n]);
    }}}}}}
    delete [] x;
    delete [] w;
    return I;
}

int main()
{
    double lambda = 10;
    double a = lambda;
    double b = - lambda;
    double alpha = 1;
    int N = 10;

    double integral = gauleg_quad(a, b, N, alpha);
    cout << integral << endl;
    return 0;
}
