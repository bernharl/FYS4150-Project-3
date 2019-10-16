# ifndef MAIN_H
# define MAIN_H

# include <iostream> 
# include <fstream>
# include <cmath> 
# include <random>
# include <ctime>
# include "weights.cpp"

double int_function(double x);
void gauss_laguerre(double *, double *, int, double);
double trapezoidal_rule ( double, double, int, double (*func)(double) );
double simpson ( double, double, int, double (*func)(double) );
void gauleg(double, double, double *, double *, int);
double gammln(double);

#endif 