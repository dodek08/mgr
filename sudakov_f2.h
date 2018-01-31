// 5.1.2018
// author: Dominik Kasperski
// kasperski.dominik@gmail.com
/// gsl
/// LHA PDF

#ifndef SUDAKOV_F2
#define SUDAKOV_F2

#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <cstring>
#include <iomanip>
#include <tuple>
#include <vector>
#include <map>
#include <algorithm>
#include <iterator>
#include <random>
#include <chrono>


#include <gsl/gsl_rng.h>
#include <utility>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>

#include "LHAPDF/LHAPDF.h"
// #include "LHAPDFWrap.h"

using namespace std;

double FL(const double &, const double &);
double fl_u(double *args, size_t dim, void *params);
double FT(const double &, const double &);
double ft_u(double *args, size_t dim, void *params);
void warm_up_f2();
void cool_down(); //uwalnia pamiec
struct pars2 {double x; double Q2; double mq2; };


double mu_2_v(const double & kt2, const double & Kt2, const double & mq2); //value of mu2 Kimber 3.15
double D1(const double & Kt2, const double & B, const double & Q2, const double & mq2); //Kimber 3.13
double D2(const double & fi, const double & B, const double & Q2, const double & mq2); //Kimber 3.14
double inv_z(const double & fi, const double & Kt2, const double & B, const double & Q2, const double & mq2); // 1/z, Kimber 3.10
double heaviside(const double & x, const double & fi, const double & Kt2, const double & B, const double & Q2, const double & mq2);

#endif