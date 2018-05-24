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

double FL_g(const double &, const double &, size_t);
double FL_g(const double &, const double &);
double fl_u_g(double *args, size_t dim, void *params);
double FT_g(const double &, const double &, size_t);
double FT_g(const double &, const double &);
double F2_q(const double &, const double &);
double ft_u_g(double *args, size_t dim, void *params);
double f2_u_q(double *args, size_t dim, void *params);
void warm_up_f2();
void cool_down(); //uwalnia pamiec
struct pars2 {double x; double Q2; double mq2; };
struct parsf2 {double x; double Q2; int pid; };


double mu_2_v(const double & kt2, const double & Kt2, const double & mq2); //value of mu2 Kimber 3.15
double D1(const double & Kt2, const double & B, const double & Q2, const double & mq2); //Kimber 3.13
double D2(const double & Kt2, const double & kt2, const double & fi, const double & B, const double & Q2, const double & mq2); //Kimber 3.14
double z(const double & Kt2, const double & kt2, const double & fi, const double & B, const double & Q2, const double & mq2); // 1/z, Kimber 3.10
double heaviside(const double & val); //theta(1-x/z)!

double fgk(const double &, const double&);
void set_pdf_name_sudakov_f2(string);



#endif