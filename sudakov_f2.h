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

double F2(const double &, const double &);
double F2_u(const double &, const double &);
double f2_u(double *args, size_t dim, void *params);
void warm_up_f2();


#endif