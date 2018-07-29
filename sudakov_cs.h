// 28.03.2018
// author: Dominik Kasperski
// kasperski.dominik@gmail.com
/// gsl
/// LHA PDF

#ifndef SUDAKOV_CS
#define SUDAKOV_CS

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

double cs_tau(const double &, const double &, const double &);
double cs_tau_f2(const double &, const double &, const double &);
void set_pdf_name_sudakov_cs(string);


class Wyniki
{
public:
	string nazwa;
	vector<double> x;
	vector<double> y;
	vector<double> Q2;
	vector<double> sigma;
	vector<double> totnoproc;
	Wyniki(string);
};

#endif