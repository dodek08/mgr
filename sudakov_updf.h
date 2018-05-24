// 6.09.2016
// author: Dominik Kasperski
// kasperski.dominik@gmail.com
// 1.11.17
// wszystkie funkcje dzialaja poprawnie 19.09.2017
/// gsl
/// LHA PDF

#ifndef SUDAKOV_UPDF
#define SUDAKOV_UPDF


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



//przyspieszenie obliczen dla rozkladu z Kimber 3.22

#define Ng 1.22
#define lambdag 0.0967
#define Bg 1.85
#define lambda0 1.
#define Cg 3. //casimir gluony

#define Cq 1.333333 //casimir kwarki, jednoczesnie Cf

//do rozkladow kwarkow
#define Ns  0.107
#define lambdas 0.0994
#define betas 8.56

#define Nu 9.40
#define alfau 0.986
#define betau 3.85

#define Nd 0.000275
#define alfad 0.000275
#define betad 5.54




using namespace std;
//czesc do rozkladow
void set_pdf_name_sudakov_updf(string);

void set_beta_g(const double &);
double get_beta_g();
void set_n_g(const double &);
double get_n_g();
void set_lambda_g(const double &);
double get_lambda_g();

void draw_gluons(const int & pid); //wylicza kilkadziesiat punktow z rozkladu gluonów
void read_Ts();
vector<double> subset_with_sort(vector<double>&);
ostream& operator<<(ostream& os,const vector<double>&);

double kt2max();
double kt2min();
double mu2min();
double mu2max();

double Tgs(const double & kt2, const double & mu2);
double Tqs(const double & kt2, const double & mu2);


double q_pid(const double&, const double&, const int&);
double fq_pid_d(const double & x, const double & la2, const double & mu2, const int & pid);
double fq_pid(const double & x, const double & kt2, const double & mu2, const int & pid);

double a(const double &, const double &);
double fad(const double &, const double &, const double &);
double fa(const double &, const double &, const double &);

double ui(const double &, const double &);
double di(const double &, const double &);
double si(const double &, const double &);

double u(const double &, const double &);
double d(const double &, const double &);
double ss(const double &, const double &);


double fuid(const double &, const double &, const double &);
double fdid(const double &, const double &, const double &);
double fsid(const double &, const double &, const double &);

double fud(const double &, const double &, const double &);
double fdd(const double &, const double &, const double &);
double fsd(const double &, const double &, const double &);

double fui(const double &, const double &, const double &);
double fdi(const double &, const double &, const double &);
double fsi(const double &, const double &, const double &);

double fu(const double &, const double &, const double &);
double fd(const double &, const double &, const double &);
double fs(const double &, const double &, const double &);

#endif