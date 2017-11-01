// 6.09.2016
// author: Dominik Kasperski
// kasperski.dominik@gmail.com
// 9.10.16
//	wszystkie funkcje dzialaja poprawnie 19.09.2017
///gsl

#ifndef SUDAKOV_G
#define SUDAKOV_G


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

double theta(const double & delta, const double & z); //wykorzystane funkcje Heaviside'a do ulepszenia
inline double Pgg(const double & z);    // Pgg Kimber 1.41
inline double totPgg(const double & z);    // Pgg Kimber 1.41
inline double totPqq(const double & delta, const double & z);    // Pqq notatki Marcina
inline double Pqg(const double & z);    // Pqg Kimber 1.43
double interpolacja(const double & kt2);
int find_index_gt(const vector<double> & vec, const double& val); // funkcja pomocnicza do wyszukiwania indeksow
double Tg(const double & kt2, const double & u2); // liczenie prawidłowe
double Tq(const double &, const double &);
void read_alphas(); //wczytuje plik z as_2(u^2) o nazwie "CTEQ10_alphas_2.dat"
void test_delta(const double &kt2, const double & u2); //testuje dzialanie theta

///czesc gsl
double ftg(double *args, size_t dim, void *params); //funkcja wymagana przez biblioteke
double ftq(double *args, size_t dim, void *params); //funkcja wymagana przez biblioteke
void warm_up(); //ustawia generator i strukture do calkowania z gsl
void cool_down(); //uwalnia pamiec
struct pars {double u2; }; // co dac jako parametr jak calkuje sie po wszystkim? nullptr! czyli 0, wtedy tylko trzeba go w funkcji zwracajacej double przekazywanej do gsl_monte_function np zrzutowac na void, zeby nie bylo, ze sie nie uzywa zmiennej. Ale to akurat nie ten przypadek, ze wzgledu na delte trzeba przekazac u2

#endif
