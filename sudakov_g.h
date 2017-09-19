// 6.09.2016
// author: Dominik Kasperski
// kasperski.dominik@gmail.com
// 9.10.16
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

//class Sudakov_g{
double theta(const double & delta, const double & z); //wykorzystane funkcje Heaviside'a do ulepszenia
inline double Pgg(const double & z);    // Pgg Kimber 1.41
inline double totPgg(const double & z);    // Pgg Kimber 1.41
inline double totPqq(const double & delta, const double & z);    // Pqq notatki Marcina
inline double Pqg(const double & z);    // Pqg Kimber 1.43
double interpolacja(const double & kt2);
int find_index_gt(const vector<double> & vec, const double& val); // funkcja pomocnicza do wyszukiwania indeksow
//public:
double Tg2(const double & kt2, const double & u2);//liczenie czyste MC
double Tg(const double & kt2, const double & u2); // liczenie prawidłowe
double Tgt(const double & kt2, const double & u2); // liczenie "kombinowane" z trapezami i rozpisane na czesci
double Integrand_as_2(const double & kt2, const double & u2); // calka z as_2 w podanych granicach liczona trapezami
//vector<double> u2s; //dostepne wartosci u2 z pliku
//vector<double> Tgs; //wektor wyliczonych Tg
//vector<double> kt2su; //wektor wykorzystanych kt2
//vector<double> u2su; //wektor wykorzystanych u2
//map<double,double> as_2; //mapa przechowujaca as_2(u^2)
void read_alphas(); //wczytuje plik z as_2(u^2) o nazwie "CTEQ10_alphas_2.dat"
void make_lattice(); // robi siatke jeszcze nie dziala
//int n=100000; //liczba iteracji
double calka_test(const double & lx, const double & ux, const double & ly, const double & uy);
void test_delta(const double &kt2, const double & u2);
void draw_gluons();

///czesc gsl
//pair<double,double>TG(const double & kt2, const double & u2); // funkcja do wywolania, (result,error) blad wyznaczenia to srednia wazona bledow
double ftg(double *args, size_t dim, void *params); //funkcja wymagana przez biblioteke
double ftq(double *args, size_t dim, void *params); //funkcja wymagana przez biblioteke
void warm_up(); //ustawia generator i strukture do calkowania z gsl
void cool_down(); //uwalnia pamiec
struct pars {double u2; }; // co dac jako parametr jak calkuje sie po wszystkim? nullptr! czyli 0, wtedy tylko trzeba go w funkcji zwracajacej double przekazywanej do gsl_monte_function np zrzutowac na void, zeby nie bylo, ze sie nie uzywa zmiennej. Ale to akurat nie ten przypadek, ze wzgledu na delte trzeba przekazac u2
//size_t calls = 10000; //dla beki liczba iteracji
//gsl_monte_function H;
//};


//czesc do rozkladow

double a(const double &, const double &);
double fad(const double &, const double &, const double &);
double fa(const double &, const double &, const double &);
//void draw(const double & x, const double & mu2, const double & kt2ul);

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


double Tq(const double &, const double &);

#endif
