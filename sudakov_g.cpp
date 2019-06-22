// 6.09.2016
// author: Dominik Kasperski
// kasperski.dominik@gmail.com

#include "sudakov_g.h"
#include "Blad.h"

unsigned seed = chrono::system_clock::now().time_since_epoch().count();
//lagged Fibonacci pseudo number generator, a bit faster than standard from cstdlib
subtract_with_carry_engine<unsigned,24,10,24> generatorFBN (seed);
constexpr double RMAX=generatorFBN.max();
const gsl_rng_type * T = gsl_rng_taus113; //generator do wyboru z gsl
gsl_rng * r=gsl_rng_alloc(T);
gsl_monte_function H; //Tg
gsl_monte_function I; //Tq
gsl_monte_vegas_state *s;

vector<double> u2s; //dostepne wartosci u2 z pliku
vector<double> kt2su; //wektor wykorzystanych kt2
vector<double> u2su; //wektor wykorzystanych u2
map<double,double> as_2;//mapa przechowujaca as_2(u^2)
size_t calls = 1000000; // liczba iteracji

// const LHAPDF::PDF* pdf = LHAPDF::mkPDF("CT10nlo", 0); //wazne!
LHAPDF::PDF* pdf;
void set_pdf_name_sudakov_g(string siatka)
{
pdf = LHAPDF::mkPDF(siatka, 0);
}

double Tq(const double & kt2, const double & u2)
{
  double kt22;
  if(sqrt(kt2)<sqrt(1.79))
  	kt22=1.79;
  else
  	kt22=kt2;
  if(kt22>=u2)
  	return 1.;
  clock_t t;
  t=clock();
  double result=0., error=0.;       // result and error
  struct pars pms={u2};
  I.params=&pms;
  s = gsl_monte_vegas_alloc(2);
  double xl[2]={kt2, 0.};
  double xu[2]={u2,1.};
  gsl_monte_vegas_integrate (&I, xl, xu, 2, calls/10, r, s, &result, &error);
  //s->stage=1;
  do
  {
  result=0.;
  error=0.;
  gsl_monte_vegas_integrate (&I, xl, xu, 2, calls/5, r, s, &result, &error);
  //cout<<"chisq<< "<<s->chisq<<endl;
  //s->stage=1;
  if(s->chisq==0)
  //gsl_monte_vegas_init(s);
  break;
  }
  while ((fabs (s->chisq - 1.0) > 0.35) ); //w celu uzyskania najwiekszej dokladnosci
  gsl_monte_vegas_free(s);

  t=clock()-t;
  cout<<"time TQ: "<<((double)t)/CLOCKS_PER_SEC<<endl;
  //cout<<result<<endl;
  return exp(-result);
}



void test_delta(const double &kt2, const double & u2)
{

	fstream file;
	file.open("wyniki_1_sposob.dat",ios::out);

	fstream file2;
	file2.open("wyniki_2_sposob.dat",ios::out);


    vector<double> kt2s1;
    vector<double> kt2s2;
    vector<double> zs1;
    vector<double> zs2;
    vector<double> res1;
    vector<double> res2;

    if(kt2>u2)
        throw Blad("kt^2 wieksze od u^2", kt2, u2, kt2, u2);
    vector<double>::iterator test=find(u2s.begin(), u2s.end(), u2);
    if((*test!=u2))
        throw Blad("nie ma alfy dla takiego u^2", kt2, u2, kt2, u2);
    else
    {
        double z=0;
        double delta=0;
        double pt2=0;
        vector<double>::iterator it=find(u2s.begin(), u2s.end(), kt2);
        if((*test!=u2))
            throw Blad("nie ma takiego kt^2", kt2, u2, kt2, u2);
        kt2su.push_back(kt2);
        u2su.push_back(u2);

        for(int i=1; i<1000; i++)
        {
            pt2=(generatorFBN()*1.0/(RMAX*1.0))*(u2-kt2)+kt2;
            kt2s1.push_back(pt2);
            delta=sqrt(pt2)/(sqrt(u2)+sqrt(pt2));
            z=(generatorFBN()*1.0/(RMAX*1.0))*(1.-2.*delta)+delta;
            zs1.push_back(z);
            res1.push_back(theta(delta,z));


            pt2=(generatorFBN()*1.0/(RMAX*1.0))*(u2-kt2)+kt2;
            kt2s2.push_back(pt2);
            delta=sqrt(pt2)/(sqrt(u2)+sqrt(pt2));
            z=(generatorFBN()*1.0/(RMAX*1.0));
            zs2.push_back(z);
            res2.push_back(theta(delta,z));
        }

        for(int i=1; i<1000; i++)
        {
            file<<kt2s1[i]<<"\t"<<zs1[i]<<"\t"<<res1[i]<<endl;
            file2<<kt2s2[i]<<"\t"<<zs2[i]<<"\t"<<res2[i]<<endl;
        }
    }
    file2.close();
    file.close();
    double z=0.5;
    double pt2= 2.4e6;
    cout<<interpolacja(2.5634955924640005e6)<<endl;
    cout<<interpolacja(322.)<<endl;
    cout<<interpolacja(pt2)/(pt2*2*M_PI)<<endl;
    cout<<totPgg(z)*theta(sqrt(pt2)/(sqrt(pt2)+sqrt(u2)),z)<<endl;
    cout<<Pqg(z)<<endl;

}


void cool_down()
{
    gsl_rng_free(r);
}


double Tg(const double & kt2, const double & u2)
{

    //zaimplementowac przyblizenie dla malych mu2
double kt22;
  if(sqrt(kt2)<sqrt(1.79))
  	kt22=1.79;
  else
  	kt22=kt2;
  if(kt22>=u2)
  	return 1.;
  double result=0., error=0.;		// result and error
  struct pars pms={u2};
  H.params=&pms;
  s = gsl_monte_vegas_alloc(2);
  double xl[2]={kt22, 0.};
  double xu[2]={u2,1.};
  gsl_monte_vegas_integrate (&H, xl, xu, 2, calls/10, r, s, &result, &error);
  //s->stage=1;
  do
  {
  result=0.;
  error=0.;
  gsl_monte_vegas_integrate (&H, xl, xu, 2, calls/5, r, s, &result, &error);
  //cout<<"chisq<< "<<s->chisq<<endl;
  //s->stage=1;
  if(s->chisq==0)
  //gsl_monte_vegas_init(s);
  break;
  }
  while ((fabs (s->chisq - 1.0) > 0.35) ); //w celu uzyskania najwiekszej dokladnosci
  gsl_monte_vegas_free(s);

  return exp(-result);
}


void warm_up()
{
  gsl_rng_env_setup ();
  gsl_rng_set(r, chrono::system_clock::now().time_since_epoch().count());
  H.f=&ftg;
  H.dim=2;
  I.f=&ftq;
  I.dim=2;
}

int find_index_gt(const vector<double> & vec, const double & val)
{
    int size = vec.size();
    int right=size+1;
    int left=0;
    int index=(right+left)/2;
    bool not_match=true;
    if (val>vec[right-2])
    {
      index = right;
      not_match = false;
      // cout<< "val>vec"<<endl;
    }
    if (val<vec[0])
    {
      index  = 1;
      not_match = false;
      // cout<<"val<vec"<<endl;
    }
    while(not_match)
    {
        index=(right+left)/2;
        // cout<<index<<"\t"<<vec[index]<<endl;
        if (vec[index]>=val and vec[index-1]<val)
            {
            not_match=false;
            break;
            }
        else
            if(val<vec[index])
                right=index;
            else
                left=index;
        if(index>size or index<1)
            throw Blad("nie dziala find_index_gt",(double)index, val, vec[index], vec[index-1]);
    }
    return index;
}

double theta(const double & delta, const double & z)
{
    if((1.-delta-z)>0)
        if((z-delta)>0)
            return 1.;
        else
        {
            return 0.;

        }
    else
    {
        return 0.;

    }
}

inline double Pgg(const double & z)
{
    double ret = 6.*(z*(1-z)+(1-z)/z+z/(1-z));
    if(ret!=ret)
        throw Blad("nan w Pgg!", z, z, z, z);
    else
        return ret;
}


inline double totPgg(const double & z)
{
    double ret = 6.0*z*(1.-z+1./(1.-z))+6.0*((1.-z)/z);
    if(ret!=ret)
        throw Blad("nan w Pgg!", z, z, z, z);
    else
        return ret;
}

inline double totPqq(const double & delta, const double & z)
{
    if((1.-delta-z)>0)
    {
        double ret = Cq*(1+z*z)/(1-z);
        if(ret!=ret)
            throw Blad("nan w totPqq!", delta, z, delta, z);
    else
        return ret;
    }
    else
        return 0.;
}

double Pqq(const double & z)
{
        double ret = Cq*(1+z*z)/(1-z);
        if(ret!=ret)
            throw Blad("nan w Pqq!", z, z, z, z);
    else
        return ret;
}

inline double Pqg(const double & z)
{
    return (z*z+(1.-z)*(1.-z))/2.; // do przerobienia na przesuwanie bitow
}

void read_alphas()
{
    fstream f;
    double tmp1,tmp2;
    f.open("CTEQ10_alphas.dat", ios::in | ios::out);
    if (!f.is_open()){ throw Blad("zly plik wejscia, nie istnieje lub zle wprowadzony");}
    while(!f.eof())
    {
        f>>tmp1;
        f>>tmp2;
        u2s.push_back(tmp1);
        as_2[tmp1]=tmp2;
    }
    f.close();
    sort(u2s.begin(), u2s.end());
}


double interpolacja(const double & kt2)
{
  double ukt2 = kt2;
    if (kt2<=1.3)
    {
      ukt2=1.3;
    }
//      interpolacja liniowa
    int index = find_index_gt(u2s,ukt2);
    double x1=u2s[index];
    double x2=u2s[index-1];
    double y1=as_2[x1];
    double y2=as_2[x2];
    double ret = y1+(y2-y1)/(x2-x1)*(ukt2-x1);
    return ret;
      // return pdf->alphasQ2(kt2); // odpowiednie as_2

}


double ftg(double *args, size_t dim, void *params)
{
    struct pars * fp = (struct pars *)params;
    //prams := {u2}
    //args := {pt2,z}
    return interpolacja(args[0])/(args[0]*2*M_PI)*(args[1]*totPgg(args[1])*theta(sqrt(args[0])/(sqrt(args[0])+sqrt(fp->u2)),args[1])+5.0*Pqg(args[1]));
}

double ftq(double *args, size_t dim, void *params)
{
    struct pars * fp = (struct pars *)params;
    //prams := {u2}
    //args := {pt2,z}
    return interpolacja(args[0])/(args[0]*2*M_PI)*totPqq(sqrt(args[0])/(sqrt(args[0])+sqrt(fp->u2)),args[1]);
}