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
double h=0.1;


double Tq(const double & kt2, const double & u2)
{
  if(kt2>=u2)
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
  if(kt2>=u2)
  return 1.;
  vector<double>::iterator test=find(u2s.begin(), u2s.end(), u2);
  // gsl_rng_set(r, chrono::system_clock::now().time_since_epoch().count());
  clock_t t;
  t=clock();
  double result=0., error=0.;		// result and error
  struct pars pms={u2};
  H.params=&pms;
  s = gsl_monte_vegas_alloc(2);
  double xl[2]={kt2, 0.};
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

  t=clock()-t;
  cout<<"time TG: "<<((double)t)/CLOCKS_PER_SEC<<endl;
  //cout<<result<<endl;
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
    //cout<<val<<" "<<vec.size()<<endl;
    bool not_match=true;
    while(not_match)
    {
        //cout<<index<<endl;
        index=(right+left)/2;
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


inline double Pqg(const double & z)
{
    return (z*z+(1.-z)*(1.-z))/2.; // do przerobienia na przesuwanie bitow
}

void read_alphas()
{
    fstream f;
    double tmp1,tmp2;
    f.open("CTEQ10_alphas_2.dat", ios::in | ios::out);
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
//      interpolacja liniowa
    int index = find_index_gt(u2s,kt2);
    double x1=u2s[index];
    double x2=u2s[index-1];
    double y1=as_2[x1];
    double y2=as_2[x2];
    return y1+(y2-y1)/(x2-x1)*(kt2-x1);
}


double ftg(double *args, size_t dim, void *params)
{
    struct pars * fp = (struct pars *)params;
    //prams := {u2}
    //args := {pt2,z}
    return interpolacja(args[0])/(args[0]*2*M_PI)*(totPgg(args[1])*theta(sqrt(args[0])/(sqrt(args[0])+sqrt(fp->u2)),args[1])+5.0*Pqg(args[1]));
}

double ftq(double *args, size_t dim, void *params)
{
    struct pars * fp = (struct pars *)params;
    //prams := {u2}
    //args := {pt2,z}
    return interpolacja(args[0])/(args[0]*2*M_PI)*totPqq(sqrt(args[0])/(sqrt(args[0])+sqrt(fp->u2)),args[1]);
}

void draw_gluons()
{
    vector<double> kt2s={
        1.689999
    };
    vector<double> xs={
    0.99,
0.9735000167,
0.9570000333
    };
    vector<double> mu2s={
    1.689999,
1.819711834,
1.959380543,
2.109769272,
2.271700818,
2.446061129,
2.633804151,
2.835957052,
3.053625835,
3.288001394,
3.540366027,
3.812100454,
4.104691368,
4.419739572,
4.758968734,
5.124234819,
5.517536247,
5.941024818,
6.397017493,
6.888009065,
7.416685813,
7.985940193,
8.598886669,
9.25887875,
9.969527336,
10.73472048,
11.55864464,
12.44580762,
13.4010632,
14.4296377,
15.53715859,

        };
    fstream save;
    string NAZWA="fragment_tg_10";

    save.open(NAZWA,ios::out);

    // for(double & x : xs)
    // {
    for(double & val : kt2s)
        {
        for(double & mu2:mu2s)
            {
            save<<val<<"\t"<<mu2<<"\t"<<Tg(val,mu2)<<endl;
            cout/*<<x*/<<"\t"<<val<<"\t"<<mu2<<endl;
            }
        }
        cout<<NAZWA<<"\t"/*<<x*/<<endl;
    // }

    save.close();
}

double a(const double & x, const double & kt2)
{
    return Ng*pow(x,-lambdag)*pow((1-x),Bg)*log(kt2/lambda0);
}

double fad(const double & x, const double & la2, const double & mu2)
{
    return a(x,la2)*Tg(la2,mu2);
}

double fa(const double & x, const double & kt2, const double & mu2)
{
    return (fad(x,kt2+h,mu2)-fad(x,kt2-h,mu2))/(2*h);
}

double ui(const double & x, const double & kt2)
{
    return 2*Ns*pow(x,-lambdas)*pow((1-x),betas)*log(kt2/lambda0);
}

double di(const double & x, const double & kt2)
{
    return 2*Ns*pow(x,-lambdas)*pow((1-x),betas)*log(kt2/lambda0);
}

double si(const double & x, const double & kt2)
{
    return Ns*pow(x,-lambdas)*pow((1-x),betas)*log(kt2/lambda0);
}

double u(const double & x, const double & kt2)
{
    return ui(x,kt2)+Nu*pow(x,alfau)*pow((1-x),betau)*log(kt2/lambda0);
}

double d(const double & x, const double & kt2)
{
    return di(x,kt2)+Nd*pow(x,alfad)*pow((1-x),betad)*log(kt2/lambda0);
}

double ss(const double & x, const double & kt2)
{
    return si(x,kt2);
}

double fuid(const double & x, const double & la2, const double & mu2)
{
    return ui(x,la2)*Tq(la2,mu2);
}

double fdid(const double & x, const double & la2, const double & mu2)
{
    return di(x,la2)*Tq(la2,mu2);
}

double fsid(const double & x, const double & la2, const double & mu2)
{
    return si(x,la2)*Tq(la2,mu2);
}

double fdd(const double & x, const double & la2, const double & mu2)
{
    return d(x,la2)*Tq(la2,mu2);
}

double fud(const double & x, const double & la2, const double & mu2)
{
    return u(x,la2)*Tq(la2,mu2);
}

double fsd(const double & x, const double & la2, const double & mu2)
{
    return ss(x,la2)*Tq(la2,mu2);
}

//double GluonDF::pochodna(double (GluonDF::*ff)(const double &, const double &, const double &), const double & x, const double & kt2, const double & mu2)
//{
//        return ((this->ff)(x,kt2+h,mu2)-(this_>ff)(x,kt2-h,mu2))/(2*h);
//}

double fu(const double & x, const double & kt2, const double & mu2)
{
    return (fud(x,kt2+h,mu2)-fud(x,kt2-h,mu2))/(2*h);
}


double fd(const double & x, const double & kt2, const double & mu2)
{
    return (fdd(x,kt2+h,mu2)-fdd(x,kt2-h,mu2))/(2*h);
}


double fs(const double & x, const double & kt2, const double & mu2)
{
    return (fsd(x,kt2+h,mu2)-fsd(x,kt2-h,mu2))/(2*h);
}


double fui(const double & x, const double & kt2, const double & mu2)
{
    return (fuid(x,kt2+h,mu2)-fuid(x,kt2-h,mu2))/(2*h);
}


double fdi(const double & x, const double & kt2, const double & mu2)
{
    return (fdid(x,kt2+h,mu2)-fdid(x,kt2-h,mu2))/(2*h);
}


double fsi(const double & x, const double & kt2, const double & mu2)
{
    return (fsid(x,kt2+h,mu2)-fsid(x,kt2-h,mu2))/(2*h);
}