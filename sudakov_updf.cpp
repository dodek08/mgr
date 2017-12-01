
#include "sudakov_g.h"
#include "sudakov_updf.h"
#include "Blad.h"


double h=0.1; //przy rozniczkowaniu numerycznym
 const LHAPDF::PDF* pdfs = LHAPDF::mkPDF("CT10nlo", 0); //wazne!
 double timetg=997.;
 double timetq=997.; 

void draw_gluons()
{
    fstream save;
    string NAZWA="siatkaTG_2_CZ_2";

    save.open(NAZWA,ios::out);
    clock_t t,t1;
    t=clock();
    t1 = clock();
     double MINX = 1.0e-6;
double MAXX = 0.99;
double MINKT2 = pow( 0.035, 2.0) ;
double MAXKT2 = pow( 1.0e4, 2.0 );
double MINMU2 = 1.79;
double MAXMU2 = pow( 1.0e4, 2.0 );
double MINLOGX, MINLOGKT2, MINLOGMU2;
double MAXLOGX, MAXLOGKT2, MAXLOGMU2;
int NX, NKT2, NMU2;
double DX, DKT2, DMU2;

MINLOGX = log( MINX );
MAXLOGX = log( MAXX );
MINLOGKT2 = log( MINKT2 );
MAXLOGKT2 = log( MAXKT2 );
MINLOGMU2 = log( MINMU2 );
MAXLOGMU2 = log( MAXMU2 );
NX = 60;
NKT2 = 140;
NMU2 = 120;
DX = ( MAXLOGX - MINLOGX ) / NX;
DKT2 = ( MAXLOGKT2 - MINLOGKT2 ) / NKT2;
DMU2 = ( MAXLOGMU2 - MINLOGMU2 ) / NMU2;
double logx;
double logkt2;
double logmu2;
double x;
double kt2;
double mu2;


//////////////////////////////////////////
/////// ATTENTION PLEASE!/////////////////
 MINLOGKT2=(MINLOGKT2+MINLOGKT2+DKT2)/2;//
//////////////////////////////////////////
// for( int ix=0; ix<NX+1; ix++ ){
//   logx = MINLOGX + ix*DX;
//   x = exp(logx);
  for( int ikt2=0; ikt2<NKT2+1; ikt2++ ){
    logkt2 = MINLOGKT2 + ikt2*DKT2;
        kt2 = exp(logkt2);
        t=clock();
    for( int imu2=0; imu2<NMU2+1; imu2++ ){
     logmu2 = MINLOGMU2 + imu2*DMU2;
        mu2 = exp(logmu2);
     
            save.precision(8);
             save<<setprecision(8)<<kt2<<"\t"<<mu2<<"\t"<<Tg(kt2,mu2)<<endl;
             cout.precision(8);
             // cout<<setprecision(8)<<kt2<<"\t"<<mu2<<"\t"<<endl;//Tg(kt2,mu2)<<endl;
        }
        cout<<setprecision(8)<<(double)(ikt2)/(NKT2)*100 <<"%\t time sigma TG \t" << (double)(clock()-t)/(CLOCKS_PER_SEC) << "sek" << endl;
    }
// }
    t1=clock()-t1;
    // cout<<"time TG: "<<((double)t1)/CLOCKS_PER_SEC<<endl;
    timetg=((double)t1)/CLOCKS_PER_SEC;
    save.close();
}

void draw_quarks_sudakov_factor()
{
    fstream save;
    string NAZWA="siatkaTQ_2_CZ_2";

        save.open(NAZWA,ios::out);
    clock_t t,t1;
    t=clock();
    t1 = clock();
     double MINX = 1.0e-6;
double MAXX = 0.99;
double MINKT2 = pow( 0.035, 2.0) ;
double MAXKT2 = pow( 1.0e4, 2.0 );
double MINMU2 = 1.79;
double MAXMU2 = pow( 1.0e4, 2.0 );
double MINLOGX, MINLOGKT2, MINLOGMU2;
double MAXLOGX, MAXLOGKT2, MAXLOGMU2;
int NX, NKT2, NMU2;
double DX, DKT2, DMU2;

MINLOGX = log( MINX );
MAXLOGX = log( MAXX );
MINLOGKT2 = log( MINKT2 );
MAXLOGKT2 = log( MAXKT2 );
MINLOGMU2 = log( MINMU2 );
MAXLOGMU2 = log( MAXMU2 );
NX = 60;
NKT2 = 140;
NMU2 = 120;
DX = ( MAXLOGX - MINLOGX ) / NX;
DKT2 = ( MAXLOGKT2 - MINLOGKT2 ) / NKT2;
DMU2 = ( MAXLOGMU2 - MINLOGMU2 ) / NMU2;
double logx;
double logkt2;
double logmu2;
double x;
double kt2;
double mu2;



//////////////////////////////////////////
/////// ATTENTION PLEASE!/////////////////
 MINLOGKT2=(MINLOGKT2+MINLOGKT2+DKT2)/2;//
//////////////////////////////////////////

// for( int ix=0; ix<NX+1; ix++ ){
//   logx = MINLOGX + ix*DX;
//   x = exp(logx);
  for( int ikt2=0; ikt2<NKT2+1; ikt2++ ){
    logkt2 = MINLOGKT2 + ikt2*DKT2;
        kt2 = exp(logkt2);
        t=clock();
    for( int imu2=0; imu2<NMU2+1; imu2++ ){
     logmu2 = MINLOGMU2 + imu2*DMU2;
        mu2 = exp(logmu2);
     
            save.precision(8);
            save<<setprecision(8)<<kt2<<"\t"<<mu2<<"\t"<<Tq(kt2,mu2)<<endl;
            // cout<<setprecision(8)<<kt2<<"\t"<<mu2<<"\t"<<Tq(kt2,mu2)<<endl;
        }
        cout<<setprecision(8)<<(double)(ikt2)/(NKT2)*100 <<"%\t time sigma TQ \t" << (double)(clock()-t)/(CLOCKS_PER_SEC) << "sek" << endl;
    }
// }
    t1=clock()-t1;
    // cout<<"time TG: "<<((double)t1)/CLOCKS_PER_SEC<<endl;
    timetq=((double)t1)/CLOCKS_PER_SEC;
    save.close();
}

void print_time()
{
    cout<<"time Tg \t"<<timetg<<endl;
    cout<<"time Tq \t"<<timetq<<endl;
}

double a(const double & x, const double & kt2)
{
    // return Ng*pow(x,-lambdag)*pow((1-x),Bg)*log(kt2/lambda0);
    return pdfs->xfxQ2(21,x,kt2)*log(kt2/lambda0); //PDF parton id (21=gluon)
    //Marcin nie mnozy przez log(kt2/lambda0) tylko od razu rozklad z biblioteki * Sudakov
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
    // return 2*Ns*pow(x,-lambdas)*pow((1-x),betas)*log(kt2/lambda0);
    return pdfs->xfxQ2(-2,x,kt2)*log(kt2/lambda0);
}

double di(const double & x, const double & kt2)
{
    // return 2*Ns*pow(x,-lambdas)*pow((1-x),betas)*log(kt2/lambda0);
    return pdfs->xfxQ2(-1,x,kt2)*log(kt2/lambda0);
}

double si(const double & x, const double & kt2)
{
    // return Ns*pow(x,-lambdas)*pow((1-x),betas)*log(kt2/lambda0);
    return pdfs->xfxQ2(-3,x,kt2)*log(kt2/lambda0);
}

double u(const double & x, const double & kt2)
{
    // return ui(x,kt2)+Nu*pow(x,alfau)*pow((1-x),betau)*log(kt2/lambda0);
    return pdfs->xfxQ2(2,x,kt2)*log(kt2/lambda0);
}

double d(const double & x, const double & kt2)
{
    // return di(x,kt2)+Nd*pow(x,alfad)*pow((1-x),betad)*log(kt2/lambda0);
    return pdfs->xfxQ2(1,x,kt2)*log(kt2/lambda0);
}

double ss(const double & x, const double & kt2)
{
    // return si(x,kt2);
    return pdfs->xfxQ2(3,x,kt2)*log(kt2/lambda0);
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