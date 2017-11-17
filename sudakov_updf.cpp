
#include "sudakov_g.h"
#include "sudakov_updf.h"
#include "Blad.h"


double h=0.1; //przy rozniczkowaniu numerycznym
 const LHAPDF::PDF* pdfs = LHAPDF::mkPDF("CT10nlo", 0); //wazne!

void draw_gluons()
{
    fstream save;
    string NAZWA="siatkaTG_2";

    save.open(NAZWA,ios::out);
    clock_t t;
    t=clock();
    // for(double & x : xs)
    // {
    for(double & val : kt2s)
        {
        for(double & mu2:mu2s)
            {
            save/*<<x*/<<"\t"<<val<<"\t"<<mu2<<"\t"<<Tg(val,mu2)<<endl;
            cout/*<<x*/<<"\t"<<val<<"\t"<<mu2<<"\t"<<Tg(val,mu2)<<endl;
            }
        }
        cout<<NAZWA/*<<"\t"<<x*/<<endl;
    // } 

    t=clock()-t;
    cout<<"time TG: "<<((double)t)/CLOCKS_PER_SEC<<endl;
    save.close();
}

void draw_quarks_sudakov_factor()
{
    fstream save;
    string NAZWA="siatkaTQ_2";

    save.open(NAZWA,ios::out);
    clock_t t;
    clock_t tt;
    t=clock();
    // for(double & x : xs)
    // {
    int licznik=0;
    for(double & val : kt2s)
        {                
        tt=clock();
        for(double & mu2:mu2s)
            {
            save/*<<x*/<<"\t"<<val<<"\t"<<mu2<<"\t"<<Tq(val,mu2)<<endl;
            cout<<licznik++<<"\t"<<val<<"\t"<<mu2<<"\t"<<Tq(val,mu2)<<endl;
            }
        tt=clock()-tt;
        cout<<"time sigma mu: "<<((double)tt)/CLOCKS_PER_SEC<<endl;
        }
        cout<<NAZWA/*<<"\t"<<x*/<<endl;
    // } 

    t=clock()-t;
    cout<<"time TQ: "<<((double)t)/CLOCKS_PER_SEC<<endl;
    save.close();
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