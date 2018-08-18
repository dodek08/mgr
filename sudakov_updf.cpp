
#include "sudakov_g.h"
#include "sudakov_updf.h"
#include "Blad.h"


#define FIT

double h=0.1; //przy rozniczkowaniu numerycznym
LHAPDF::PDF* pdfs;// = LHAPDF::mkPDF("CT14nlo", 0); //wazne!

double mu0;// = LHAPDF::mkPDF("CT14nlo", 0)->q2Min()+0.1;

multimap<tuple<double,double>,double> TG;
multimap<tuple<double,double>,double> TQ;
multimap<tuple<double,double,double>,double> ccfm;

vector<double> ccfm_X;
vector<double> ccfm_Mu;
vector<double> ccfm_Kt2;

vector<double> ccfm_Xred;
vector<double> ccfm_Mured; //(Mureduced czyli z unikalnymi wartosciami
vector<double> ccfm_Kt2red;

vector<double> Mu;
vector<double> Kt2;

vector<double> Mured; //(Mureduced czyli z unikalnymi wartosciami
vector<double> Kt2red;

vector<double> qMu;
vector<double> qKt2;

vector<double> qMured; //(Mureduced czyli z unikalnymi wartosciami
vector<double> qKt2red;

double MU2MAX, MU2MIN, KT2MAX, KT2MIN, KTRANGE;
double Xindexgt;

double beta_g = 0;
double n_g = 1;
double lambda_g = 0;

void set_beta_g(const double & val)
{
    beta_g = val;
}

double get_beta_g()
{
    return beta_g;
}


void set_n_g(const double & val)
{
    n_g = val;
}

double get_n_g()
{
    return n_g;
}

void set_lambda_g(const double & val)
{
    lambda_g = val;
}

double get_lambda_g()
{
    return lambda_g;
}

void set_pdf_name_sudakov_updf(string siatka)
{
 pdfs = LHAPDF::mkPDF(siatka, 0);
 mu0 =  LHAPDF::mkPDF(siatka, 0)->q2Min()+0.1;
}

double kt2max()
{
	return KT2MAX;
}

double kt2min()
{
	return KT2MIN;
}

double mu2max()
{
	return MU2MAX;
}

double mu2min()
{
	return MU2MIN;
}


void read_ccfm()
{
/*    cout<<setprecision(8)<<pdfs->alphasQ2(2*2)<<endl;
    cout<<setprecision(8)<<pdfs->alphasQ(2)<<endl;
    cout<<pdfs->xfxQ2(1, 0.1, 2*2)<<endl;*/ //cos jak zlota regula Fermiego
    fstream f;
    double tmp0,tmp1,tmp2,tmp3;
    f.open("ccfm.dat", ios::in | ios::out);
    if (!f.is_open()){ throw Blad("zly plik wejscia ccfm, nie istnieje lub zle wprowadzony");}

    while(!f.eof())
    {
        f>>tmp0; //x
        f>>tmp1; //kt2
        f>>tmp2;    //mu2
        f>>tmp3;    //TG
        ccfm.insert(pair<tuple<double,double,double>,double>(make_tuple(exp(tmp0),exp(tmp1),exp(tmp2)*exp(tmp2)),tmp3));

        ccfm_X.push_back(exp(tmp0));
        ccfm_Kt2.push_back(exp(tmp1));
        ccfm_Mu.push_back(exp(tmp2)*exp(tmp2));
    }
    f.close();

    ccfm_Kt2red=subset_with_sort(ccfm_Kt2);
    ccfm_Xred=subset_with_sort(ccfm_X);
    ccfm_Mured=subset_with_sort(ccfm_Mu);
    cout<<ccfm_Xred<<endl;
    cout<<endl;
    cout<<ccfm_Kt2red<<endl;
    cout<<endl;
    cout<<ccfm_Mured<<endl;
}



double ccfm_s(const double & x, const double & kt2, const double & mu2)
{
    // if(mu2<MU2MIN or mu2>MU2MAX)
        // throw Blad("out of range w Tq",mu2,kt2,0.,0.);
    int tab[3]={0}; //table of indexes
    //unsigned int i=1;

    multimap<tuple<double,double>,double,double>::iterator it=TQ.find(make_tuple(kt2,mu2));

    if(it!=TQ.end() or (qMu.back()==mu2 and qKt2.back()==kt2))
    {
        if (it->second<0)
            return 0.0;
        else
            return it->second;
    }
    else
    {
        double xl, xh;
        int x_ind = find_index_gt(ccfm_Xred, x);
        xl = ccfm_Xred[x_ind-1];
        xh = ccfm_Xred[x_ind];
        double kth, ktl, muh, mul;
        tab[0]=find_index_gt(ccfm_Mured,mu2);
        tab[1]=tab[0]-1;
        tab[2]=find_index_gt(ccfm_Kt2red,kt2);
        tab[3]=tab[2]-1;
        muh = ccfm_Mured[tab[0]];
        mul = ccfm_Mured[tab[1]];
        kth = ccfm_Kt2red[tab[2]];
        ktl = ccfm_Kt2red[tab[3]];
        //punkt pierwszy
        double lf[] = {mul, ktl, ccfm.find(make_tuple(xl, ktl, mul))->second};
        double rf[] = {muh, ktl, ccfm.find(make_tuple(xl, ktl, muh))->second};
        double lr[] = {mul, kth, ccfm.find(make_tuple(xl, kth, mul))->second};
        double rr[] = {muh, kth, ccfm.find(make_tuple(xl, kth, muh))->second};//interpolation
        double a,b,c,d,f,g; //y=ax+b y1=cx+d y2=fx+g
        a = (rf[2]-rr[2])/(rf[1]-rr[1]);
        c = (lf[2]-lr[2])/(rf[1]-rr[1]);
        b = rf[2]-a*rf[1];
        d = lf[2]-c*lf[1];
        double y = a*kt2+b;
        double y1 = c*kt2+d;
        f = (y-y1)/(rf[0]-lf[0]);
        g = y1-f*lf[0];
        double ret1= f*mu2+g;
        //punkt drugi
        double lf2[] = {mul, ktl, ccfm.find(make_tuple(xh, ktl, mul))->second};
        double rf2[] = {muh, ktl, ccfm.find(make_tuple(xh, ktl, muh))->second};
        double lr2[] = {mul, kth, ccfm.find(make_tuple(xh, kth, mul))->second};
        double rr2[] = {muh, kth, ccfm.find(make_tuple(xh, kth, muh))->second};//interpolation
        //y=ax+b y1=cx+d y2=fx+g
        a = (rf2[2]-rr2[2])/(rf2[1]-rr2[1]);
        c = (lf2[2]-lr2[2])/(rf2[1]-rr2[1]);
        b = rf2[2]-a*rf2[1];
        d = lf2[2]-c*lf2[1];
        y = a*kt2+b;
        y1 = c*kt2+d;
        f = (y-y1)/(rf2[0]-lf2[0]);
        g = y1-f*lf2[0];
        double ret2= f*mu2+g;
        //interpolacja miedzy xl a xh
        a = (ret2 - ret1)/(xh - xl);
        b = ret1 - xl*a;
        double ret = a*x + b;

        if (ret<0)
            return 0.0;
        else
          //  cout<<difftime(stop,start)<<endl;
        return ret;
        }
}


void read_Ts()
{
/*    cout<<setprecision(8)<<pdfs->alphasQ2(2*2)<<endl;
    cout<<setprecision(8)<<pdfs->alphasQ(2)<<endl;
    cout<<pdfs->xfxQ2(1, 0.1, 2*2)<<endl;*/ //cos jak zlota regula Fermiego
    fstream f;
    double tmp1,tmp2,tmp3;
    f.open("COMPLETE_TG_", ios::in | ios::out);
    if (!f.is_open()){ throw Blad("zly plik wejscia, nie istnieje lub zle wprowadzony");}

    while(!f.eof())
    {
        f>>tmp1; //kt2
        f>>tmp2;    //mu2
        f>>tmp3;    //TG
        TG.insert(pair<tuple<double,double>,double>(make_tuple(tmp1,tmp2),tmp3));

        Kt2.push_back(tmp1);
        Mu.push_back(tmp2);
    }
    f.close();

    KT2MAX=(*max_element(Kt2.begin(),Kt2.end()));
    KT2MIN=(*min_element(Kt2.begin(),Kt2.end()));
    Kt2red=subset_with_sort(Kt2);
    KTRANGE=KT2MAX-KT2MIN;
    Mured=subset_with_sort(Mu);
    MU2MAX=(*max_element(Mured.begin(),Mured.end()));
    MU2MIN=(*min_element(Mured.begin(),Mured.end()));

    f.open("COMPLETE_TQ_", ios::in | ios::out);
    if (!f.is_open()){ throw Blad("zly plik wejscia, nie istnieje lub zle wprowadzony");}

    while(!f.eof())
    {
        f>>tmp1; //kt2
        f>>tmp2;    //mu2
        f>>tmp3;    //TQ
        TQ.insert(pair<tuple<double,double>,double>(make_tuple(tmp1,tmp2),tmp3));

        qKt2.push_back(tmp1);
        qMu.push_back(tmp2);
    }
    f.close();

    qKt2red=subset_with_sort(qKt2);
    qMured=subset_with_sort(qMu);
}


double Tqs(const double & kt2, const double & mu2)
{
    // if(mu2<MU2MIN or mu2>MU2MAX)
        // throw Blad("out of range w Tq",mu2,kt2,0.,0.);
    int tab[3]={0}; //table of indexes
    //unsigned int i=1;

    multimap<tuple<double,double>,double,double>::iterator it=TQ.find(make_tuple(kt2,mu2));

    if(it!=TQ.end() or (qMu.back()==mu2 and qKt2.back()==kt2))
    {
        if (it->second<0)
            return 0.0;
        else
            return it->second;
    }
    else
    {
        double kth, ktl, muh, mul;
        tab[0]=find_index_gt(qMured,mu2);
        tab[1]=tab[0]-1;
        tab[2]=find_index_gt(qKt2red,kt2);
        tab[3]=tab[2]-1;
        muh = qMured[tab[0]];
        mul = qMured[tab[1]];
        kth = qKt2red[tab[2]];
        ktl = qKt2red[tab[3]];
        //tables of point around the choosen by the user
        double lf[] = {mul, ktl, TQ.find(make_tuple(ktl, mul))->second};
        double rf[] = {muh, ktl, TQ.find(make_tuple(ktl, muh))->second};
        double lr[] = {mul, kth, TQ.find(make_tuple(kth, mul))->second};
        double rr[] = {muh, kth, TQ.find(make_tuple(kth, muh))->second};//interpolation
        double a,b,c,d,f,g; //y=ax+b y1=cx+d y2=fx+g
        a = (rf[2]-rr[2])/(rf[1]-rr[1]);
        c = (lf[2]-lr[2])/(rf[1]-rr[1]);
        b = rf[2]-a*rf[1];
        d = lf[2]-c*lf[1];
        double y = a*kt2+b;
        double y1 = c*kt2+d;
        f = (y-y1)/(rf[0]-lf[0]);
        g = y1-f*lf[0];
        double ret= f*mu2+g;
        if (ret<0)
            return 0.0;
        else
          //  cout<<difftime(stop,start)<<endl;
        return ret;
        }
}


double Tgs(const double & kt2, const double & mu2)
{
    // if(mu2<MU2MIN or mu2>MU2MAX)
    //     throw Blad("out of range w Tg",mu2,kt2,MU2MIN,MU2MAX);
    int tab[3]={0}; //table of indexes
    //unsigned int i=1;

    multimap<tuple<double,double>,double,double>::iterator it=TG.find(make_tuple(kt2,mu2));

    if(it!=TG.end() or (Mu.back()==mu2 and Kt2.back()==kt2))
    {
        if (it->second<0)
            return 0.0;
        else
            return it->second;
    }
    else
    {

        double kth, ktl, muh, mul;
        tab[0]=find_index_gt(Mured,mu2);
        tab[1]=tab[0]-1;
        tab[2]=find_index_gt(Kt2red,kt2);
        tab[3]=tab[2]-1;
        muh = Mured[tab[0]];
        mul = Mured[tab[1]];
        kth = Kt2red[tab[2]];
        ktl = Kt2red[tab[3]];
        //tables of point around the choosen by the user
        double lf[] = {mul, ktl, TG.find(make_tuple(ktl, mul))->second};
        double rf[] = {muh, ktl, TG.find(make_tuple(ktl, muh))->second};
        double lr[] = {mul, kth, TG.find(make_tuple(kth, mul))->second};
        double rr[] = {muh, kth, TG.find(make_tuple(kth, muh))->second};

        // cout<<"LF \t"<<mul<<"\t"<<ktl<<"\t"<<TG.find(make_tuple(ktl, mul))->second<<"\t"<<Tg(ktl,mul)<<endl;
        // cout<<"RF \t"<<muh<<"\t"<<ktl<<"\t"<<TG.find(make_tuple(ktl, muh))->second<<"\t"<<Tg(ktl,muh)<<endl;
        // cout<<"LR \t"<<mul<<"\t"<<kth<<"\t"<<TG.find(make_tuple(kth, mul))->second<<"\t"<<Tg(kth,mul)<<endl;
        // cout<<"RR \t"<<muh<<"\t"<<kth<<"\t"<<TG.find(make_tuple(kth, muh))->second<<"\t"<<Tg(kth,muh)<<endl;
        //interpolation
        double a,b,c,d,f,g; //y=ax+b y1=cx+d y2=fx+g
        a = (rf[2]-rr[2])/(rf[1]-rr[1]);
        c = (lf[2]-lr[2])/(rf[1]-rr[1]);
        b = rf[2]-a*rf[1];
        d = lf[2]-c*lf[1];
        double y = a*kt2+b;
        double y1 = c*kt2+d;
        f = (y-y1)/(rf[0]-lf[0]);
        g = y1-f*lf[0];
        double ret= f*mu2+g;
        if (ret<0)
            return 0.0;
        else
          //  cout<<difftime(stop,start)<<endl;
        return ret;
        }
}

ostream& operator<<(ostream& os, const vector<double>& v)
{
    for (double val : v)
    {
        os<<val<<"\n";
    }
    return os;
}

vector<double> subset_with_sort(vector<double>& v)
{
    vector<double> ret;
    sort(v.begin(),v.end());
    vector<double>::iterator it=v.begin();
    ret.push_back(*it++);
    while(it!=v.end())
    {
        if(*it!=*(it-1))
        {
                ret.push_back(*it);
                // cout<<Kt2red.back()<<endl;
        }
        it++;
    }

    return ret;
}

void draw_gluons(const int & pid)
{
    
    
   // cout<<mu0<<endl;
    string NAZWA = "siatki/rozklad_CT10nlo_pid_" + to_string(pid);
    fstream save;
    save.open(NAZWA,ios::out);
    if (!save.is_open()){ throw Blad("zly plik wejscia, nie istnieje lub zle wprowadzony");}

//kod od marcina generujacy siatke

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
        clock_t t = clock();

if(pid != 21)
{
    for( int ix=0; ix<NX+1; ix++ ){
  logx = MINLOGX + ix*DX;
  x = exp(logx);
  for( int ikt2=0; ikt2<NKT2+1; ikt2++ ){
    logkt2 = MINLOGKT2 + ikt2*DKT2;
        kt2 = exp(logkt2);
        // t=clock();
    for( int imu2=0; imu2<NMU2+1; imu2++ ){
     logmu2 = MINLOGMU2 + imu2*DMU2;
        mu2 = exp(logmu2);
     
            save<<x<<"\t"<<kt2<<"\t"<<mu2<<"\t"<<fq_pid(x,kt2,mu2,pid)<<endl;
        }
        // cout << "time sigma kt2 " << (double)(clock()-t)/(CLOCKS_PER_SEC) << "sek" << endl;
    }
}
}
else
{
    vector<double> xss = {0.00001,0.0001,0.001,0.01,0.1};
    vector<double> mu2s = {10,100,1000,100000,1000000};
    for( double x : xss)
    {
  for( int ikt2=0; ikt2<NKT2+1; ikt2++ )
  {
    logkt2 = MINLOGKT2 + ikt2*DKT2;
        kt2 = exp(logkt2);
        // t=clock();
    for(double mu2 : mu2s)
    {
     
            save<<x<<"\t"<<kt2<<"\t"<<mu2<<"\t"<<fa(x,kt2,mu2)<<endl;
    }
    }
    }
// for( int ix=0; ix<NX+1; ix++ ){
//   logx = MINLOGX + ix*DX;
//   x = exp(logx);
//   for( int ikt2=0; ikt2<NKT2+1; ikt2++ ){
//     logkt2 = MINLOGKT2 + ikt2*DKT2;
//         kt2 = exp(logkt2);
//         // t=clock();
//     for( int imu2=0; imu2<NMU2+1; imu2++ ){
//      logmu2 = MINLOGMU2 + imu2*DMU2;
//         mu2 = exp(logmu2);
     
//             save<<x<<"\t"<<kt2<<"\t"<<mu2<<"\t"<<fa(x,kt2,mu2)<<endl;
//         }
        // cout << "time sigma kt2 " << (double)(clock()-t)/(CLOCKS_PER_SEC) << "sek" << endl;
    }


    save.close();
     cout << "time sigma \t"<<pid<<"\t" << (double)(clock()-t)/(CLOCKS_PER_SEC) << "sek" << endl;

}

double q_pid(const double & x, const double & kt2, const int & pid)
{
    return pdfs->xfxQ2(pid,x,kt2);
}

double fq_pid_d(const double & x, const double & la2, const double & mu2, const int & pid)
{
    // if (pid == 21)
    // {
    //     return pdfs->xfxQ2(pid,x,la2)*Tgs(la2,mu2);
    // }
    // else
    return pdfs->xfxQ2(pid,x,la2)*Tqs(la2,mu2);
}

double fq_pid(const double & x, const double & kt2, const double & mu2, const int & pid)
{
    if(sqrt(kt2)<sqrt(mu0))
    {
        // if(pid == 21)
        //     return pdfs->xfxQ2(pid,x,mu0)*Tgs(mu0,mu2)/mu0;
        // else
        return pdfs->xfxQ2(pid,x,mu0)*Tqs(mu0,mu2)/mu0;
    }
    else
    return (fq_pid_d(x,kt2+h,mu2,pid)-fq_pid_d(x,kt2-h,mu2,pid))/(2*h);
}

double a(const double & x, const double & kt2)
{
    // return Ng*pow(x,-lambdag)*pow((1-x),Bg)*log(kt2/lambda0);
    return pdfs->xfxQ2(21,x,kt2);//*log(kt2/lambda0); //PDF parton id (21=gluon)
    //Marcin nie mnozy przez log(kt2/lambda0) tylko od razu rozklad z biblioteki * Sudakov
}

double fad(const double & x, const double & la2, const double & mu2)
{
    return pdfs->xfxQ2(21,x,la2)*Tgs(la2,mu2);
}

double ui(const double & x, const double & kt2)
{
    // return 2*Ns*pow(x,-lambdas)*pow((1-x),betas)*log(kt2/lambda0);
    return pdfs->xfxQ2(-2,x,kt2);//*log(kt2/lambda0);
}

double di(const double & x, const double & kt2)
{
    // return 2*Ns*pow(x,-lambdas)*pow((1-x),betas)*log(kt2/lambda0);
    return pdfs->xfxQ2(-1,x,kt2);//*log(kt2/lambda0);
}

double si(const double & x, const double & kt2)
{
    // return Ns*pow(x,-lambdas)*pow((1-x),betas)*log(kt2/lambda0);
    return pdfs->xfxQ2(-3,x,kt2);//*log(kt2/lambda0);
}

double u(const double & x, const double & kt2)
{
    // return ui(x,kt2)+Nu*pow(x,alfau)*pow((1-x),betau)*log(kt2/lambda0);
    return pdfs->xfxQ2(2,x,kt2);//*log(kt2/lambda0);
}

double d(const double & x, const double & kt2)
{
    // return di(x,kt2)+Nd*pow(x,alfad)*pow((1-x),betad)*log(kt2/lambda0);
    return pdfs->xfxQ2(1,x,kt2);//*log(kt2/lambda0);
}

double ss(const double & x, const double & kt2)
{
    // return si(x,kt2);
    return pdfs->xfxQ2(3,x,kt2);//*log(kt2/lambda0);
}

double fuid(const double & x, const double & la2, const double & mu2)
{
    return pdfs->xfxQ2(-2,x,la2)*Tqs(la2,mu2);
}

double fdid(const double & x, const double & la2, const double & mu2)
{
    return pdfs->xfxQ2(-1,x,la2)*Tqs(la2,mu2);
}

double fsid(const double & x, const double & la2, const double & mu2)
{
    return pdfs->xfxQ2(-3,x,la2)*Tqs(la2,mu2);
}

double fdd(const double & x, const double & la2, const double & mu2)
{
    return pdfs->xfxQ2(1,x,la2)*Tqs(la2,mu2);
}

double fud(const double & x, const double & la2, const double & mu2)
{
    return pdfs->xfxQ2(2,x,la2)*Tqs(la2,mu2);
}

double fsd(const double & x, const double & la2, const double & mu2)
{
    return pdfs->xfxQ2(3,x,la2)*Tqs(la2,mu2);
}

double fa(const double & x, const double & kt2, const double & mu2)
{
    if(sqrt(kt2)<sqrt(mu0))
    {
        return pdfs->xfxQ2(21,x,mu0)*Tgs(mu0,mu2)/mu0;
    }
    else
    #ifdef FIT
    return (n_g*pow(x,lambda_g)*pow((1.-x),beta_g))*(fad(x,kt2+h,mu2)-fad(x,kt2-h,mu2))/(2*h);
    #else
    return (fad(x,kt2+h,mu2)-fad(x,kt2-h,mu2))/(2*h);
    #endif

}

double fu(const double & x, const double & kt2, const double & mu2)
{
    if(sqrt(kt2)<sqrt(mu0))
	{
		return pdfs->xfxQ2(2,x,mu0)*Tqs(mu0,mu2)/mu0;
	}
	else
    return (fud(x,kt2+h,mu2)-fud(x,kt2-h,mu2))/(2*h);
}


double fd(const double & x, const double & kt2, const double & mu2)
{
    if(sqrt(kt2)<sqrt(mu0))
	{
		return pdfs->xfxQ2(1,x,mu0)*Tqs(mu0,mu2)/mu0;
	}
	else
    return (fdd(x,kt2+h,mu2)-fdd(x,kt2-h,mu2))/(2*h);
}


double fs(const double & x, const double & kt2, const double & mu2)
{
	if(sqrt(kt2)<sqrt(mu0))
	{
		return pdfs->xfxQ2(3,x,mu0)*Tqs(mu0,mu2)/mu0;
	}
	else
    return (fsd(x,kt2+h,mu2)-fsd(x,kt2-h,mu2))/(2*h);
}


double fui(const double & x, const double & kt2, const double & mu2)
{
    if(sqrt(kt2)<sqrt(mu0))
	{
		return pdfs->xfxQ2(-2,x,mu0)*Tqs(mu0,mu2)/mu0;
	}
	else
    return (fuid(x,kt2+h,mu2)-fuid(x,kt2-h,mu2))/(2*h);
}


double fdi(const double & x, const double & kt2, const double & mu2)
{
    if(sqrt(kt2)<sqrt(mu0))
	{
		return pdfs->xfxQ2(-1,x,mu0)*Tqs(mu0,mu2)/mu0;
	}
	else
    return (fdid(x,kt2+h,mu2)-fdid(x,kt2-h,mu2))/(2*h);
}


double fsi(const double & x, const double & kt2, const double & mu2)
{
	double ret;
	if(sqrt(kt2)<sqrt(mu0))
	{
		ret= pdfs->xfxQ2(-3,x,mu0)*Tqs(mu0,mu2)/mu0;
	}
	else
	{
    ret= (fsid(x,kt2+h,mu2)-fsid(x,kt2-h,mu2))/(2*h);
	}
	return ret;
}

