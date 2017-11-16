
#include "sudakov_g.h"
#include "sudakov_updf.h"
#include "Blad.h"


double h=0.001; //przy rozniczkowaniu numerycznym
const LHAPDF::PDF* pdfs = LHAPDF::mkPDF("CT10nlo", 0); //wazne!

double mu0 = 1.79;/*LHAPDF::mkPDF("CT10nlo", 0)->q2Min();*/

multimap<tuple<double,double>,double> TG;
multimap<tuple<double,double>,double> TQ;

vector<double> Mu;
vector<double> Kt2;

vector<double> Mured; //(Mureduced czyli z unikalnymi wartosciami X
vector<double> Kt2red;

vector<double> qMu;
vector<double> qKt2;

vector<double> qMured; //(Mureduced czyli z unikalnymi wartosciami X
vector<double> qKt2red;

double XMAX, XMIN, KT2MAX, KT2MIN, KTRANGE;
double Xindexgt;


void read_Ts()
{
/*    cout<<setprecision(8)<<pdfs->alphasQ2(2*2)<<endl;
    cout<<setprecision(8)<<pdfs->alphasQ(2)<<endl;
    cout<<pdfs->xfxQ2(1, 0.1, 2*2)<<endl;*/ //cos jak zlota regula Fermiego
    fstream f;
    double tmp1,tmp2,tmp3;
    f.open("siatkaTG", ios::in | ios::out);
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
    sort(Kt2.begin(),Kt2.end());
    vector<double>::iterator it=Kt2.begin();
    Kt2red.push_back(*it++);
    while(it!=Kt2.end())
    {
        if(*it!=*(it-1))
        {
                Kt2red.push_back(*it);
                // cout<<Kt2red.back()<<endl;
        }
        it++;
    }
    KTRANGE=KT2MAX-KT2MIN;
   
    XMAX=(*max_element(Mu.begin(),Mu.end()));
    XMIN=(*min_element(Mu.begin(),Mu.end()));
    sort(Mu.begin(),Mu.end());
    vector<double>::iterator it2=Mu.begin();
    Mured.push_back(*it2++);
    while(it2!=Mu.end())
    {
        if(*it2!=*(it2-1))
        {
                Mured.push_back(*it2);
                // cout<<Mured.back()<<endl;
        }
        it2++;
    }
    XMAX=(*max_element(Mured.begin(),Mured.end()));
    XMIN=(*min_element(Mured.begin(),Mured.end()));

    f.open("siatkaTQ", ios::in | ios::out);
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
    sort(qKt2.begin(),qKt2.end());
    vector<double>::iterator it3=qKt2.begin();
    qKt2red.push_back(*it++);
    while(it3!=qKt2.end())
    {
        if(*it3!=*(it3-1))
        {
                qKt2red.push_back(*it3);
                // cout<<Kt2red.back()<<endl;
        }
        it3++;
    }

    sort(qMu.begin(),qMu.end());
    vector<double>::iterator it21=qMu.begin();
    qMured.push_back(*it21++);
    while(it21!=qMu.end())
    {
        if(*it21!=*(it21-1))
        {
                qMured.push_back(*it21);
                // cout<<Mured.back()<<endl;
        }
        it21++;
    }
}


double Tqs(const double & kt2, const double & mu2)
{
    if(mu2<XMIN or mu2>XMAX)
        throw Blad("out of range w Tq",mu2,kt2,0.,0.);
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
        tab[0]=find_index_gt(qMured,mu2);
        tab[1]=tab[0]-1;
        tab[2]=find_index_gt(qKt2red,kt2);
        tab[3]=tab[2]-1;
        //tables of point around the choosen by the user
        double lf[] = {qMured[tab[1]], qKt2red[tab[3]], TQ.find(make_tuple(qKt2red[tab[1]], qMured[tab[3]]))->second};
        double rf[] = {qMured[tab[0]], qKt2red[tab[3]], TQ.find(make_tuple(qKt2red[tab[0]], qMured[tab[3]]))->second};
        double lr[] = {qMured[tab[1]], qKt2red[tab[2]], TQ.find(make_tuple(qKt2red[tab[1]], qMured[tab[2]]))->second};
        double rr[] = {qMured[tab[0]], qKt2red[tab[2]], TQ.find(make_tuple(qKt2red[tab[0]], qMured[tab[2]]))->second};
        //interpolation
        double a,b,c,d,f,g; //y=ax+b y1=cx+d y2=fx+g
        a = (lf[2]-rf[2])/(lf[0]-rf[0]);
        c = (lr[2]-rr[2])/(lr[0]-rr[0]);
        if(a==0 and c==0)
            return 0.0;
        b = lf[2]-a*lf[0];
        d = lr[2]-a*lr[0];
        double y = a*mu2+b;
        //double y1 = c*x+d;
        f = (y-(c*mu2+d))/(lf[1]-lr[1]);
        g = y-f*lf[1];
        double ret= f*kt2+g;
        if (ret<0)
            return 0.0;
        else
          //  cout<<difftime(stop,start)<<endl;
        return ret;
        }
}


double Tgs(const double & kt2, const double & mu2)
{
    if(mu2<XMIN or mu2>XMAX)
        throw Blad("out of range w fg",mu2,kt2,0.,0.);
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
        tab[0]=find_index_gt(Mured,mu2);
        tab[1]=tab[0]-1;
        tab[2]=find_index_gt(Kt2red,kt2);
        tab[3]=tab[2]-1;
        //tables of point around the choosen by the user
        double lf[] = {Mured[tab[1]], Kt2red[tab[3]], TG.find(make_tuple(Kt2red[tab[1]], Mured[tab[3]]))->second};
        double rf[] = {Mured[tab[0]], Kt2red[tab[3]], TG.find(make_tuple(Kt2red[tab[0]], Mured[tab[3]]))->second};
        double lr[] = {Mured[tab[1]], Kt2red[tab[2]], TG.find(make_tuple(Kt2red[tab[1]], Mured[tab[2]]))->second};
        double rr[] = {Mured[tab[0]], Kt2red[tab[2]], TG.find(make_tuple(Kt2red[tab[0]], Mured[tab[2]]))->second};
        //interpolation
        double a,b,c,d,f,g; //y=ax+b y1=cx+d y2=fx+g
        a = (lf[2]-rf[2])/(lf[0]-rf[0]);
        c = (lr[2]-rr[2])/(lr[0]-rr[0]);
        if(a==0 and c==0)
            return 0.0;
        b = lf[2]-a*lf[0];
        d = lr[2]-a*lr[0];
        double y = a*mu2+b;
        //double y1 = c*x+d;
        f = (y-(c*mu2+d))/(lf[1]-lr[1]);
        g = y-f*lf[1];
        double ret= f*kt2+g;
        if (ret<0)
            return 0.0;
        else
          //  cout<<difftime(stop,start)<<endl;
        return ret;
        }
}


void draw_gluons()
{
    /*vector<double> kt2s={
        0.001225,
0.00146581,
0.00175396,
0.00209875,
0.00251133,
0.00300501,
0.00359573,
0.00430258,
0.00514839,
0.00616046,
0.00737148,
0.00882057,
0.0105545,
0.0126293,
0.015112,
0.0180827,
0.0216375,
0.025891,
0.0309806,
0.0370708,
0.0443582,
0.0530781,
0.0635123,
0.0759975,
0.0909371,
0.108814,
0.130204,
0.1558,
0.186427,
0.223075,
0.266927,
0.3194,
0.382187,
0.457318,
0.547218,
0.65479,
0.783509,
0.937531,
1.12183,
1.34236,
1.60624,
1.922,
2.29983,
2.75193,
3.2929,
3.94022,
4.71479,
5.64163,
6.75067,
8.07771,
9.66563,
11.5657,
13.8393,
16.5598,
19.8152,
23.7104,
28.3715,
33.9487,
40.6224,
48.608,
58.1633,
69.5971,
83.2785,
99.6494,
119.239,
142.679,
170.726,
204.288,
244.447,
292.5,
350,
418.803,
501.132,
599.644,
717.523,
858.574,
1027.35,
1229.31,
1470.97,
1760.13,
2106.14,
2520.16,
3015.58,
3608.38,
4317.72,
5166.5,
6182.13,
7397.42,
8851.6,
10591.7,
12673.8,
15165.2,
18146.4,
21713.6,
25982,
31089.6,
37201.2,
44514.2,
53264.9,
63735.7,
76264.9,
91257.1,
109196,
130662,
156348,
187083,
223860,
267866,
320523,
383532,
458927,
549143,
657094,
786266,
940829,
1.12578e+06,
1.34708e+06,
1.61189e+06,
1.92876e+06,
2.30792e+06,
2.76161e+06,
3.30449e+06,
3.95409e+06,
4.73138e+06,
5.66148e+06,
6.77441e+06,
8.10613e+06,
9.69964e+06,
1.16064e+07,
1.3888e+07,
1.66181e+07,
1.98849e+07,
2.37939e+07,
2.84713e+07,
3.40682e+07,
4.07653e+07,
4.87789e+07,
5.83679e+07,
6.98419e+07,
8.35715e+07,
// 1e+08
    };
    vector<double> xs={
    
    5e-06,
6.38122e-06,
8.14398e-06,
1.03937e-05,
1.32649e-05,
1.69292e-05,
2.16058e-05,
2.75742e-05,
3.51914e-05,
4.49127e-05,
5.73196e-05,
7.31537e-05,
9.33619e-05,
0.000119152,
0.000152067,
0.000194075,
0.000247687,
0.000316108,
0.000403431,
0.000514876,
0.000657107,
0.000838628,
0.00107029,
0.00136595,
0.00174329,
0.00222486,
0.00283946,
0.00362384,
0.0046249,
0.0059025,
0.00753302,
0.00961396,
0.0122697,
0.0156592,
0.0199849,
0.0255056,
0.0325513,
0.0415434,
0.0530195,
0.0676657,
0.0863579,
0.110214,
0.140659,
0.179516,
0.229105,
0.292394,
0.373166,
0.47625,
0.607811,
0.775714,
0.99
    
    };
    vector<double> mu2s={
    1.79,
2.07689,
2.40975,
2.79596,
3.24408,
3.76401,
4.36727,
5.06722,
5.87935,
6.82164,
7.91496,
9.1835,
10.6553,
12.3631,
14.3445,
16.6436,
19.3111,
22.4061,
25.9971,
30.1637,
34.9981,
40.6073,
47.1154,
54.6667,
63.4282,
73.5939,
85.3889,
99.0743,
114.953,
133.377,
154.753,
179.556,
208.333,
241.723,
280.465,
325.415,
377.57,
438.083,
508.295,
589.761,
684.282,
793.953,
921.201,
1068.84,
1240.15,
1438.91,
1669.52,
1937.1,
2247.56,
2607.78,
3025.74,
3510.67,
4073.33,
4726.17,
5483.64,
6362.51,
7382.24,
8565.4,
9938.19,
11531,
13379.1,
15523.4,
18011.3,
20898,
24247.4,
28133.5,
32642.5,
37874.2,
43944.3,
50987.4,
59159.2,
68640.7,
79641.9,
92406.2,
107216,
124400,
144338,
167471,
194312,
225454,
261588,
303513,
352157,
408598,
474085,
550067,
638227,
740516,
859200,
996905,
1.15668e+06,
1.34206e+06,
1.55716e+06,
1.80672e+06,
2.09629e+06,
2.43227e+06,
2.82209e+06,
3.27439e+06,
3.79918e+06,
4.40808e+06,
5.11457e+06,
5.93428e+06,
6.88538e+06,
7.98891e+06,
9.26931e+06,
1.07549e+07,
1.24786e+07,
1.44786e+07,
1.67991e+07,
1.94915e+07,
2.26154e+07,
2.624e+07,
3.04455e+07,
3.53251e+07,
4.09867e+07,
4.75556e+07,
5.51775e+07,
6.40208e+07,
7.42815e+07,
8.61867e+07,
1e+08
        };
    */       
   cout<<mu0<<endl;
    string NAZWA = "siatki/GLUON_test_";
    stringstream stream;
    stream << fixed << setprecision(2) << "marcina" <<endl;
    string s = NAZWA + stream.str();
    cout<<s<<endl;
    fstream save;
    save.open(s,ios::out);
    if (!save.is_open()){ throw Blad("zly plik wejscia, nie istnieje lub zle wprowadzony");}
    clock_t t;
    t=clock();

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

for( int ix=0; ix<NX+1; ix++ ){
  double logx = MINLOGX + ix*DX;
  cout<< setprecision(8) <<logx<<"\t"<<exp(logx)<<endl;
  // for( int ikt2=0; ikt2<NKT2+1; ikt2++ ){
  //   double logkt2 = MINLOGKT2 + ikt2*DKT2;
  //   for( int imu2=0; imu2<NMU2+1; imu2++ ){
  //    double logmu2 = MINLOGMU2 + imu2*DMU2;
        // save<<logx<<"\t"<<logkt2<<"\t"<<logmu2<<"\t"<<fa(log(logx),log(logkt2),log(logmu2))<<endl;
        // cout<<logx<<"\t"<<logkt2<<"\t"<<logmu2<<endl;
        // }
    // }
}

double XTMP=exp(-13.815511);
double KT2TMP=exp(-6.7048144);

 cout<<Tg(KT2TMP, 0.58221562)<<endl;
 cout<<Tg(KT2TMP, 2.3660621)<<endl;
 cout<<Tg(KT2TMP, 4.1499086)<<endl;













    /*for(double & x : xs)
    {

    for(double & val : kt2s)
    {
        for(double & mu2 : mu2s)
        {
            
            save<<x<<"\t"<<val<<"\t"<<mu2<<"\t"<<fa(x,val,mu2)<<endl;
        }
    }
    cout<<s<<endl;
    t=clock()-t;
    cout<<"time sigma x: "<<((double)t)/CLOCKS_PER_SEC<<endl;
    } */
    save.close();


/*    string NAZWA="siatki/d/d_";    
    clock_t t;
    for(double & x : xs)
    {
    stringstream stream;
    stream << fixed << setprecision(2) << x <<endl;
    string s = NAZWA + stream.str();
    cout<<s<<endl;
    fstream save;
    save.open(s,ios::out);
    if (!save.is_open()){ throw Blad("zly plik wejscia, nie istnieje lub zle wprowadzony");}
    t=clock();

    for(double & val : kt2s)
    {
        for(double & mu2 : mu2s)
        {
            
            save<<x<<"\t"<<val<<"\t"<<mu2<<"\t"<<fd(x,val,mu2)<<endl;
        }
    }
    cout<<s<<endl;
    t=clock()-t;
    cout<<"time sigma x: "<<((double)t)/CLOCKS_PER_SEC<<endl;
    save.close();
    } 

    NAZWA="siatki/u/u_";    

    for(double & x : xs)
    {
    stringstream stream;
    stream << fixed << setprecision(2) << x <<endl;
    string s = NAZWA + stream.str();
    cout<<s<<endl;
    fstream save;
    save.open(s,ios::out);
    if (!save.is_open()){ throw Blad("zly plik wejscia, nie istnieje lub zle wprowadzony");}
    t=clock();

    for(double & val : kt2s)
    {
        for(double & mu2 : mu2s)
        {
            
            save<<x<<"\t"<<val<<"\t"<<mu2<<"\t"<<fu(x,val,mu2)<<endl;
        }
    }
    cout<<s<<endl;
    t=clock()-t;
    cout<<"time sigma x: "<<((double)t)/CLOCKS_PER_SEC<<endl;
    save.close();
    } 

 NAZWA="siatki/s/s_";    

    for(double & x : xs)
    {
    stringstream stream;
    stream << fixed << setprecision(2) << x <<endl;
    string s = NAZWA + stream.str();
    cout<<s<<endl;
    fstream save;
    save.open(s,ios::out);
    if (!save.is_open()){ throw Blad("zly plik wejscia, nie istnieje lub zle wprowadzony");}
    t=clock();

    for(double & val : kt2s)
    {
        for(double & mu2 : mu2s)
        {
            
            save<<x<<"\t"<<val<<"\t"<<mu2<<"\t"<<fs(x,val,mu2)<<endl;
        }
    }
    cout<<s<<endl;
    t=clock()-t;
    cout<<"time sigma x: "<<((double)t)/CLOCKS_PER_SEC<<endl;
    save.close();
    } 
    NAZWA="siatki/di/di_";    

    for(double & x : xs)
    {
    stringstream stream;
    stream << fixed << setprecision(2) << x <<endl;
    string s = NAZWA + stream.str();
    cout<<s<<endl;
    fstream save;
    save.open(s,ios::out);
    if (!save.is_open()){ throw Blad("zly plik wejscia, nie istnieje lub zle wprowadzony");}
    t=clock();

    for(double & val : kt2s)
    {
        for(double & mu2 : mu2s)
        {
            
            save<<x<<"\t"<<val<<"\t"<<mu2<<"\t"<<fdi(x,val,mu2)<<endl;
        }
    }
    cout<<s<<endl;
    t=clock()-t;
    cout<<"time sigma x: "<<((double)t)/CLOCKS_PER_SEC<<endl;
    save.close();
    } 

    NAZWA="siatki/ui/ui_";    


    for(double & x : xs)
    {
    stringstream stream;
    stream << fixed << setprecision(2) << x <<endl;
    string s = NAZWA + stream.str();
    cout<<s<<endl;
    fstream save;
    save.open(s,ios::out);
    if (!save.is_open()){ throw Blad("zly plik wejscia, nie istnieje lub zle wprowadzony");}
    t=clock();

    for(double & val : kt2s)
    {
        for(double & mu2 : mu2s)
        {
            
            save<<x<<"\t"<<val<<"\t"<<mu2<<"\t"<<fui(x,val,mu2)<<endl;
        }
    }
    cout<<s<<endl;
    t=clock()-t;
    cout<<"time sigma x: "<<((double)t)/CLOCKS_PER_SEC<<endl;
    save.close();
    } 

    NAZWA="siatki/si/si_";    


    for(double & x : xs)
    {
    stringstream stream;
    stream << fixed << setprecision(2) << x <<endl;
    string s = NAZWA + stream.str();
    cout<<s<<endl;
    fstream save;
    save.open(s,ios::out);
    if (!save.is_open()){ throw Blad("zly plik wejscia, nie istnieje lub zle wprowadzony");}
    t=clock();

    for(double & val : kt2s)
    {
        for(double & mu2 : mu2s)
        {
            
            save<<x<<"\t"<<val<<"\t"<<mu2<<"\t"<<fsi(x,val,mu2)<<endl;
        }
    }
    cout<<s<<endl;
    t=clock()-t;
    cout<<"time sigma x: "<<((double)t)/CLOCKS_PER_SEC<<endl;
    save.close();
    } 
*/
}

double a(const double & x, const double & kt2)
{
    // return Ng*pow(x,-lambdag)*pow((1-x),Bg)*log(kt2/lambda0);
    return pdfs->xfxQ2(21,x,kt2)*log(kt2/lambda0); //PDF parton id (21=gluon)
    //Marcin nie mnozy przez log(kt2/lambda0) tylko od razu rozklad z biblioteki * Sudakov
}

double fad(const double & x, const double & la2, const double & mu2)
{
    return pdfs->xfxQ2(21,x,la2)*Tgs(la2,mu2);
}

double fa(const double & x, const double & kt2, const double & mu2)
{
    if(sqrt(kt2)<sqrt(mu0))
	{
		return pdfs->xfxQ2(21,x,mu0)*Tqs(mu0,mu2);
	}
	else
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

double fu(const double & x, const double & kt2, const double & mu2)
{
    if(sqrt(kt2)<sqrt(mu0))
	{
		return pdfs->xfxQ2(2,x,mu0)*Tqs(mu0,mu2);
	}
	else
    return (fud(x,kt2+h,mu2)-fud(x,kt2-h,mu2))/(2*h);
}


double fd(const double & x, const double & kt2, const double & mu2)
{
    if(sqrt(kt2)<sqrt(mu0))
	{
		return pdfs->xfxQ2(1,x,mu0)*Tqs(mu0,mu2);
	}
	else
    return (fdd(x,kt2+h,mu2)-fdd(x,kt2-h,mu2))/(2*h);
}


double fs(const double & x, const double & kt2, const double & mu2)
{
	if(sqrt(kt2)<sqrt(mu0))
	{
		return pdfs->xfxQ2(3,x,mu0)*Tqs(mu0,mu2);
	}
	else
    return (fsd(x,kt2+h,mu2)-fsd(x,kt2-h,mu2))/(2*h);
}


double fui(const double & x, const double & kt2, const double & mu2)
{
    if(sqrt(kt2)<sqrt(mu0))
	{
		return pdfs->xfxQ2(-2,x,mu0)*Tqs(mu0,mu2);
	}
	else
    return (fuid(x,kt2+h,mu2)-fuid(x,kt2-h,mu2))/(2*h);
}


double fdi(const double & x, const double & kt2, const double & mu2)
{
    if(sqrt(kt2)<sqrt(mu0))
	{
		return pdfs->xfxQ2(-1,x,mu0)*Tqs(mu0,mu2);
	}
	else
    return (fdid(x,kt2+h,mu2)-fdid(x,kt2-h,mu2))/(2*h);
}


double fsi(const double & x, const double & kt2, const double & mu2)
{
	double ret;
	if(sqrt(kt2)<sqrt(mu0))
	{
		ret= pdfs->xfxQ2(-3,x,mu0)*Tqs(mu0,mu2);
	}
	else
	{
    ret= (fsid(x,kt2+h,mu2)-fsid(x,kt2-h,mu2))/(2*h);
	}
	return ret;
}