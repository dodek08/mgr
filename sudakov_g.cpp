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
gsl_monte_vegas_state *s;

vector<double> u2s; //dostepne wartosci u2 z pliku
vector<double> Tgs; //wektor wyliczonych Tg
vector<double> kt2su; //wektor wykorzystanych kt2
vector<double> u2su; //wektor wykorzystanych u2
map<double,double> as_2;//mapa przechowujaca as_2(u^2)
size_t calls = 1000000; //dla beki liczba iteracji
int n=100000000; //liczba iteracji
double h=0.1;

//ftg na dole bo sie pluje


double Tq(const double &, const double &)
{
    return 1.;
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
        double ret=0;
        double ret1=0;
        double ret2=0;
        double z=0;
        double delta=0;
        double tmp=0;
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

double nTg(const double & kt2, const double & u2)
{
    if(kt2>u2)
        throw Blad("kt^2 wieksze od u^2", kt2, u2, kt2, u2);
    vector<double>::iterator test=find(u2s.begin(), u2s.end(), u2);
    if((*test!=u2))
        throw Blad("nie ma alfy dla takiego u^2", kt2, u2, kt2, u2);
    else
    {
        double ret=0;
        double ret1=0;
        double ret2=0;
        double z=0;
        double delta=0;
        double tmp=0;
        vector<double>::iterator it=find(u2s.begin(), u2s.end(), kt2);
        if((*test!=u2))
            throw Blad("nie ma takiego kt^2", kt2, u2, kt2, u2);
        kt2su.push_back(kt2);
        u2su.push_back(u2);

        for(int i=0; i<n; i=1)


        Tgs.push_back(exp(-ret));
        return Tgs.back();
    }
}


double Tg(const double & kt2, const double & u2)
{
	  if(kt2>u2)
	  return 1.;
    vector<double>::iterator test=find(u2s.begin(), u2s.end(), u2);
    gsl_rng_set(r, chrono::system_clock::now().time_since_epoch().count());
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
//            cout<<"chisq<< "<<s->chisq<<endl;
            //s->stage=1;
            if(s->chisq==0)
                //gsl_monte_vegas_init(s);
                break;
        }
        while ((fabs (s->chisq - 1.0) > 0.35) ); //w celu uzyskania najwiekszej dokladnosci
    gsl_monte_vegas_free(s);

  t=clock()-t;
  cout<<"time TG: "<<((double)t)/CLOCKS_PER_SEC<<endl;
//  cout<<result<<endl;
  return exp(-result);
}


void warm_up()
{
  gsl_rng_env_setup ();
  H.f=&ftg;
  H.dim=2;
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
        //cout<<"Heaviside"<<endl;
            return 0.;

        }
    else
    {
    //cout<<"Heaviside"<<endl;
        return 0.;

    }
}

inline double Pgg(const double & z)
{
    double ret = 6*(z*(1-z)+(1-z)/z+z/(1-z));
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
        //as_2.insert(pair<double,double>(tmp1,tmp2));
        as_2[tmp1]=tmp2;
        //cout<<u2s.back()<<" "<<as_2[tmp1]<<endl;
    }
    f.close();
    sort(u2s.begin(), u2s.end());
}

void make_lattice()
{
    double tmp=0;
    vector<double>::iterator U2=u2s.begin();
    vector<double>::iterator KT2=u2s.begin();
    while(U2!=u2s.end())
    {
        while(KT2!=U2)
        {
            tmp=Tg(*KT2++,*U2);
        }
        U2++;
        KT2=u2s.begin();
    }
    vector<double>::iterator it1=kt2su.begin();
    vector<double>::iterator it2=u2su.begin();
    vector<double>::iterator it3=Tgs.begin();
    fstream save;
    save.open("wyniki.dat", ios::out );
    while(it3!=Tgs.end())
    {
            save<<*it1++<<"\t"<<*it2++<<"\t"<<*it3++<<endl;
    }
    save.close();
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
//    int index = find_index_gt(u2s,kt2);
//    if(index<2)
//        index=2;
//    double x2=u2s[index];
//    double x3=u2s[index+1];
//    double x0=u2s[index-2];
//    double x1=u2s[index-1];
//    double y0=as_2[x0];
//    double y1=as_2[x1];
//    double y2=as_2[x2];
//    double y3=as_2[x3];
//    return y0*(kt2-x1)/(x0-x1)*(kt2-x2)/(x0-x2)*(kt2-x3)/(x0-x3)+y1*(kt2-x0)/(x1-x0)*(kt2-x2)/(x1-x2)*(kt2-x3)/(x1-x3)+y2*(kt2-x1)/(x2-x1)*(kt2-x1)/(x2-x1)*(kt2-x3)/(x2-x3)+y3*(kt2-x0)/(x3-x0)*(kt2-x2)/(x3-x2)*(kt2-x1)/(x3-x1);

}

double Tg2(const double & kt2, const double & u2)
{
    if(kt2>u2)
        throw Blad("kt^2 wieksze od u^2", kt2, u2, kt2, u2);
    vector<double>::iterator test=find(u2s.begin(), u2s.end(), u2);
    if((*test!=u2))
        throw Blad("nie ma alfy dla takiego u^2", kt2, u2, kt2, u2);
    else
    {
        double ret=0;
        double z=0;
        double delta=0;
        double pt2;
        vector<double>::iterator it=find(u2s.begin(), u2s.end(), kt2);
        if((*test!=u2))
            throw Blad("nie ma takiego kt^2", kt2, u2, kt2, u2);
        kt2su.push_back(kt2);
        u2su.push_back(u2);

        for(int i=0; i<n; i++)
        {
            pt2=(generatorFBN()*1.0/(RMAX*1.0))*(u2-kt2)+kt2;
            //cout<<pt2<<"\t"<<interpolacja(pt2)<<endl;
            delta=sqrt(pt2)/(sqrt(u2)+sqrt(pt2));
            for(int j=0; j<100; j++)
            {

                z=(generatorFBN()*1.0/(RMAX*1.0))*(1.-2.*delta)+delta;
            ret+=interpolacja(pt2)/(pt2*2*M_PI)*totPgg(z)*theta(delta,z)+Pqg(z);
            }
            //cout<<(u2-kt2)*(1-2*delta)/(1.*n)*ret<<endl;
        }
        ret=(u2-kt2)/(100.*n)*ret;
        Tgs.push_back(exp(-ret));
        return Tgs.back();
    }
}

///
double calka_test(const double & lx, const double & ux, const double & ly, const double & uy)
{

        double ret2=0;
        double ret=0;
        double z=0;
        double x=0;
        for(int i=0; i<n; i++)
        {
        z=(generatorFBN()*1.0/(RMAX*1.0))*(uy-ly)+ly;
        x=(generatorFBN()*1.0/(RMAX*1.0))*(ux-lx)+lx;
        ret2+=x*(z*z+z); //jeden w pamieci z calki po Pqg
       // cout<<ret2<<endl;
        }
        ret=ret2*((ux-lx)*(uy-ly))/(n*1.); //wartosc calki po dy
                //cout<<ret<<endl;
        return ret;
}


double Tgt(const double & kt2, const double & u2)
{
    if(kt2>u2)
        throw Blad("kt^2 wieksze od u^2", kt2, u2, kt2, u2);
    vector<double>::iterator test=find(u2s.begin(), u2s.end(), kt2);
    if((*test!=kt2))
            throw Blad("nie ma takiego kt^2", kt2, u2, kt2, u2);
    test=find(u2s.begin(), u2s.end(), u2);
    if((*test!=u2))
        throw Blad("nie ma alfy dla takiego u^2", kt2, u2, kt2, u2);
    else
    {
        test--;
        double ret=0;
        double ret1=0;
        double ret2=0;
        double z=0;
        double delta=0;
        double tmp=0;
        vector<double>::iterator it=find(u2s.begin(), u2s.end(), kt2);
        kt2su.push_back(kt2);
        u2su.push_back(u2);
        while(*it!=u2)
        {
            delta=sqrt(*it)/(sqrt(u2)+sqrt(*it));
            ret1=as_2[*it]/(*it*4*M_PI);// f(a)/2
            //cout<<delta<<endl;
                for(int i=0; i<n; i++)
                {
                    z=(generatorFBN()*1.0/(RMAX*1.0))*(1.-2.*delta)+delta;
                    ret2+=totPgg(z); //jeden w pamieci z calki po Pqg
                    //z=(generatorFBN()*1.0/(RMAX*1.0));
                    //ret2+=Pqg(z);
                }
                tmp=*it++;
                ret1+=as_2[*it]/(*it*4*M_PI);//f(b)/2
                ret2=ret2/(n*1.)*(1.-2.*delta); //wartosc calki po dz
                ret+=ret1*abs(tmp-*it)*ret2; //elementarny trapez *dx * calka po dz
                ret1=0;
                ret2=0;
                //cout<<ret<<endl;
        }
        ret+=Integrand_as_2(kt2,u2);
        Tgs.push_back(exp(-ret));
        return Tgs.back();
    }
}

double Integrand_as_2(const double & kt2, const double & u2)
{
    if(kt2>u2)
        throw Blad("kt^2 wieksze od u^2", kt2, u2, kt2, u2);
    vector<double>::iterator test=find(u2s.begin(), u2s.end(), kt2);
    if((*test!=kt2))
            throw Blad("nie ma takiego kt^2", kt2, u2, kt2, u2);
    test=find(u2s.begin(), u2s.end(), u2);
    if((*test!=u2))
        throw Blad("nie ma alfy dla takiego u^2", kt2, u2, kt2, u2);
    else
    {
        test--;
        double ret=0;
        double ret1=0;
        double tmp=0;
        vector<double>::iterator it=find(u2s.begin(), u2s.end(), kt2);
        while(*it!=u2)
        {
            ret1=as_2[*it]/(*it*2*M_PI);// f(a)
            tmp=*it++;
            ret1+=as_2[*it]/(*it*2*M_PI);//f(b)
            ret+=ret1*abs(tmp-*it)/2.; //*dx/2
        }
//        for(int i=0; i<n; i++)
//        {
//            tmp=(generatorFBN()*1.0/(RMAX*1.0))*(u2-kt2)+kt2;
//            ret+=interpolacja(tmp)/(tmp*2*M_PI);
//
//            //cout<<(u2-kt2)*(1-2*delta)/(1.*n)*ret<<endl;
//        }
        cout<<"as_2 scalkowane od "<<kt2<<" do "<<u2<<" = "<<ret<<endl;
        return ret;
    }
}

double ftg(double *args, size_t dim, void *params)
{
    struct pars * fp = (struct pars *)params;
    //prams := {u2}
    //args := {pt2,z}
    return interpolacja(args[0])/(args[0]*2*M_PI)*(totPgg(args[1])*theta(sqrt(args[0])/(sqrt(args[0])+sqrt(fp->u2)),args[1])+Pqg(args[1]));
}

void draw_gluons()
{
    vector<double> kt2s={
    /*0.001225,*/
714285.715501964,
1428571.42977893,
2142857.14405589
/*,
2857142.85833286,
3571428.57260982,
4285714.28688679,
5000000.00116375,
5714285.71544071,
6428571.42971768,
7142857.14399464,
7857142.8582716,
8571428.57254857,
9285714.28682553,
10000000.0011025,
10714285.7153795,
11428571.4296564,
12142857.1439334,
12857142.8582104,
13571428.5724873,
14285714.2867643,
15000000.0010412,
15714285.7153182,
16428571.4295952,
17142857.1438721,
17857142.8581491,
18571428.5724261,
19285714.286703,
20000000.00098,
20714285.715257,
21428571.4295339,
22142857.1438109,
22857142.8580879,
23571428.5723648,
24285714.2866418,
25000000.0009188,
25714285.7151957,
26428571.4294727,
27142857.1437497,
27857142.8580266,
28571428.5723036,
29285714.2865805,
30000000.0008575,
30714285.7151345,
31428571.4294114,
32142857.1436884,
32857142.8579654,
33571428.5722423,
34285714.2865193,
35000000.0007963,
35714285.7150732,
36428571.4293502,
37142857.1436271,
37857142.8579041,
38571428.5721811,
39285714.286458,
40000000.000735,
40714285.715012,
41428571.4292889,
42142857.1435659,
42857142.8578429,
43571428.5721198,
44285714.2863968,
45000000.0006737,
45714285.7149507,
46428571.4292277,
47142857.1435046,
47857142.8577816,
48571428.5720585,
49285714.2863355,
50000000.0006125,
50714285.7148894,
51428571.4291664,
52142857.1434434,
52857142.8577203,
53571428.5719973,
54285714.2862742,
55000000.0005512,
55714285.7148282,
56428571.4291051,
57142857.1433821,
57857142.8576591,
58571428.571936,
59285714.286213,
60000000.0004899,
60714285.7147669,
61428571.4290439,
62142857.1433208,
62857142.8575978,
63571428.5718747,
64285714.2861517,
65000000.0004287,
65714285.7147056,
66428571.4289826,
67142857.1432596,
67857142.8575365,
68571428.5718135,
69285714.2860905,
70000000.0003674,
70714285.7146444,
71428571.4289214,
72142857.1431983,
72857142.8574753,
73571428.5717523,
74285714.2860292,
75000000.0003062,
75714285.7145832,
76428571.4288602,
77142857.1431371,
77857142.8574141,
78571428.5716911,
79285714.285968,
80000000.000245,
80714285.714522,
81428571.4287989,
82142857.1430759,
82857142.8573529,
83571428.5716299,
84285714.2859068,
85000000.0001838,
85714285.7144608,
86428571.4287377,
87142857.1430147,
87857142.8572917,
88571428.5715686,
89285714.2858456,
90000000.0001226,
90714285.7143995,
91428571.4286765,
92142857.1429535,
92857142.8572305,
93571428.5715074,
94285714.2857844,
95000000.0000614,
95714285.7143383,
96428571.4286153,
97142857.1428923,
97857142.8571692,
98571428.5714462,
99285714.2857232,
100000000 */
    };
    vector<double> xs={
    0.99,
0.9735000167,
0.9570000333/*,
0.94050005,
0.9240000667,
0.9075000833,
0.8910001,
0.8745001167,
0.8580001333,
0.84150015,
0.8250001667,
0.8085001833,
0.7920002,
0.7755002167,
0.7590002333,
0.74250025,
0.7260002667,
0.7095002833,
0.6930003,
0.6765003167,
0.6600003333,
0.64350035,
0.6270003667,
0.6105003833,
0.5940004,
0.5775004167,
0.5610004333,
0.54450045,
0.5280004667,
0.5115004833,
0.4950005,
0.4785005167,
0.4620005333,
0.44550055,
0.4290005667,
0.4125005833,
0.3960006,
0.3795006167,
0.3630006333,
0.34650065,
0.3300006667,
0.3135006833,
0.2970007,
0.2805007167,
0.2640007333,
0.24750075,
0.2310007667,
0.2145007833,
0.1980008,
0.1815008167,
0.1650008333,
0.14850085,
0.1320008667,
0.1155008833,
0.0990009,
0.0825009167,
0.0660009333,
0.04950095,
0.0330009667,
0.0165009833,
0.000001,*/
    };
    vector<double> mu2s={
    1.69,
833335.00925,
1666668.3285
/*,
2500001.64775,
3333334.967,
4166668.28625,
5000001.6055,
5833334.92475,
6666668.244,
7500001.56325,
8333334.8825,
9166668.20175,
10000001.521,
10833334.84025,
11666668.1595,
12500001.47875,
13333334.798,
14166668.11725,
15000001.4365,
15833334.75575,
16666668.075,
17500001.39425,
18333334.7135,
19166668.03275,
20000001.352,
20833334.67125,
21666667.9905,
22500001.30975,
23333334.629,
24166667.94825,
25000001.2675,
25833334.58675,
26666667.906,
27500001.22525,
28333334.5445,
29166667.86375,
30000001.183,
30833334.50225,
31666667.8215,
32500001.14075,
33333334.46,
34166667.77925,
35000001.0985,
35833334.41775,
36666667.737,
37500001.05625,
38333334.3755,
39166667.69475,
40000001.014,
40833334.33325,
41666667.6525,
42500000.97175,
43333334.291,
44166667.61025,
45000000.9295,
45833334.24875,
46666667.568,
47500000.88725,
48333334.2065,
49166667.52575,
50000000.845,
50833334.16425,
51666667.4835,
52500000.80275,
53333334.122,
54166667.44125,
55000000.7605,
55833334.07975,
56666667.399,
57500000.7182501,
58333334.0375001,
59166667.3567501,
60000000.6760001,
60833333.9952501,
61666667.3145001,
62500000.6337501,
63333333.9530001,
64166667.2722501,
65000000.5915001,
65833333.9107501,
66666667.2300001,
67500000.5492501,
68333333.8685001,
69166667.1877501,
70000000.5070001,
70833333.8262501,
71666667.1455001,
72500000.4647501,
73333333.7840001,
74166667.1032501,
75000000.4225001,
75833333.7417501,
76666667.0610001,
77500000.3802501,
78333333.6995001,
79166667.0187501,
80000000.3380001,
80833333.6572501,
81666666.9765001,
82500000.2957501,
83333333.6150001,
84166666.9342501,
85000000.2535001,
85833333.5727501,
86666666.8920001,
87500000.2112501,
88333333.5305001,
89166666.8497502,
90000000.1690002,
90833333.4882502,
91666666.8075002,
92500000.1267502,
93333333.4460002,
94166666.7652502,
95000000.0845002,
95833333.4037502,
96666666.7230002,
97500000.0422502,
98333333.3615002,
99166666.6807502,
100000000 */
    };
    fstream save;
    string NAZWA="fragment_gluon";

    save.open(NAZWA,ios::out);

    for(double & x : xs)
    {
    for(double & val : kt2s)
        {
        for(double & mu2:mu2s)
            {
            save<<x<<"\t"<<val<<"\t"<<mu2<<"\t"<<fa(x,val,mu2)<<endl;
            //save<<val<<"\t"<<val*val<<"\t"<<((val+h)*(val+h)-(val-h)*(val-h))/(2*h)<<"\t"<<2*val<<endl;
            cout<<x<<"\t"<<val<<"\t"<<mu2<<endl;
            }
        }
        cout<<NAZWA<<"\t"<<x<<endl;
    }

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
