
#include "sudakov_g.h"
#include "sudakov_updf.h"
#include "Blad.h"


double h=0.1; //przy rozniczkowaniu numerycznym
 const LHAPDF::PDF* pdfs = LHAPDF::mkPDF("CT10nlo", 0); //wazne!
multimap<tuple<double,double>,double> TG;
multimap<tuple<double,double>,double> TQ;

vector<double> Mu;
vector<double> Kt2;

vector<double> Mured; //(Mureduced czyli X z unikalnymi wartosciami X
vector<double> Kt2red;

double XMAX, XMIN, KT2MAX, KT2MIN, KTRANGE;
double Xindexgt;


void read_Ts()
{
    fstream f;
    double tmp1,tmp2,tmp3;
    f.open("siatkaTG", ios::in | ios::out);
    if (!f.is_open()){ throw Blad("zly plik wejscia, nie istnieje lub zle wprowadzony");}

    int i=-1;
    while(!f.eof())
    {
        i++;
        f>>tmp1; //kt2
        f>>tmp2;    //mu2
        f>>tmp3;    //TG
        TG.insert(pair<tuple<double,double>,double>(make_tuple(tmp1,tmp2),tmp3));

        Kt2.push_back(tmp1);
        Mu.push_back(tmp2);
    }

    cout<<"to były mu2!"<<endl;
    cout<<endl;
    cout<<endl;
    cout<<endl;

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
    cout<<endl;
    cout<<endl;
    cout<<endl;
    cout<<Kt2red.size()<<endl;
    cout<<endl;
    cout<<endl;
    cout<<endl;

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
    cout<<Mured.size()<<endl;
    // f.close();
    // f.open("siatkaTQ", ios::in | ios::out);
    // if (!f.is_open()){ throw Blad("zly plik wejscia, nie istnieje lub zle wprowadzony");}
    // while(!f.eof())
    // {
    //     f>>tmp1;
    //     f>>tmp2;
    //     f>>tmp3;
    //     TQ.insert(pair<tuple<double,double>,double>(make_tuple(tmp1,tmp2),tmp3);
    // }
    // f.close();
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
    vector<double> kt2s={
        //0.001225,
714285.715501964,
1428571.42977893,
2142857.14405589,
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
// 10000000
    };
    vector<double> xs={
    
    0.99,
0.9735000167,
0.9570000333,
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
0.000001

    };
    vector<double> mu2s={
    1.69,
833335.0092,
1666668.3285,
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
100000000


        };
   
    string NAZWA="siatki/gluon_";    
    clock_t t;
    t=clock();
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
    cout<<find_index_gt(kt2s,kt2s.back())<<endl;
    cout<<"znalazlo sie"<<endl;
    cout<<find_index_gt(Kt2red,100000000)<<endl;
    cout<<"znalazlo sie"<<endl;
    for(double & val : kt2s)
    {
        for(double & mu2 : mu2s)
        {
            if (val==100000000)
            {
                save<<x<<"\t"<<val<<"\t"<<mu2<<"\t"<<fa(x,100000000,mu2)<<endl; //wstyd mi za to rozwiazanie
            }
            else
            save<<x<<"\t"<<val<<"\t"<<mu2<<"\t"<<fa(x,val,mu2)<<endl;
        // cout/*<<x*/<<"\t"<<val<<"\t"<<mu2<<"\t"<<Tg(val,mu2)<<endl;
        }
    }
    cout<<s<<endl;
    t=clock()-t;
    cout<<"time sigma x: "<<((double)t)/CLOCKS_PER_SEC<<endl;
    save.close();
    } 

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