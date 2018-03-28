#include "sudakov_g.h"
#include "sudakov_updf.h"
#include "sudakov_f2.h"
#include "Blad.h"

int main()
{
    try
    {
        warm_up();
        warm_up_f2();
        read_Ts();
        clock_t start = clock();
        cout<<setprecision(10)<<endl;


vector<double> xs ={
    0.001,
    0.0016,
    0.0025,
    0.004,
    0.0063,
    0.01
};

vector<double> Q2s ={
    0.85,
    2,
    4.5,
    8.5,
    10,
    12,
    15,
    28,
    35,
    45
};

vector<double> v = {
    0.260152, 
0.378166, 
0.514264, 
0.6337, 
0.668525, 
0.705416, 
0.751432, 
0.886809, 
0.936697, 
0.989183, 
0.241045, 
0.347602, 
0.464871, 
0.571415, 
0.601343, 
0.634476, 
0.674635, 
0.793811, 
0.834733, 
0.888612, 
0.224397,
0.318606, 
0.424675, 
0.516548, 
0.543148, 
0.573198, 
0.610129, 
0.715383,
0.753623, 
0.799966, 
0.20682, 
0.291718, 
0.3831, 
0.466582, 
0.487928,
0.5128, 
0.546862, 
0.638276, 
0.669872, 
0.709097, 
0.191236, 
0.266765, 
0.347036, 
0.418379, 
0.440429, 
0.46199, 
0.489463, 
0.573636, 
0.60204,
0.636123, 
0.175609, 
0.24243, 
0.314371, 
0.378344, 
0.394634, 
0.412573, 
0.437676, 
0.510581, 
0.537518, 
0.562973
};
vector<double>::iterator iter = v.begin();
    for(double & x : xs)
    {

    for(double & Q2 : Q2s)
    {       
            // cout<<x<<" "<<Q2<<" "<<fgk(x,Q2)<<endl;
            cout<<x<<"\t"<<Q2<<"\t"<<FL(x,Q2)+FT(x,Q2)<<"\t"<<*iter++<<endl;
    }
    }




        cout << "Czas wykonywania: " << (double)(clock()-start)/(CLOCKS_PER_SEC) << "sek" << endl;
        cool_down();
    }
    catch(Blad B)
    {
        B.what();
        cout<<B.a1<<" "<<B.a2<<" "<<B.a3<<" "<<B.a4<<" "<<endl;
    }
return 0;
}
