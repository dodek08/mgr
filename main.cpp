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


vector<double> xs =
{
    -13.815511,
-12.434965,
-11.054419,
-9.673872,
-8.293326,
-6.912780,
-5.532234,
-4.151688,
-2.771142,
-1.390596,
-0.010050
};

vector<double> kt2s =
{
    -6.704814,
-4.192265,
-1.679715,
0.832834,
3.345384,
5.857933,
8.370483,
10.883032,
13.395582,
15.908131,
18.420681
};
double x1,val1,mu21;

    for(double & x : xs)
    {
        x1=exp(x);

    for(double & val : kt2s)
    {
        val1=exp(val);
       
            cout<<x1<<" "<<val1<<" "<<FL(x1,val1)+FT(x1,val1)<<endl;
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
