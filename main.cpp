#include "sudakov_g.h"
#include "sudakov_updf.h"
#include "Blad.h"

int main()
{
    try
    {
        warm_up();
        read_Ts();
        clock_t start = clock();
        cout<<setprecision(10)<<endl;
        draw_gluons();
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
