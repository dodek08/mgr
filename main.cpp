#include "sudakov_g.h"
#include "sudakov_updf.h"
#include "Blad.h"

int main()
{
    try
    {
        warm_up();
        clock_t start = clock();
        cout<<setprecision(8)<<endl;
        draw_gluons();
        draw_quarks_sudakov_factor();
        cout << "Czas wykonywania: " << (double)(clock()-start)/(CLOCKS_PER_SEC) << "sek" << endl;
        cool_down();
        print_time();

    }
    catch(Blad B)
    {
        B.what();
        cout<<B.a1<<" "<<B.a2<<" "<<B.a3<<" "<<B.a4<<" "<<endl;
    }
return 0;
}
