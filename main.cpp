#include "sudakov_g.h"
#include "Blad.h"

int main()
{
    try
    {
        /*Sudakov_g Obj;
        Obj.*/warm_up();
        /*Obj.*/read_alphas();
        clock_t start = clock();
        cout<<setprecision(10)<<endl;
        //cout<<Obj.Tg2(322.72326025,9.966698770556024e7)<<endl; //1.6900000000000002
        //cout<<Obj.Tg2(1.6900000000000002, 322.72326025)<<endl;
        //cout<<Tg2(2.5634955924640005e6, 7.5692940628561e7)<<endl;
        //cout<<endl;
        cout<<Tq(322.72326025,9.966698770556024e7)<<endl; //1.6900000000000002
        cout<<Tq(1.6900000000000002, 322.72326025)<<endl;
        //cout<<Tg(2.5634955924640005e6, 7.5692940628561e7)<<endl;
        cout<<endl;
        //cout<<Tgt(322.72326025,9.966698770556024e7)<<endl; //1.6900000000000002
        //cout<<Tgt(1.6900000000000002, 322.72326025)<<endl;
        //cout<<Tgt(2.5634955924640005e6, 7.5692940628561e7)<<endl;
        test_delta(2.5634955924640005e6, 7.5692940628561e7);
//        cout<<TG(2.5634955924640005e6, 7.5692940628561e7).first<<endl;
//        cout<<TG(322.72326025,9.966698770556024e7).first<<endl; //1.6900000000000002
        draw_gluons();
        //cout<<TG(1.6900000000000002, 322.72326025).first<<endl;
        //cout<<TG(2.5634955924640005e6, 7.5692940628561e7).first<<endl;
        cout << "Czas wykonywania: " << (double)(clock()-start)/(CLOCKS_PER_SEC) << "sek" << endl;
        cool_down();
        //Obj.make_lattice();
        //cout<<Obj.calka_test(0,2,0,4)<<endl;
        //cout << "Czas wykonywania: " << (double)(clock()-start)/(CLOCKS_PER_SEC) << "sek" << endl;
    }
    catch(Blad B)
    {
        B.what();
        cout<<B.a1<<" "<<B.a2<<" "<<B.a3<<" "<<B.a4<<" "<<endl;
    }
return 0;
}
