#include "sudakov_g.h"
#include "sudakov_updf.h"
#include "sudakov_f2.h"
#include "sudakov_cs.h"
#include "Blad.h"

int main(int argc, char* argv[])
{
    try
    {
        warm_up();
        warm_up_f2();
        read_Ts();
        // read_ccfm();
        // read_pb();
        read_alphas();
        string siatka = argv[1];
        // int lp = atoi(argv[2]);
        // int l2 = atoi(argv[3]);
        // int l3 = atoi(argv[4]);
        // double n_g = atof(argv[2]);
        // double lambda_g = atof(argv[3]);
        // double beta_g = atof(argv[4]);
        set_pdf_name_sudakov_g(siatka);
        set_pdf_name_sudakov_cs(siatka);
        set_pdf_name_sudakov_f2(siatka);
        set_pdf_name_sudakov_updf(siatka);
        clock_t start = clock();



// double fl = 0;
// 	double ft = 0;
// 	double f2=0;

// std::vector<double> xs = {0.01, 0.001};
// std::vector<double> Q2s = {1,10,50};

// draw_gluons(21);


////Wyliczanie FX
// std::vector<double> xs = {0.001};
// std::vector<double> Q2s = {50};
// std::vector<double> calls = {10000, 100000, 1000000};
double call = 100000;
std::vector<double> factor = {10,16,25,40,63};
std::vector<double> Q2factor = {2.5,16,63};
cout<<setprecision(6);

	string NAZWA = siatka+".dat";
    fstream save;
    save.open(NAZWA,ios::out);
    if (!save.is_open()){ throw Blad("zly plik wejscia, nie istnieje lub zle wprowadzony");}
    save<<setprecision(6)<<endl;

    save<<"x"<<"\t"<<"Q2"<<"\t"<<"F2"<<"\t"<<"F2_q"<<"\t"<<"FL"<<"\t"<<"FT"<<endl;
    for (auto const& f : factor)
    {
    	/* code */
    
    for(int i = 2; i<7; i++)
    {
    	double x = f*pow(10,-i);
    	for(int j = 1; j<4; j++)
    	{
    		// for (auto const& fac : factor)
    		// {
    		// double Q2 = fac*pow(10,j-2);
    		for (auto const& fac : Q2factor)
    		{
    		double Q2 = fac;
    		double fl = FL_g(x,Q2);
		 	double ft = FT_g(x,Q2);
			double f2 = F2_q(x,Q2);
    		save<<x<<"\t"<<Q2<<"\t"<<fl+ft+f2<<"\t"<<f2<<"\t"<<fl<<"\t"<<ft<<endl;
			}    	
    	}
    }
}
//// wyliczanie rozkładów partonów dla wykresów
 //    save<<"x"<<"\t"<<"kt2"<<"\t"<<"Q2"<<"\t"<<"gluon"<<"\t"<<"d"<<"\t"<<"di"<<"\t"<<"u"<<"\t"<<"ui"<<"\t"<<"s"<<"\t"<<"si"<<"\t"<<"c"<<"\t"<<"ci"<<endl;
 //    for(auto const& x : xs)
 //    {
 //        for(auto const& Q2: Q2s )
 //        {
	// 	for(int i = 0; i<500; i ++)
	// 		{
	// 			double kt2 = i*0.2;
	// 			save<<x<<"\t"<<kt2<<"\t"<<Q2<<"\t"<<fq_pid(x,kt2,Q2,21)<<"\t"<<fq_pid(x,kt2,Q2,1)<<"\t"<<fq_pid(x,kt2,Q2,-1)<<"\t"<<fq_pid(x,kt2,Q2,2)<<"\t"<<fq_pid(x,kt2,Q2,-2)<<"\t"<<fq_pid(x,kt2,Q2,3)<<"\t"<<fq_pid(x,kt2,Q2,-3)<<"\t"<<fq_pid(x,kt2,Q2,4)<<"\t"<<fq_pid(x,kt2,Q2,-4)<<endl;
	// 		}
	// 	}
	// }
	// save.close();


//// Kompromis między czasem wykonywania a dokładnoscia
// cout<<"x"<<"\t"<<"Q2"<<"\t"<<"F2"<</*"\t"<<"F2_q"<<*/"\t"<<"Fl"<<"\t"<<"Ft"<<"\t"<<"calls"<<"\t"<<"t"<<endl;

//  for (auto const& call : calls)
//  {   
//     // for (int i = 0; i < 100; ++i)
//     // {
//     for(auto const& x : xs)
//     {
//         for(auto const& Q2: Q2s )
//         {
//         clock_t start = clock();
// 		fl = FL_g(x,Q2,call);
// 		ft = FT_g(x,Q2,call);
//         // f2= F2_q(x,Q2);
// 		cout<<x<<"\t"<<Q2<<"\t"<<fl+ft/*+f2<<"\t"<<f2*/<<"\t"<<fl<<"\t"<<ft<<"\t"<<call<<"\t"<<(double)(clock()-start)/(CLOCKS_PER_SEC)<<endl;
//     }
//     }
//  // }
// }
////Walidacja
    // for (int i = 0; i < 100; ++i)
    // {
 //    for(auto const& x : xs)
 //    {
 //        for(auto const& Q2: Q2s )
 //        {
 //            cout<<x<<"\t"<<Q2<<"\t";
 //            for (auto const& call : calls)
 //        { 
 //        clock_t start = clock();
 //        fl = FL_g(x,Q2,call);
 //        ft = FT_g(x,Q2,call);
 //        // f2= F2_q(x,Q2);
 //        cout<<fl+ft/*+f2<<"\t"<<f2*/<<"\t";//<<fl<<"\t"<<ft<<"\t"<<call<<"\t"<<(double)(clock()-start)/(CLOCKS_PER_SEC)<<endl;
 //    }
 //    cout<<endl;
 //    }
 // // }
 //    }


//// Fit i stare dzieje
// // string NAZWA = siatka+to_string(n_g)+to_string(lambda_g)+to_string(beta_g);
// //     fstream save;
// //     save.open(NAZWA,ios::out);
// //     if (!save.is_open()){ throw Blad("zly plik wejscia, nie istnieje lub zle wprowadzony");}

// Wyniki HERA("HERA1+2_NCep_920.dat");
// cout<<HERA.nazwa<<endl;

// // // draw_gluons(21);
// // // cout<<lp<<"\t"<</*l2<<*/endl;
// // // cout<<"n_g"<<"\t"<<"lambda_g"<<"\t"<<"beta_g"<<"\t"<<"X2"<<"\t"<<"av_my_err"<</*"\t"<<"f2_full"<<*/endl;
// // // cout<<"x"<<"\t"<<"Q2"<<"\t"<<"cs_calc"<<"\t"<<"cs_data"<<"\t"<<"err"<<endl;
// // cout<<"x"<<"\t"<<"Q2"<<"\t"<<"cs_data"<<"\t"<<"err"<<"\t"<<"cs_calc"<<"\t"<<"F2_q"<<"\t"<<"Fl"<<"\t"<<"Ft"<<endl;
// // double alpha =7.297e-3;
// // double xmp0 = 0.93827;
// vector<double>::iterator x = HERA.x.begin();
// vector<double>::iterator Q2 = HERA.Q2.begin();
// vector<double>::iterator y = HERA.y.begin();
// vector<double>::iterator sigma = HERA.sigma.begin();
// vector<double>::iterator err = HERA.totnoproc.begin();
// // double av_my_err=0;
// double X2 =0;
// // int licznik =1;


// set_n_g(n_g);
// set_lambda_g(lambda_g);
// set_beta_g(beta_g);
// // cout<<get_n_g()<<"\t"<<get_lambda_g()<<"\t"<<get_beta_g()<<endl;



// // for(double n_tmp = 0; n_tmp<0.1; n_tmp+=0.02)
// // {
// // 	set_n_g(n_g+n_tmp);
// // for(double lambda_tmp = 0; lambda_tmp<1; lambda_tmp+=0.2)
// // {set_lambda_g(lambda_g+lambda_tmp);
// 	// for (double beta_tmp = 0; beta_tmp<1; beta_tmp+=0.2)
// 	// {
// 	// 	set_beta_g(beta_g+beta_tmp);


// while(x!=HERA.x.end())
// {
	// if(*Q2>=2 and *x<0.01)
	// {
		// double f2 = F2_q(*x,*Q2);
		// double fl = FL_g(*x,*Q2);
	// double cs = cs_tau_f2(*x++,*Q2++,*y++);
	// double cs = (fl+ft+f2);//-(*y*(*y))/(1.+(1.-*y)*(1.-*y))*fl;
	// double cs = f2-(*y*(*y))/(1.+(1.-*y*(*y)))*fl;
	// double f2_full = fl+ft+f2;
	// double f2_data = *Q2/(4.*M_PI*M_PI*0.2)*(*sigma);
	// cout<<get_n_g()<<"\t"<<get_lambda_g()<<"\t"<<get_beta_g()<<"\t";
	// cout<<(*x)<<"\t"<<(*Q2)<<"\t"<<cs<<"\t"<<(*sigma)<<"\t"<<(*err)<<endl;
	// double cs =1.0;// xD
	// cout<<(*x)<<"\t"<<(*Q2)<<"\t"<<(*Q2)*(*Q2)<<"\t"<<ccfm_s((*x),(*Q2),(*Q2)*(*Q2))<<"\t"<<fa((*x),(*Q2),(*Q2)*(*Q2))<<endl;
	// double fac = pow((*Q2),2)*(1-(*x))/
// (4*pow(M_PI,2)*alpha*((*Q2)+pow(2*(*x)*xmp0,2)));
	// double fac = (*Q2)/(4*M_PI*M_PI*alpha);
 //            double units = 10.0 * pow(197.3271,2);
 //            fac = fac /(units * 1.e-3);
 //            double F2_data = fac*(*sigma);
 //            double err_data = fac*(*sigma)*(*err)*0.01;
	// cout<<get_n_g()<<"\t"<<get_lambda_g()<<"\t"<<get_beta_g()<<endl;
	// X2 += ((cs-*sigma)*(cs-*sigma))/((*err*0.01)*(*sigma)*(*err*0.01)*(*sigma));
	// av_my_err += abs(cs/(*sigma)-1.);
//     cout<<setprecision(6)<<(*err)<<"\t"<<(*x)<<"\t"<<(*Q2)<<"\t"<<*y<<"\t"<<*sigma<<"\t"<<cs<<"\t"<<f2<<"\t"<<fl<</*"\t"<<ft<<*/endl;
// 	x++;
// 	y++;
// 	Q2++;
// 	sigma++;
// 	err++;
// 	// licznik++;
// 	}
// 	else
// 	{
// 	x++;
// 	y++;
// 	Q2++;
// 	sigma++;
// 	err++;
// 	}
// }

// string NAZWA = "ch2.dat";
// // X2=n_g*n_g+lambda_g*lambda_g+beta_g*beta_g+1.;


    // fstream save;
    // save.open(NAZWA,ios::out);
    // if (!save.is_open()){ throw Blad("zly plik wejscia, nie istnieje lub zle wprowadzony");}

// x = HERA.x.begin();
// Q2 = HERA.Q2.begin();
// y = HERA.y.begin();
// sigma = HERA.sigma.begin();
// err = HERA.totnoproc.begin();
// // save<<get_n_g()<<"\t"<<get_lambda_g()<<"\t"<<get_beta_g()<<"\t";
// // save<<X2<<"\t"<<av_my_err/(licznik*1.)<<endl;
// save<<X2<<endl;
// licznik = 0;
// X2=0;

// // }
// // }
// // }
    save.close();

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
