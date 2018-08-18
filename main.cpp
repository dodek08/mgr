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
        string siatka = argv[1];
        // int lp = atoi(argv[2]);
        // int l2 = atoi(argv[3]);
        // int l3 = atoi(argv[4]);
        double n_g = atof(argv[2]);
        double lambda_g = atof(argv[3]);
        double beta_g = atof(argv[4]);
        set_pdf_name_sudakov_g(siatka);
        set_pdf_name_sudakov_cs(siatka);
        set_pdf_name_sudakov_f2(siatka);
        set_pdf_name_sudakov_updf(siatka);
        clock_t start = clock();
        cout<<setprecision(10)<<endl;

Wyniki HERA("HERA1+2_NCep_920.dat");
cout<<HERA.nazwa<<endl;
// cout<<lp<<"\t"<</*l2<<*/endl;
// cout<<"n_g"<<"\t"<<"lambda_g"<<"\t"<<"beta_g"<<"\t"<<"X2"<<"\t"<<"av_my_err"<</*"\t"<<"f2_full"<<*/endl;
// cout<<"x"<<"\t"<<"Q2"<<"\t"<<"cs_calc"<<"\t"<<"cs_data"<<"\t"<<"err"<<endl;
cout<<"x"<<"\t"<<"Q2"<<"\t"<<"cs_data"<<"\t"<<"err"<<"\t"<<"cs_calc"<<"\t"<<"F2_q"<<"\t"<<"Fl"<<"\t"<<"Ft"<<endl;
// double alpha =7.297e-3;
// double xmp0 = 0.93827;
vector<double>::iterator x = HERA.x.begin();
vector<double>::iterator Q2 = HERA.Q2.begin();
vector<double>::iterator y = HERA.y.begin();
vector<double>::iterator sigma = HERA.sigma.begin();
vector<double>::iterator err = HERA.totnoproc.begin();
double av_my_err=0;
double X2 =0;
int licznik =1;

draw_gluons(21);

// set_n_g(-31);
// set_beta_g(0.5);
// set_lambda_g(0.5);
// set_n_g(-0.6);
// set_beta_g(4.5);
// set_lambda_g(-0.9);


// for(double beta_g = lp-59; beta_g>lp-60; beta_g-=0.1)
// {
	// for(double lambda_g = l2-1; lambda_g>(l2-2); lambda_g-=0.1)
	// for(double lambda_g = -0.4; lambda_g>-1.4; lambda_g-=0.1)
	// {
     // double n_g = lp*0.3-6;
     // double beta_g = l2*0.6-2;
		// for(double beta_g = 0; beta_g<2; beta_g+=0.5)
		// {
// double beta_g = lp*0.5;
//  for(double n_g = -.5; n_g>-0.65; n_g-=0.05)
// {
//  for(double lambda_g = -0.8; lambda_g>-0.95; lambda_g-=0.05)
// 	{
set_n_g(n_g);
set_lambda_g(lambda_g);
set_beta_g(beta_g);
// cout<<get_n_g()<<"\t"<<get_lambda_g()<<"\t"<<get_beta_g()<<endl;

while(x!=HERA.x.end())
{
	if(*Q2>=1.5 and *x<0.01)
	{
		// cout<<licznik<<"\t";
	// double cs = cs_tau_f2(*x++,*Q2++,*y++);
	// double fl = FT_g(*x,*Q2);
	// double ft = FL_g(*x,*Q2);
	// double f2 = F2_q(*x,*Q2);
	// double cs = (fl+ft+f2);//-(*y*(*y))/(1.+(1.-*y)*(1.-*y))*fl;
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
	X2 += ((cs-*sigma)*(cs-*sigma))/((*err*0.01)*(*sigma)*(*err*0.01)*(*sigma));
	// av_my_err += abs(cs/(*sigma)-1.);
    // cout<<(*x)<<"\t"<<(*Q2)<<"\t"<<*sigma<<"\t"<<*err<<"\t"<<cs<<"\t"<<f2<<"\t"<<fl<<"\t"<<ft<<endl;
	x++;
	y++;
	Q2++;
	sigma++;
	err++;
	licznik++;
	}
	else
	{
	x++;
	y++;
	Q2++;
	sigma++;
	err++;
	}
}

string NAZWA = "ch2.dat";
    fstream save;
    save.open(NAZWA,ios::out);
    if (!save.is_open()){ throw Blad("zly plik wejscia, nie istnieje lub zle wprowadzony");}

// X2=n_g*n_g+lambda_g*lambda_g+beta_g*beta_g+1.;

x = HERA.x.begin();
Q2 = HERA.Q2.begin();
y = HERA.y.begin();
sigma = HERA.sigma.begin();
err = HERA.totnoproc.begin();
// save<<get_n_g()<<"\t"<<get_lambda_g()<<"\t"<<get_beta_g()<<"\t";
// save<<X2<<"\t"<<av_my_err/(licznik*1.)<<endl;
save<<X2<<endl;
licznik = 0;
X2=0;
// av_my_err=0;

// }
// }
// vector<double>::iterator iter = v.begin();
//     for(double & x : xs)
//     {

//     for(double & Q2 : Q2s)
//     {       

//             cout<<x<<"\t"<<Q2<<"\t"<<FL_g(x,Q2)+FT_g(x,Q2)<<"\t"<<F2_q(x,Q2)<<"\t"<<*iter++<<endl;
//     }
//     }



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
