#include "sudakov_g.h"
#include "sudakov_updf.h"
#include "Blad.h"
#include "sudakov_f2.h"
#include "sudakov_cs.h"

LHAPDF::PDF* pdf3;// = LHAPDF::mkPDF("CT14nlo", 0);

void set_pdf_name_sudakov_cs(string siatka)
{
 pdf3 = LHAPDF::mkPDF(siatka, 0);
}


double cs_tau(const double & x, const double & Q2, const double & y)
{
	double fl = FT_g(x,Q2);
	double ft = FL_g(x,Q2);
	return (fl+ft)-(y*y)/(1.+(1.-y)*(1.-y))*fl;
}

double cs_tau_f2(const double & x, const double & Q2, const double & y)
{
	double fl = FT_g(x,Q2);
	double ft = FL_g(x,Q2);
	double f2 = F2_q(x,Q2);
	return (fl+ft+f2)-(y*y)/(1.+(1.-y)*(1.-y))*fl;
}

Wyniki::Wyniki(string name) {
	nazwa = name;
	fstream file;
    file.open(name, ios::in);
    if(file.good() == true){

        cout << "Uzyskano dostep do pliku \n";
        double tmp, tmp_x, tmp_q2, tmp_y, tmp_cs, tmp_totnoproc;
        while(!file.eof()) {
            file >> tmp_q2;
            file >> tmp_x;
            file >> tmp_y;
            file >> tmp_cs;
            // cout<<tmp_q2<<"\t"<<tmp_x<<"\t"<<tmp_y<<"\t"<<tmp_cs<<endl;
            for(int i=0; i<172; ++i) 
            {
            	file >> tmp;
                if(i==164)
            	tmp_totnoproc=tmp;
            }
            Q2.push_back(tmp_q2);
            x.push_back(tmp_x);
            y.push_back(tmp_y);
            sigma.push_back(tmp_cs);
            totnoproc.push_back(tmp_totnoproc);
        }
        file.close();

    } else cout << "Dostep do pliku zostal zabroniony! \n";
}


// double cs_tau(const double & x, const double & Q2, const double & s)
// {
// 	return 
// }