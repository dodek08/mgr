#include "sudakov_g.h"
#include "sudakov_updf.h"
#include "Blad.h"
#include "sudakov_f2.h"

const gsl_rng_type * T2 = gsl_rng_taus113; //generator do wyboru z gsl
gsl_rng *r2=gsl_rng_alloc(T2);
gsl_monte_function FL_u; //FL
gsl_monte_function FT_u;
gsl_monte_vegas_state *s2;
map<char,double> mq2;
map<char,double> cq2;
char quarks[] = {'u','d','s'};
const LHAPDF::PDF* pdf2 = LHAPDF::mkPDF("CT10nlo", 0);
size_t calls2 = 1000000;

double FL(const double & x, const double & Q2)
{
  double result=0., error=0.;   // result and error
  double xu[4]={kt2max(), kt2max(), 1., 2*M_PI};
  double xl[4]={kt2min(), kt2min(), 0, 0};
  double final_result=0;
  for(char q : {'u','d','s'})
  {
    result=0;
    error=0;
    // char q = quarks[i];
    // cout<<mq2[q]<<endl;
    struct pars2 pms={x, Q2, mq2[q]};
    FL_u.params=&pms;
    s2 = gsl_monte_vegas_alloc(4);
    //args :={kt2, Kt2, B, fi}
    gsl_monte_vegas_integrate(&FL_u, xl, xu, 4, calls2/10, r2, s2, &result, &error);
    // do
    // {
     result=0.;
     error=0.;
     gsl_monte_vegas_integrate(&FL_u, xl, xu, 4, calls2, r2, s2, &result, &error);
     // cout<<s2->chisq<<"\t"<<result<<"\t"<<error<<endl;
     // if(s2->chisq==0)
      // break;
    // }
    // while ((fabs (s2->chisq - 1.0) > 0.35) ); //more accurate
    gsl_monte_vegas_free(s2);
    final_result+=cq2[q]*result;
  }
  return final_result=Q2/(4*M_PI)*final_result;

}

double FT(const double & x, const double & Q2)
{
  double result=0., error=0.;   // result and error
  double xu[4]={kt2max(), kt2max(), 1., 2*M_PI};
  double xl[4]={kt2min(), kt2min(), 0, 0};
  double final_result=0;
  for(char q : {'u','d','s'})
  {
    result=0;
    error=0;
    // char q = quarks[i];
    // cout<<mq2[q]<<endl;
    struct pars2 pms={x, Q2, mq2[q]};
    FT_u.params=&pms;
    s2 = gsl_monte_vegas_alloc(4);
    //args :={kt2, Kt2, B, fi}
    gsl_monte_vegas_integrate(&FT_u, xl, xu, 4, calls2/10, r2, s2, &result, &error);
    // do
    // {
     result=0.;
     error=0.;
     gsl_monte_vegas_integrate(&FT_u, xl, xu, 4, calls2, r2, s2, &result, &error);
     // cout<<s2->chisq<<"\t"<<result<<"\t"<<error<<endl;
     // if(s2->chisq==0)
      // break;
    // }
    // while ((fabs (s2->chisq - 1.0) > 0.35) ); //more accurate
    gsl_monte_vegas_free(s2);
    final_result+=cq2[q]*result;
  }
  return final_result=Q2/(4*M_PI)*final_result;

}

void warm_up_f2()
{
  gsl_rng_env_setup ();
  gsl_rng_set(r2, chrono::system_clock::now().time_since_epoch().count());
  mq2['u'] = 0.0022*0.0022; //ALL UNITS IN GeV!
  cq2['u'] = 4./9.;
  mq2['d'] = 0.0047*0.0047; 
  cq2['d'] = 1./9.;
  mq2['s'] = 0.096*0.096;
  cq2['s'] = 1./9.;
  FL_u.f=&fl_u;
  FL_u.dim=4;
  FT_u.f=&ft_u;
  FT_u.dim=4;
}

double mu_2_v(const double & kt2, const double & Kt2, const double & mq2) //value of mu2 Kimber 3.15
{
  return kt2 + Kt2 + mq2;
}

double D1(const double & Kt2, const double & B, const double & Q2, const double & mq2) //Kimber 3.13
{
  return Kt2 + B*(1.-B)*Q2 + mq2;
}

double D2(const double & fi, const double & B, const double & Q2, const double & mq2) //Kimber 3.14
{
  return fi*fi + B*(1.-B)*Q2 + mq2;
}

double inv_z(const double & fi, const double & Kt2, const double & B, const double & Q2, const double & mq2) // 1/z, Kimber 3.10
{
  // return 1. + ((Kt2+mq2)*B+(fi*fi+mq2)*(1.-B))/(B*Q2*(1.-B));
  return 1. + (Kt2+mq2)/((1.-B)*Q2) + (fi*fi+mq2)/(B*Q2);
}

double heaviside(const double & x, const double & fi, const double & Kt2, const double & B, const double & Q2, const double & mq2)
{
  double ret = 1 - x*inv_z(fi,Kt2,B,Q2,mq2);
  if (ret>0)
  {
    ret = 1.;
  }
  else
  {
    ret = 0.;
  }

  return ret;
}

double fl_u(double *args, size_t dim, void *params)
{
    struct pars2 * fp = (struct pars2 *)params;
    //params := {x,Q2,mq2}
    //args :={kt2, Kt2, B, fi}
  double ret = 1. - fp->x*inv_z(args[3], args[1], args[2] , fp->Q2, fp->mq2);
  if (ret>0)
  {
    double d1 = D1(args[1], args[2], fp->Q2, fp->mq2);
    double d2 = D2(args[3], args[2], fp->Q2, fp->mq2);
    double mu2 = mu_2_v(args[1], args[2], fp->mq2);
    double invz = inv_z(args[3], args[1], args[2], fp->Q2, fp->mq2);
    ret = pdf2->alphasQ2(mu2)*fa(fp->x*invz,fp->x,mu2)*4.*fp->Q2*args[2]*args[2]*(1.-args[2])*(1.-args[2])*(1./d1-1./d2)*(1./d1-1./d2)/(args[0]*args[0]*2*M_PI);
  }
  else
  {
    ret = 0.;
  }
    return ret;
  }

double ft_u(double *args, size_t dim, void *params)
{
    struct pars2 * fp = (struct pars2 *)params;
    //params := {x,Q2,mq2}
    //args :={kt2, Kt2, B, fi}
  double ret = 1. - fp->x*inv_z(args[3], args[1], args[2] , fp->Q2, fp->mq2);
  if (ret>0)
  {
    double d1 = D1(args[1], args[2], fp->Q2, fp->mq2);
    double d2 = D2(args[3], args[2], fp->Q2, fp->mq2);
    double mu2 = mu_2_v(args[1], args[2], fp->mq2);
    double invz = inv_z(args[3], args[1], args[2], fp->Q2, fp->mq2);
    ret = pdf2->alphasQ2(mu2)*fa(fp->x*invz,fp->x,mu2)*((args[2]*args[2]+(1.-args[2])*(1.-args[2]))*(args[1]/d1-args[3]/d2)*(args[1]/d1-args[3]/d2)+fp->mq2*(1./d1-1./d2)*(1./d1-1./d2))/(args[0]*args[0]*2*M_PI);
  }
  else
  {
    ret = 0.;
  }
    return ret;
}