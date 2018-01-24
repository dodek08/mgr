#include "sudakov_g.h"
#include "sudakov_updf.h"
#include "Blad.h"

const gsl_rng_type * T = gsl_rng_taus113; //generator do wyboru z gsl
gsl_rng * r=gsl_rng_alloc(T);
gsl_rng * r=gsl_rng_alloc(T);
gsl_monte_function U; //Tg
gsl_monte_vegas_state *s;

double F2(const double & x, const double & Q2)
{
	return F2_u(x,Q2);
}

double F2_u(const double & x, const double &Q2)
{

}

void warm_up()
{
  gsl_rng_env_setup ();
  gsl_rng_set(r, chrono::system_clock::now().time_since_epoch().count());
  U.f=&f2_u;
  U.dim=2;
}

double f2_u(double *args, size_t dim, void *params)
{
    struct pars * fp = (struct pars *)params;
    //prams := {u2}
    //args := {pt2,z}
    return interpolacja(args[0])/(args[0]*2*M_PI)*(args[1]*totPgg(args[1])*theta(sqrt(args[0])/(sqrt(args[0])+sqrt(fp->u2)),args[1])+5.0*Pqg(args[1]));

    //params := {x,Q2}
    //args :={kappa, kt, z}
    
}