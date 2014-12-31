//
// Analytical functions in cosmology, e.g., linear growth rate
//  * Flat Lambda CDM assumed
//

#include <math.h>
#include <assert.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include "cosmology.h"

static double omega_m0;
static double growth_fac;

static double growth_integrand(double a, void* param);
static double growth_unnormalised(const double a);

void cosmology_init(const double omega_m0_)
{
  omega_m0= omega_m0_;
  growth_fac= 1.0/growth_unnormalised(1.0); // D_growth=1 at a=1
}

double cosmology_D_growth(const double a)
{
  // Linear growth factor D
  assert(growth_fac > 0.0); // initialised?
  return growth_fac*growth_unnormalised(a);
}

double cosmology_f_growth_rate(const double a)
{
  // Linear growth rate f=dlnD/dlna
  const double d_un= growth_unnormalised(a);
  const double hf= cosmology_hubble_function(a);

  return 1.0/(d_un*a*a*hf*hf) - 1.5*omega_m0/(hf*hf*a*a*a);   
}  


void cosmology_growth(const double a,
		      double* const D_result, double* const f_result)
{
  // Both linear growth factor D(a) and growth rate f=dlnD/dlna
  
  const double d_un= growth_unnormalised(a);
  const double hf= cosmology_hubble_function(a);

  *D_result= growth_fac*d_un;
  *f_result= 1.0/(d_un*a*a*hf*hf) - 1.5*omega_m0/(hf*hf*a*a*a);   
}

double cosmology_hubble_function(const double a)
{
  // H/H0= sqrt(Omega_m0*a^-3 + Omega_Lambda)
  return sqrt(omega_m0/(a*a*a) + (1 - omega_m0));
}

double growth_integrand(double a, void* param)
{
  // sqrt[(a*H/H0)^3]
  //return pow(a/(omega_m0 + (1 - Omega)*a*a*a), 1.5);
  const double aHinv= 1.0/sqrt(omega_m0/a + (1 - omega_m0)*(a*a));
  return aHinv*aHinv*aHinv;
}

double growth_unnormalised(const double a)
{
  // D(a) \propto \int_0^a (a H(a)/H0)^-3 da
  const size_t worksize= 1000;

  gsl_integration_workspace *workspace=
    gsl_integration_workspace_alloc(worksize);

  gsl_function F;
  F.function = &growth_integrand;

  double result, abserr;
  gsl_integration_qag(&F, 0, a, 0, 0.5e-8, worksize, GSL_INTEG_GAUSS41,
		      workspace, &result, &abserr);

  gsl_integration_workspace_free(workspace);

  return cosmology_hubble_function(a) * result;
}

