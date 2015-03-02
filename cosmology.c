///
/// \file  cosmology.c
/// \brief Analytical functions in cosmology, e.g. linear growth rate
///
/// * Flat Lambda CDM assumed
///

#include <math.h>
#include <assert.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include "msg.h"
#include "cosmology.h"

static double omega_m0;
static double growth_normalisation;

static double growth_integrand(double a, void* param);
static double growth_unnormalised(const double a);

void cosmology_init(const double omega_m0_)
{
  omega_m0= omega_m0_;
  growth_normalisation= 1.0/growth_unnormalised(1.0); // D_growth=1 at a=1
}

void cosmology_check(void)
{
  // Check if this module is initilised
  if(growth_normalisation == 0.0)
    msg_abort("Error: cosmology module not initialised.\n");
  assert(growth_normalisation > 0.0);
}
    
double cosmology_D_growth(const double a)
{
  cosmology_check();
  // Linear growth factor D
  if(a == 0.0) return 0.0;
  
  return growth_normalisation*growth_unnormalised(a);
}

double cosmology_D2_growth(const double a, const double D)
{
  // 2nd-order growth factor D2
  if(a == 0.0) return 0.0;
  
  return -3.0/7.0*D*D*pow(cosmology_omega(a), -1.0/143.0);
}

double cosmology_f_growth_rate(const double a)
{
  cosmology_check();

  if(a == 0.0) return 1.0;
  
  // Linear growth rate f=dlnD/dlna
  const double d_un= growth_unnormalised(a);
  const double hf= cosmology_hubble_function(a);

  return 1.0/(d_un*a*a*hf*hf) - 1.5*omega_m0/(hf*hf*a*a*a);   
}  


void cosmology_growth(const double a,
		      double* const D_result, double* const f_result)
{
  cosmology_check();
  // Both linear growth factor D(a) and growth rate f=dlnD/dlna

  if(a == 0.0) {
    *D_result= 0.0;
    *f_result= 1.0;
    return;
  }
  
  const double d_un= growth_unnormalised(a);
  const double hf= cosmology_hubble_function(a);

  *D_result= growth_normalisation*d_un;
  *f_result= 1.0/(d_un*a*a*hf*hf) - 1.5*omega_m0/(hf*hf*a*a*a);   
}

double cosmology_hubble_function(const double a)
{
  // H/H0= sqrt(Omega_m0*a^-3 + Omega_Lambda)
  return sqrt(omega_m0/(a*a*a) + (1 - omega_m0));
}

double cosmology_omega(const double a)
{
  // Omega_m(a)
  return omega_m0/(omega_m0 + (1 - omega_m0)*(a*a*a));
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

