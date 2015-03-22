#include <math.h>
#include <assert.h>
#include <mpi.h>

#include <gsl/gsl_integration.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_sf_hyperg.h> 
#include <gsl/gsl_errno.h>

#include "particle.h"
#include "msg.h"
#include "cola.h"
#include "cosmology.h"
#include "write.h"

static float Om= -1.0f;
static const float nLPT= -2.5f;

double Sq(double ai, double af, double aRef);
float Qfactor(const float a);

void cola_kick(Particles* const particles, const float avel1)
{
  const float ai=  particles->a_v;  // t - 0.5*dt
  const float a=   particles->a_x;  // t
  const float af=  avel1;           // t + 0.5*dt

  Om= particles->omega_m;
  msg_printf(msg_info, "Kick %g -> %g\n", ai, avel1);

  const float Om143= pow(Om/(Om + (1 - Om)*a*a*a), 1.0/143.0);
  const float kick_factor= (pow(af, nLPT) - pow(ai, nLPT))/
                           (nLPT*pow(a, nLPT)*sqrt(Om/a+(1.0-Om)*a*a));
  const float growth1= cosmology_D_growth(a);
	
  msg_printf(msg_debug, "growth factor %g\n", growth1);

  // Need to review/rederive !!!
  const float q2=-3.0/7.0*pow(Om, -1.0/143.0)*1.5*Om*growth1*growth1*(1.0 + 7.0/3.0*Om143);
  const float q1=1.5*Om*growth1;

  Particle* const p= particles->p;
  const int np= particles->np_local;
  float3* const f= particles->force;

  // Kick using acceleration at scale factor a
  // Assume forces at a is in particles->force
#ifdef _OPENMP
  #pragma omp parallel for default(shared)
#endif
  for(int i=0; i<np; i++) {
    float ax= -1.5*Om*f[i][0] - (p[i].dx1[0]*q1 + p[i].dx2[0]*q2);
    float ay= -1.5*Om*f[i][1] - (p[i].dx1[1]*q1 + p[i].dx2[1]*q2);
    float az= -1.5*Om*f[i][2] - (p[i].dx1[2]*q1 + p[i].dx2[2]*q2);

    /*
    fprintf(fp, "%e %e %e %e %e %e %e\n",
	    p[i].v[0],
	    -1.5*Om*f[i][0],
	    p[i].dx1[0],
	    p[i].dx2[0],
	    p[i].dx1[0]*q1 + p[i].dx2[0]*q2,
	    ax,
	    ax*kick_factor);
    */

    p[i].v[0] += ax*kick_factor;
    p[i].v[1] += ay*kick_factor;
    p[i].v[2] += az*kick_factor;

  }
  //fclose(fp); abort();
  
  //velocity is now at a= avel1
  particles->a_v= avel1;
}

void cola_drift(Particles* const particles, const float apos1)
{
  const float ai= particles->a_x;
  const float af= apos1;
  Om= particles->omega_m;
  
  Particle* const P= particles->p;
  const int np= particles->np_local;

  
  const float dyyy=Sq(ai, af, particles->a_v);

  const double growth_i= cosmology_D_growth(ai);
  const double growth_f= cosmology_D_growth(af);
  const float_t da1= growth_f - growth_i;

  const float_t da2= cosmology_D2_growth(af, growth_f) -
                     cosmology_D2_growth(ai, growth_i);

  msg_printf(msg_info, "Drift %g -> %g\n", ai, af);
    
  // Drift
#ifdef _OPENMP
  #pragma omp parallel for default(shared)
#endif
  for(int i=0; i<np; i++) {
    P[i].x[0] += P[i].v[0]*dyyy + 
                 (P[i].dx1[0]*da1 + P[i].dx2[0]*da2);
    P[i].x[1] += P[i].v[1]*dyyy +
                 (P[i].dx1[1]*da1 + P[i].dx2[1]*da2);
    P[i].x[2] += P[i].v[2]*dyyy + 
                 (P[i].dx1[2]*da1 + P[i].dx2[2]*da2);
  }
    
  particles->a_x= af;
}

double fun (double a, void * params) {
  return pow(a, nLPT)/(sqrt(Om/(a*a*a) + 1.0 - Om)*a*a*a);
}

double Sq(double ai, double af, double av) {
  gsl_integration_workspace * w 
    = gsl_integration_workspace_alloc (5000);
  
  double result, error;
  
  gsl_function F;
  F.function = &fun;
  F.params = 0;
  
  gsl_integration_qag (&F, ai, af, 0, 1e-5, 5000,6,
		       w, &result, &error); 
  
  gsl_integration_workspace_free (w);
     
  return result/pow(av, nLPT);
}

