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
//static const int fullT= 1; // velocity growth model

//
//float growthD(const float a);
//float growthD2(const float a);
//double Sphi(double ai, double af, double aRef);
double Sq(double ai, double af, double aRef);
float Qfactor(const float a);
// Leap frog time integration

void cola_kick(Particles* const particles, const float avel1)
{
  const float ai=  particles->a_v;  // t - 0.5*dt
  const float a=   particles->a_x;  // t
  const float af=  avel1;           // t + 0.5*dt

  Om= particles->omega_m;
  msg_printf(msg_info, "Kick %g -> %g\n", ai, avel1);

  const float Om143= pow(Om/(Om + (1 - Om)*a*a*a), 1.0/143.0);
  //const float dda= Sphi(AI, AF, A);
  //const float dda_new= Sphi_new(AI, AF, A);
  //printf("Sphi(%e %e %e) -> %e\n", AI, AF, A, dda);
  //const float growth1=growthD(A);
  const float kick_factor= (pow(af, nLPT) - pow(ai, nLPT))/
                           (nLPT*pow(a, nLPT)*sqrt(Om/a+(1.0-Om)*a*a));
  const float growth1= cosmology_D_growth(a);
  //const float growth2= cosmology_D2_growth(A, growth1);
  //fprintf(stderr, "%e %e %e\n", dda, kick_factor, fabs(dda-kick_factor)/dda);
	
  msg_printf(msg_debug, "growth factor %g\n", growth1);

  // Need to review/rederive
  const float q2=-3.0/7.0*pow(Om, -1.0/143.0)*1.5*Om*growth1*growth1*(1.0 + 7.0/3.0*Om143);
  const float q1=1.5*Om*growth1;

  //fprintf(stderr, "check %e %e\n", growth1*growth1*(1.0 + 7.0/3.0*Om143),
  //growth2);

  //fprintf(stderr, "new %e %e %e %e\n", A, q1, q2, dda);

  //printf("kick_factor %e %e %e\n", kick_factor, q1, q2);
  //printf("Om %e\n", Om);
  //abort();

  Particle* const p= particles->p;
  const int np= particles->np_local;
  float3* const f= particles->force;

  //FILE* fp= fopen("kick.txt", "w");
  
  // Kick using acceleration at a= A
  // Assume forces at a=A is in particles->force
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

    //fprintf(stderr, "%e %e %e\n", A, -1.5*Om*f[i][0],
    //P[i].dx1[0]*q1 + P[i].dx2[0]*q2);
    // check here, I don't understand why not -1.5*Om*(f[i][0] - P[i].dx1[0]
    // !!!
    
    p[i].v[0] += ax*kick_factor;
    p[i].v[1] += ay*kick_factor;
    p[i].v[2] += az*kick_factor;

    //printf("%e\n", p[i].v[0]);
  }
  //fclose(fp); abort();
  
  //velocity is now at a= avel1
  particles->a_v= avel1;

  //write_particles_txt("particles_kicked.txt", particles); abort();
}

void cola_drift(Particles* const particles, const float apos1)
{
  const float ai=  particles->a_x; // A:  t
  //const float av= particles->a_v; // AC: t + 0.5*dt
  const float af= apos1;          // AF: t + dt
  Om= particles->omega_m;
  
  Particle* const P= particles->p;
  const int np= particles->np_local;

  
  const float dyyy=Sq(ai, af, particles->a_v);

  //const float da1= growthD(AF) - growthD(A);    // change in D_{1lpt}
  //const float da2= growthD2(AF) - growthD2(A);  // change in D_{2lpt}

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

    //printf("%e\n", P[i].v[0]*dyyy + (P[i].dx1[0]*da1 + P[i].dx2[0]*da2));
  }
    
  particles->a_x= af;
}

/*

float growthDtemp(const float a){
    // Decided to use the analytic expression for LCDM. More transparent if I change this to numerical integration?
    float x=-Om/(Om - 1.0)/(a*a*a);
    
    
    float hyperP=0,hyperM=0;
    
    if (fabs(x-1.0) < 1.e-3) {
      hyperP= 0.859596768064608 - 0.1016599912520404*(-1.0 + x) + 0.025791094277821357*pow(-1.0 + x,2) - 0.008194025861121475*pow(-1.0 + x,3) + 0.0029076305993447644*pow(-1.0 + x,4) - 0.0011025426387159761*pow(-1.0 + x,5) + 0.00043707304964624546*pow(-1.0 + x,6) - 0.0001788889964687831*pow(-1.0 + x,7);
      hyperM= 1.1765206505266006 + 0.15846194123099624*(-1.0 + x) - 0.014200487494738975*pow(-1.0 + x,2) + 0.002801728034399257*pow(-1.0 + x,3) - 0.0007268267888593511*pow(-1.0 + x,4) + 0.00021801569226706922*pow(-1.0 + x,5) - 0.00007163321597397065*pow(-1.0 + x,6) +    0.000025063737576245116*pow(-1.0 + x,7);
    }
    else {
      if (x < 1.0) {
	hyperP=gsl_sf_hyperg_2F1(1.0/2.0,2.0/3.0,5.0/3.0,-x);
	hyperM=gsl_sf_hyperg_2F1(-1.0/2.0,2.0/3.0,5.0/3.0,-x);
      }
      x=1.0/x;
      if ((x < 1.0) && (x>1.0/30)) {
	
	hyperP=gsl_sf_hyperg_2F1(-1.0/6.0,0.5,5.0/6.0,-x);
	hyperP*=4*sqrt(x);
	hyperP+=-3.4494794123063873799*pow(x,2.0/3.0);
        
	hyperM=gsl_sf_hyperg_2F1(-7.0/6.0,-0.5,-1.0/6.0,-x);
	hyperM*=4.0/7.0/sqrt(x);
	hyperM+=pow(x,2.0/3.0)*(-1.4783483195598803057); //-(Gamma[-7/6]*Gamma[5/3])/(2*sqrt[Pi])
      }
      if (x<=1.0/30.0){
	hyperP=3.9999999999999996*sqrt(x) - 3.4494794123063865*pow(x,0.6666666666666666) + 0.3999999999999999*pow(x,1.5) -    0.13636363636363635*pow(x,2.5) + 0.07352941176470587*pow(x,3.5) - 0.04755434782608695*pow(x,4.5) +    0.033943965517241374*pow(x,5.5) - 0.02578125*pow(x,6.5) + 0.020436356707317072*pow(x,7.5) -    0.01671324384973404*pow(x,8.5) + 0.013997779702240564*pow(x,9.5) - 0.011945562847590041*pow(x,10.5) + 0.01035003662109375*pow(x,11.5) - 0.009080577904069926*pow(x,12.5);
	hyperM=0.5714285714285715/sqrt(x) + 2.000000000000001*sqrt(x) - 1.4783483195598794*pow(x,0.66666666666666666) +    0.10000000000000002*pow(x,1.5) - 0.022727272727272735*pow(x,2.5) + 0.009191176470588237*pow(x,3.5) -    0.004755434782608697*pow(x,4.5) + 0.002828663793103449*pow(x,5.5) - 0.0018415178571428578*pow(x,6.5) +    0.0012772722942073172*pow(x,7.5) - 0.0009285135472074472*pow(x,8.5) + 0.0006998889851120285*pow(x,9.5) -    0.0005429801294359111*pow(x,10.5) + 0.0004312515258789064*pow(x,11.5) - 0.00034925299631038194*pow(x,12.5);
      }
    }
    
    
    if (a > 0.2) 
      return sqrt(1.0 + (-1.0 + pow(a,-3))*Om)*(3.4494794123063873799*pow(-1.0 + 1.0/Om,0.666666666666666666666666666) + (hyperP*(4*pow(a,3)*(-1.0 + Om) - Om) - 7.0*pow(a,3)*hyperM*(-1.0 + Om))/(pow(a,5)*(-1.0+ Om) - pow(a,2)*Om));

    return (a*pow(1 - Om,1.5)*(1291467969*pow(a,12)*pow(-1 + Om,4) + 1956769650*pow(a,9)*pow(-1 + Om,3)*Om + 8000000000*pow(a,3)*(-1 + Om)*pow(Om,3) + 37490640625*pow(Om,4)))/(1.5625e10*pow(Om,5));    
}

float growthD(const float a) { // growth factor for LCDM
    return growthDtemp(a)/growthDtemp(1.0);
}


float Qfactor(const float a) { // Q\equiv a^3 H(a)/H0.
    return sqrt(Om/(a*a*a)+1.0-Om)*a*a*a;
}




float growthD2temp(const float a){
    float d= growthD(a);
    float omega=Om/(Om + (1.0 - Om)*a*a*a);
    return d*d*pow(omega, -1.0/143.);
}

float growthD2(const float a) {// Second order growth factor
  return growthD2temp(a)/growthD2temp(1.0); // **???
}


float growthD2v(const float a){ // explanation is in main()
    float d2= growthD2(a);
    float omega=Om/(Om + (1.0 - Om)*a*a*a);
    return Qfactor(a)*(d2/a)*2.0*pow(omega, 6.0/11.);
}

float decayD(float a){ // D_{-}, the decaying mode
    return sqrt(Om/(a*a*a)+1.0-Om);
}

double DprimeQ(double a,float nGrowth)
{ // returns Q*d(D_{+}^nGrowth*D_{-}^nDecay)/da, where Q=Qfactor(a)
  float nDecay=0.0;// not interested in decay modes in this code.
  float Nn=6.0*pow(1.0 - Om,1.5)/growthDtemp(1.0);
  return (pow(decayD(a),-1.0 + nDecay)*pow(growthD(a),-1.0 + nGrowth)*(nGrowth*Nn- (3.0*(nDecay + nGrowth)*Om*growthD(a))/(2.*a)));  
}


*/ 

//
// Functions for our modified time-stepping (used when StdDA=0):
//

//double gpQ(double a) { 
//  return pow(a, nLPT);
//}

double fun (double a, void * params) {
  //double f;
  /*
  if (fullT==1) f = gpQ(a)/Qfactor(a); 
  else f = 1.0/Qfactor(a);
  */

  return pow(a, nLPT)/(sqrt(Om/(a*a*a) + 1.0 - Om)*a*a*a);
  
  //return f;
}


/*     
      When StdDA=0, one needs to set fullT and nLPT.
         fullT=0 assumes time dependence for velocity = A + B a^nLPT, with A>>B a^nLPT. (A and B are irrelevant)
         fullT=1 assumes time dep. for velocity = B a^nLPT
         nLPT is a real number. Sane values lie in the range (-4,3.5). Cannot be 0, but of course can be -> 0 (say 0.001).
         See Section A.3 of TZE.
*/

double Sq(double ai, double af, double av) {
  gsl_integration_workspace * w 
    = gsl_integration_workspace_alloc (5000);
  
  double result, error;
  double alpha=0;
  
  gsl_function F;
  F.function = &fun;
  F.params = &alpha;
  
  gsl_integration_qag (&F, ai, af, 0, 1e-5, 5000,6,
		       w, &result, &error); 
  
  gsl_integration_workspace_free (w);
     
  //if (fullT==1)
  //  return result/gpQ(aRef);
  //return result;
  return result/pow(av, nLPT);
}

/*
double DERgpQ(double a) { // This must return d(gpQ)/da
  return nLPT*pow(a, nLPT-1);
}

double Sphi(double ai, double af, double aRef) {
  double result;
  result=(gpQ(af)-gpQ(ai))*aRef/Qfactor(aRef)/DERgpQ(aRef);
  
  return result;
}

*/
