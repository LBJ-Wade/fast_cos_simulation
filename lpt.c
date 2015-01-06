#include <math.h>
#include <assert.h>
#include <gsl/gsl_rng.h>
#include "msg.h"
#include "mem.h"
#include "config.h"
#include "cosmology.h"
#include "power.h"
#include "particle.h"
#include "fft.h"
#include "lpt.h"

static unsigned int* seedtable;
static double boxsize;

//static int ix0= Local_x_start
//static int nx -- Local_nx

static size_t nc;
static size_t local_nx;
static size_t local_ix0;

FFT* fft_psi[3];    // Zeldovichi displacement Psi_i
FFT* fft_psi_ij[6]; // derivative Psi_i,j= dPsi_i/dq_j
FFT* fft_psi2[3];   // 2nd order displacement Psi(2)
FFT* fft_div_psi2;  // divergence of Psi(2)


static void set_seedtable(const int nc, gsl_rng* random_generator,
			  unsigned int* const stable);
static void lpt_generate_psi_k(const unsigned long seed,
			       PowerSpectrum* const ps);
static void lpt_compute_psi2_k(void);

void lpt_init(const int nc_, const double boxsize_, Mem* mem)
{
  boxsize= boxsize_;
  nc= nc_;

  msg_printf(msg_debug, "lpt_init(nc= %d, boxsize= %.1lf)\n", nc, boxsize);
  
  //if(mem == 0)
  //  mem= mem_init("mem_lpt");

  if(mem != NULL)
    mem_use_from_zero(mem, 0);
  
  for(int i=0; i<3; i++)
    fft_psi[i]= fft_alloc("Psi_i", nc, mem, 0);

  for(int i=0; i<6; i++)
    fft_psi_ij[i]= fft_alloc("Psi_ij", nc, mem, 0);

  for(int i=0; i<3; i++)
    fft_psi2[i]= fft_psi_ij[i];

  fft_div_psi2= fft_psi_ij[3];
  
  seedtable = malloc(nc*nc*sizeof(unsigned int)); assert(seedtable);

  // checks
  local_nx= fft_psi[0]->local_nx;
  local_ix0= fft_psi[0]->local_ix0;
  
  for(int i=0; i<3; i++) {
    assert(fft_psi[i]->nc == nc);
    assert(fft_psi[i]->local_nx == local_nx);
    assert(fft_psi[i]->local_ix0 == local_ix0);
  }
  for(int i=0; i<6; i++) {
    assert(fft_psi_ij[i]->nc == nc);
    assert(fft_psi_ij[i]->local_nx == local_nx);
    assert(fft_psi_ij[i]->local_ix0 == local_ix0);
  }
}

void set_seedtable(const int nc, gsl_rng* random_generator,
		   unsigned int* const stable)
{
  // from N-GenIC
  assert(stable);
  

  for(int i=0; i<nc/2; i++) {
    for(int j=0; j<i; j++)
      stable[i*nc + j] = 0x7fffffff * gsl_rng_uniform(random_generator);

    for(int j=0; j<i+1; j++)
      stable[j*nc + i] = 0x7fffffff * gsl_rng_uniform(random_generator);

    for(int j=0; j<i; j++)
      stable[(nc - 1 - i)*nc + j] = 
	0x7fffffff * gsl_rng_uniform(random_generator);

    for(int j=0; j<i+1; j++)
      stable[(nc - 1 - j)*nc + i] = 
	0x7fffffff * gsl_rng_uniform(random_generator);

    for(int j=0; j<i; j++)
      stable[i*nc + (nc - 1 - j)] = 
	0x7fffffff * gsl_rng_uniform(random_generator);

    for(int j=0; j<i+1; j++)
      stable[j*nc + (nc - 1 - i)] = 
	0x7fffffff * gsl_rng_uniform(random_generator);

    for(int j=0; j<i; j++)
      stable[(nc - 1 - i)*nc + (nc - 1 - j)] = 
	0x7fffffff * gsl_rng_uniform(random_generator);

    for(int j=0; j<i+1; j++)
      stable[(nc - 1 - j)*nc + (nc - 1 - i)] = 
	0x7fffffff * gsl_rng_uniform(random_generator);
  }
}

void lpt_generate_psi_k(const unsigned long seed, PowerSpectrum* const ps)
{
  // Generates 1LPT (Zeldovich) displacements, Psi_k
  // from N-GenIC by Volker Springel
  msg_printf(msg_verbose, "Generating delta_k...\n");
  msg_printf(msg_info, "Random Seed = %lu\n", seed);

  for(int i=0; i<3; i++) assert(fft_psi[i]);
  
  complex_t* psi_k[]= {fft_psi[0]->fk, fft_psi[1]->fk, fft_psi[2]->fk};

  const size_t nckz= nc/2 + 1;
  const double dk= 2.0*M_PI/boxsize;
  const double knq= nc*M_PI/boxsize; // Nyquist frequency
  const double fac= pow(2*M_PI/boxsize, 1.5);
  const double fac_2pi3= 1.0/(8.0*M_PI*M_PI*M_PI);
  
  gsl_rng* random_generator = gsl_rng_alloc(gsl_rng_ranlxd1);
  gsl_rng_set(random_generator, seed);
  set_seedtable(nc, random_generator, seedtable);

  
  // clean the delta_k grid
  for(size_t ix=0; ix<local_nx; ix++)
   for(size_t iy=0; iy<nc; iy++)
    for(size_t iz=0; iz<nckz; iz++)
      for(int i=0; i<3; i++) {
	size_t index= (ix*nc + iy)*nckz + iz;
	psi_k[i][index][0] = 0;
	psi_k[i][index][1] = 0;
      }

  double kvec[3];
  for(size_t ix=0; ix<nc; ix++) {
    size_t iix = nc - ix;
    if(iix == nc)
      iix = 0;

    if(!((local_ix0 <= ix  && ix  < (local_ix0 + local_nx)) ||
	 (local_ix0 <= iix && iix < (local_ix0 + local_nx))))
      continue;
    
    for(size_t iy=0; iy<nc; iy++) {
      gsl_rng_set(random_generator, seedtable[ix*nc + iy]);
      
      for(size_t iz=0; iz<nc/2; iz++) {
	double phase= gsl_rng_uniform(random_generator)*2*M_PI;
	double ampl;
	do
	  ampl = gsl_rng_uniform(random_generator);
	while(ampl == 0.0);
	
	if(ix == nc/2 || iy == nc/2 || iz == nc/2)
	  continue;
	if(ix == 0 && iy == 0 && iz == 0)
	  continue;
	
	if(ix < nc/2)
	  kvec[0]= dk*ix;
	else
	  kvec[0]= -dk*(nc - ix);
	
	if(iy < nc/2)
	  kvec[1]= dk*iy;
	else
	  kvec[1]= -dk*(nc - iy);
	
	if(iz < nc/2)
	  kvec[2]= dk*iz;
	else
	  kvec[2]= -dk*(nc - iz);
	
	double kmag2 = kvec[0]*kvec[0] + kvec[1]*kvec[1] + kvec[2]*kvec[2];
	double kmag = sqrt(kmag2);
	
#ifdef SPHEREMODE
	// select a sphere in k-space
	if(kmag > knq)
	  continue;
#else
	if(fabs(kvec[0]) > knq)
	  continue;
	if(fabs(kvec[1]) > knq)
	  continue;
	if(fabs(kvec[2]) > knq)
	  continue;
#endif
	
	//double p_of_k = PowerSpec(kmag); // = 1/(2pi)^3*P(k)
	//double pk= fac_2pi3*power_spectrum(ps, kmag);	
	//p_of_k *= -log(ampl);

	double delta2= -log(ampl)*fac_2pi3*power_spectrum(ps, kmag);
	
	double delta_k_mag= fac*sqrt(delta2);
	// delta_k_mag -- |delta_k| extrapolated to a=1
	// Displacement is extrapolated to a=1
	
	if(iz > 0) {
	  if(local_ix0 <= ix && ix < (local_ix0 + local_nx)) {
	    size_t index= ((ix - local_ix0)*nc + iy)*nckz + iz;
	    for(int i=0; i<3; i++) {
	      psi_k[i][index][0]= -kvec[i]/kmag2*delta_k_mag*sin(phase);
	      psi_k[i][index][1]=  kvec[i]/kmag2*delta_k_mag*cos(phase);
	    }
	  }
	}
	else { // k=0 plane needs special treatment
	  if(ix == 0) {
	    if(iy >= nc/2)
	      continue;
	    else {
	      if(local_ix0 <= ix && ix < (local_ix0 + local_nx)) {
		size_t iiy= nc - iy; // note: j!=0 surely holds at this point
		size_t index= ((ix - local_ix0)*nc + iy)*nckz + iz;				size_t iindex= ((ix - local_ix0)*nc + iiy)*nckz + iz;
		
		for(int i=0; i<3; i++) {
		  psi_k[i][index][0]=  -kvec[i]/kmag2*delta_k_mag*sin(phase);
		  psi_k[i][index][1]=   kvec[i]/kmag2*delta_k_mag*cos(phase);
		  
		  psi_k[i][iindex][0]= -kvec[i]/kmag2*delta_k_mag*sin(phase);
		  psi_k[i][iindex][1]= -kvec[i]/kmag2*delta_k_mag*cos(phase);
		}
	      }
	    }
	  }
	  else { // here comes i!=0 : conjugate can be on other processor!
	    if(ix >= nc/2)
	      continue;
	    else {
	      iix = nc - ix;
	      if(iix == nc)
		iix = 0;
	      int iiy = nc - iy;
	      if(iiy == nc)
		iiy = 0;
	      
	      if(local_ix0 <= ix && ix < (local_ix0 + local_nx)) {
		size_t index= ((ix - local_ix0)*nc + iy)*nckz + iz;
		for(int i=0; i<3; i++) {
		  psi_k[i][index][0]= -kvec[i]/kmag2*delta_k_mag*sin(phase);
		  psi_k[i][index][1]=  kvec[i]/kmag2*delta_k_mag*cos(phase);
		}
	      }
	      
	      if(local_ix0 <= iix && iix < (local_ix0 + local_nx)) {
		size_t index= ((iix - local_ix0)*nc + iiy)*nckz + iz;
		for(int i=0; i<3; i++) {
		  psi_k[i][index][0]= -kvec[i]/kmag2*delta_k_mag*sin(phase);
		  psi_k[i][index][1]= -kvec[i]/kmag2*delta_k_mag*cos(phase);
		}
	      }
	    }
	  }
	}
      }
    }
  }

  gsl_rng_free(random_generator);  
}

void lpt_compute_psi2_k(void)
{
  // Compute 2nd order Psi(2) from 1st order Psi
  //   Prerequest Psi_k  in fft_psi[]->fk
  //   Result     Psi2_k in fft_psi2[]->fk (Fourier space)
  
  msg_printf(msg_verbose, "Computing 2LPT displacement fields...\n");

  const size_t nckz= nc/2 + 1;
  const double dk= 2.0*M_PI/boxsize;

  complex_t* psi_k[]= {fft_psi[0]->fk, fft_psi[1]->fk, fft_psi[2]->fk};

  //const double fac = pow(2*M_PI/boxsize, 1.5);
  double kvec[3];

  //
  // 2nd order LPT
  //
  complex_t* psi_ij_k[6];
  for(int i=0; i<6; i++) {
    assert(fft_psi_ij[i]);
    psi_ij_k[i]= fft_psi_ij[i]->fk;
  }

  // Take derivative dPsi_i/dq_j in Fourier space
  for(size_t ix=0; ix<local_nx; ix++) {
    for(size_t iy=0; iy<nc; iy++) {
      for(size_t iz=0; iz<nckz; iz++) {
	size_t index= (ix*nc + iy)*nckz + iz;
	if((ix + local_ix0) < nc/2)
	  kvec[0]= dk*(ix + local_ix0);
	else
	  kvec[0]= -dk*(nc - (ix + local_ix0));
	      
	if(iy < nc/2)
	  kvec[1]= dk*iy;
	else
	  kvec[1]= -dk*(nc - iy);
	      
	if(iz < nc/2)
	  kvec[2]= dk*iz;
	else
	  kvec[2]= -dk*(nc - iz);
	      
	// Derivatives of ZA displacements
	// dPsi_i/dq_j -> sqrt(-1) k_j Psi_i(k)
	psi_ij_k[0][index][0]= -psi_k[0][index][1]*kvec[0]; // Psi_1,1
	psi_ij_k[0][index][1]=  psi_k[0][index][0]*kvec[0];

	psi_ij_k[1][index][0]= -psi_k[0][index][1]*kvec[1]; // Psi_1,2
	psi_ij_k[1][index][1]=  psi_k[0][index][0]*kvec[1];

	psi_ij_k[2][index][0]= -psi_k[0][index][1]*kvec[2]; // Psi_1,3
	psi_ij_k[2][index][1]=  psi_k[0][index][0]*kvec[2];
	      
	psi_ij_k[3][index][0]= -psi_k[1][index][1]*kvec[1]; // Psi_2,2
	psi_ij_k[3][index][1]=  psi_k[1][index][0]*kvec[1];

	psi_ij_k[4][index][0]= -psi_k[1][index][1]*kvec[2]; // Psi_2,3
	psi_ij_k[4][index][1]=  psi_k[1][index][0]*kvec[2];

	psi_ij_k[5][index][0]= -psi_k[2][index][1]*kvec[2]; // Psi_3,3
	psi_ij_k[5][index][1]=  psi_k[2][index][0]*kvec[2];
      }
    }
  }

  // Second-order displacement Psi(2)
  // div.Psi(2) = Sum_{i<j} [ Psi_i,j Psi_i,j - Psi_i,i Psi_j,j ]
  // in realspace

  msg_printf(msg_verbose, "Fourier transforming displacement gradient...\n");
  for(int i=0; i<6; i++) 
    fft_execute_inverse(fft_psi_ij[i]);

  float_t* psi_ij[6]; for(int i=0; i<6; i++) psi_ij[i]= fft_psi_ij[i]->fx;
  float_t* const div_psi2= fft_div_psi2->fx;

  size_t nczr= 2*(nc/2 + 1);
  for(size_t ix=0; ix<local_nx; ix++) {
    for(size_t iy=0; iy<nc; iy++) {
      for(size_t iz=0; iz<nc; iz++) {
	size_t index= (ix*nc + iy)*nczr + iz;

	div_psi2[index]=
	    psi_ij[0][index]*(psi_ij[3][index] + psi_ij[5][index])
	  + psi_ij[3][index]*psi_ij[5][index]
          - psi_ij[1][index]*psi_ij[1][index]
          - psi_ij[2][index]*psi_ij[2][index]
	  - psi_ij[4][index]*psi_ij[4][index];
      }
    }
  }

  // Solve Poisson eq. for div.Psi(2) in Fourier space
  msg_printf(msg_verbose, "Fourier transforming second order source...\n");
  
  fft_execute_forward(fft_div_psi2);
  complex_t* div_psi2_k= fft_div_psi2->fk;
  complex_t* psi2_k[]= {fft_psi2[0]->fk, fft_psi2[1]->fk, fft_psi2[2]->fk};

  if(local_ix0 == 0) {
    for(int i=0; i<3; i++)
      psi2_k[i][0][0]= psi2_k[i][0][1] = 0.0;
    // avoid zero division kmag2 = 0
  }
	    
  for(size_t ix=0; ix<local_nx; ix++) {
    for(size_t iy=0; iy<nc; iy++) {
      int iz0= (ix + local_ix0 == 0) && (iy == 0); // skip kvec=(0,0,0)
      for(size_t iz=iz0; iz<nckz; iz++) {
	size_t index= (ix*nc + iy)*nckz + iz;
	if((ix + local_ix0) < nc/2)
	  kvec[0]=  dk*(ix + local_ix0);
	else
	  kvec[0]= -dk*(nc - (ix + local_ix0));
	
	if(iy < nc/2)
	  kvec[1]= dk*iy;
	else
	  kvec[1]= -dk*(nc - iy);
	
	if(iz < nc/2)
	  kvec[2] = dk*iz;
	else
	  kvec[2] = -dk*(nc - iz);
	
	double kmag2= kvec[0]*kvec[0] + kvec[1]*kvec[1] + kvec[2]*kvec[2];
	assert(kmag2 > 0); // !!! Heavy assert. Remove later.
	    
	// Psi(2)_k = div.Psi(2)_k * k / (sqrt(-1) k^2)
	for(int i=0; i<3; i++) {
	  psi2_k[i][index][0]=  div_psi2_k[index][1]*kvec[i]/kmag2;
	  psi2_k[i][index][1]= -div_psi2_k[index][0]*kvec[i]/kmag2;
	}
      }
    }
  }
}

void lpt_set_displacements(const unsigned long seed, PowerSpectrum* const ps,
			   const double a, Particles* particles)
{
  assert(particles);
  size_t np_local= local_nx*nc*nc;
  if(particles->np_allocated < np_local)
    msg_abort("Error: Not enough particles allocated to put initial particles\n"
	      "np_allocated= %lu < required %lu\n",
	      particles->np_allocated, np_local);
 
  lpt_generate_psi_k(seed, ps);
  lpt_compute_psi2_k();

  // Prerequest: psi_k in fft_psi[]->fk and psi2_k in fft_psi2[]->fk
  
  // Convert Psi_k Psi2_k to realspace
  msg_printf(msg_verbose, "Fourier transforming 2LPT displacements\n");
  for(int i=0; i<3; i++) {
    fft_execute_inverse(fft_psi[i]);
    fft_execute_inverse(fft_psi2[i]);
  }

  float_t* psi[]=  {fft_psi[0]->fx, fft_psi[1]->fx, fft_psi[2]->fx};
  float_t* psi2[]= {fft_psi2[0]->fx, fft_psi2[1]->fx, fft_psi2[2]->fx};
  

  msg_printf(msg_verbose, "Setting particle grid and displacements\n");

  const size_t nczr= 2*(nc/2 + 1);
  const float_t dx= boxsize/nc;
  Particle* p= particles->p;

  double nmesh3_inv= 1.0/pow((double)nc, 3.0);
  uint64_t id= (uint64_t) local_ix0*nc*nc + 1;

  const float_t D1= cosmology_D_growth(a);
  const float_t D2= cosmology_D2_growth(a, D1);
  
  float_t x[3];
  for(size_t ix=0; ix<local_nx; ix++) {
   x[0]= (local_ix0 + ix + 0.5f)*dx; /// !!! should we have this 0.5?
   for(size_t iy=0; iy<nc; iy++) {
    x[1]= (iy + 0.5f)*dx;
    for(int iz=0; iz<nc; iz++) {
     x[2]= (iz + 0.5f)*dx;

     size_t index= (ix*nc + iy)*nczr + iz;
     for(int i=0; i<3; i++) {
       float_t dis=  psi[i][index];
       float_t dis2= nmesh3_inv*psi2[i][index];
       // psi2 had two inverse Fourier transofroms, giving additional nmesh3
       
       p->x[i]= x[i] + D1*dis + D2*dis2;
       p->dx1[i]= dis;              // 1LPT extrapolated to a=1
       p->dx2[i]= dis2;             // 2LPT displacement
                                    // multiply by cosmology_D2_growth() for a
       p->v[i]= 0;                  // velocity in comoving 2LPT	  
     }
     p->id= id++;
     p++;
    }
   }
  }

  msg_printf(msg_verbose, "2LPT displacements calculated.\n");
  
  particles->np_local= np_local; 
  particles->a_x= a;
  particles->a_v= 0.0;
}

