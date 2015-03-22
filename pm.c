#include <math.h>
#include <assert.h>
#include <gsl/gsl_rng.h>
#include "msg.h"
#include "mem.h"
#include "config.h"
#include "cosmology.h"
#include "comm.h"
#include "particle.h"
#include "fft.h"

static int pm_factor;
static size_t nc, nzpad;
static float_t boxsize;

static FFT* fft_pm;
static complex_t* delta_k;

static inline void grid_assign(float_t * const d, 
	    const size_t ix, const size_t iy, const size_t iz, const float_t f)
{
#ifdef _OPENMP
  #pragma omp atomic
#endif
  d[(ix*nc + iy)*nzpad + iz] += f;
}

static inline float grid_val(float_t const * const d,
			const size_t ix, const size_t iy, const size_t iz)
{
  return d[(ix*nc + iy)*nzpad + iz];
}


static size_t send_buffer_positions(Particles* const particles);
static void pm_assign_cic_density(Particles* particles, size_t np);
static void check_total_density(float_t const * const density);
static void compute_delta_k(void);
static void compute_force_mesh(const int k);
static void force_at_particle_locations(
		 Particles* const particles, const int np, const int axis);
static void add_buffer_forces(Particles* const particles, const size_t np);

//
// Public functions
//
void pm_init(const int nc_pm, const int pm_factor_,
	     Mem* const mem_pm, Mem* const mem_density,
	     const float_t boxsize_)
{
  msg_printf(msg_verbose, "PM module init\n");
  nc= nc_pm;
  pm_factor= pm_factor_;
  nzpad= 2*(nc/2 + 1);
  boxsize= boxsize_;

  const size_t nckz= nc/2 + 1;
  
  mem_use_from_zero(mem_pm, 0);
  fft_pm= fft_alloc("PM", nc, mem_pm, 1);


  delta_k= mem_use_from_zero(mem_density,
			     nc*(fft_pm->local_nky)*nckz*sizeof(complex_t));
  //assert(mem_pm != mem_density);
  //assert(mem_pm->buf != mem_density->buf);
  //assert(mem_pm->buf == fft_pm->fk);
  //assert(delta_k != fft_pm->fk);
}

void pm_compute_forces(Particles* particles)
{
  // Main routine of this source file
  msg_printf(msg_verbose, "PM force computation...\n");
  
  size_t np_plus_buffer= send_buffer_positions(particles);
  pm_assign_cic_density(particles, np_plus_buffer);
  check_total_density(fft_pm->fx);

  compute_delta_k();

  for(int axis=0; axis<3; axis++) {
    // delta(k) -> f(x_i)
    compute_force_mesh(axis);

    force_at_particle_locations(particles, np_plus_buffer, axis);


    //force_at_particle_locations(particles->p, np_plus_buffer, axes,
    //(float*) fftdata, particles->force);
  }
  add_buffer_forces(particles, np_plus_buffer);
}


//
// Private (static) functions
//

size_t send_buffer_positions(Particles* const particles)
{
  assert(boxsize > 0);
  // ToDo !!! send with MPI !!!
  const size_t np= particles->np_local;
  Particle* const p= particles->p;
  const int nbuf= particles->np_allocated;
  
  const float_t eps= boxsize/nc;
  const float_t x_left= eps;
  const float_t x_right= boxsize - eps;


  size_t ibuf= np;

  // Periodic wrap up and make a buffer copy near x-edges
  for(size_t i=0; i<np; i++) {
    if(p[i].x[0] < 0) p[i].x[0] += boxsize;
    else if(p[i].x[0] >= boxsize) p[i].x[0] -= boxsize;

    if(p[i].x[1] < 0) p[i].x[1] += boxsize;
    else if(p[i].x[1] >= boxsize) p[i].x[1] -= boxsize;
    
    if(p[i].x[2] < 0) p[i].x[2] += boxsize;
    else if(p[i].x[2] >= boxsize) p[i].x[2] -= boxsize;
    
    if(p[i].x[0] < x_left) {
      if(ibuf >= nbuf)
	msg_abort("Error: not enough space for buffer particles. "
		  "%lu %lu\n", ibuf, nbuf);
      
      p[ibuf].x[0]= p[i].x[0] + boxsize;
      p[ibuf].x[1]= p[i].x[1];
      p[ibuf].x[2]= p[i].x[2];
      p[ibuf].id= p[i].id;

      ibuf++;
    }
    else if(p[i].x[0] > x_right) {
      if(ibuf >= nbuf)
	msg_abort("Error: not enough space for buffer particles. "
		  "%ud %ud\n", ibuf, nbuf);

      p[ibuf].x[0]= p[i].x[0] - boxsize;
      p[ibuf].x[1]= p[i].x[1];
      p[ibuf].x[2]= p[i].x[2];
      p[ibuf].id= p[i].id;
      ibuf++;
    }


#ifdef CHECK
    if(!(p[i].x[0] >= 0 && p[i].x[0] <= boxsize))
      printf("%lu %e\n", i, p[i].x[0]);
    assert(p[i].x[0] >= 0 && p[i].x[0] <= boxsize);
    assert(p[i].x[1] >= 0 && p[i].x[1] <= boxsize);
    assert(p[i].x[2] >= 0 && p[i].x[2] <= boxsize);
#endif

  }

  return ibuf;
}


  
void pm_assign_cic_density(Particles* particles, size_t np) 
{
  // Input:  particle positions in particles->p.x
  // Result: density field delta(x) in fft_pm->fx

  // particles are assumed to be periodiclly wraped up in y,z direction
  // and np is the number of particles including buffer particles
  

  float_t* const density= (float*) fft_pm->fx;
  Particle* p= particles->p;
  const size_t local_nx= fft_pm->local_nx;
  const size_t local_ix0= fft_pm->local_ix0;
 
  msg_printf(msg_verbose, "particle position -> density mesh\n");

  const float_t dx_inv= nc/boxsize;
  
  const float_t fac= pm_factor*pm_factor*pm_factor;


#ifdef _OPENMP
  #pragma omp parallel for default(shared)
#endif
  for(size_t ix = 0; ix < local_nx; ix++)
    for(size_t iy = 0; iy < nc; iy++)
      for(size_t iz = 0; iz < nc; iz++)
	density[(ix*nc + iy)*nzpad + iz] = -1;

#ifdef _OPENMP
  #pragma omp parallel for default(shared)
#endif
  for(size_t i=0; i<np; i++) {
    float x=p[i].x[0]*dx_inv;
    float y=p[i].x[1]*dx_inv;
    float z=p[i].x[2]*dx_inv;

    //iI, J, K
    int ix0= (int) floorf(x); // without floor, -1 < X < 0 is mapped to iI=0
    int iy0= (int) y;         // assuming y,z are positive
    int iz0= (int) z;

    // CIC weight on left grid
    float_t wx1= x - ix0; // D1 D2 D3
    float_t wy1= y - iy0;
    float_t wz1= z - iz0;

    // CIC weight on right grid
    float_t wx0= 1 - wx1;    // T1 T2 T3
    float_t wy0= 1 - wy1;
    float_t wz0= 1 - wz1;

    //float_tfloat T2W =T2*WPAR;
    //float D2W =D2*WPAR;

#ifdef CHECK
    assert(y >= 0.0f && z >= 0.0f);
#endif
            
    // No periodic wrapup in x direction. 
    // Buffer particles are copied from adjacent nodes, instead
    if(iy0 >= nc) iy0= 0; 
    if(iz0 >= nc) iz0= 0;

    // I1 J1 K1
    int ix1= ix0 + 1;
    int iy1= iy0 + 1; if(iy1 >= nc) iy1= 0; // assumes y,z < boxsize
    int iz1= iz0 + 1; if(iz1 >= nc) iz1= 0;

    ix0 -= local_ix0;
    ix1 -= local_ix0;

    if(0 <= ix0 && ix0 < local_nx) {
      grid_assign(density, ix0, iy0, iz0, fac*wx0*wy0*wz0); //T3*T1*T2W
      grid_assign(density, ix0, iy0, iz1, fac*wx0*wy0*wz1); //D3*T1*T2W);
      grid_assign(density, ix0, iy1, iz0, fac*wx0*wy1*wz0); //T3*T1*D2W);
      grid_assign(density, ix0, iy1, iz1, fac*wx0*wy1*wz1); //D3*T1*D2W);
    }

    if(0 <= ix1 && ix1 < local_nx) {
      grid_assign(density, ix1, iy0, iz0, fac*wx1*wy0*wz0); // T3*D1*T2W);
      grid_assign(density, ix1, iy0, iz1, fac*wx1*wy0*wz1); // D3*D1*T2W);
      grid_assign(density, ix1, iy1, iz0, fac*wx1*wy1*wz0); // T3*D1*D2W);
      grid_assign(density, ix1, iy1, iz1, fac*wx1*wy1*wz1); // D3*D1*D2W);
    }
  }

  /*
  FILE* fp= fopen("particle.txt", "w");
  for(size_t i = 0; i < particles->np_local; i++)
    fprintf(fp, "%e %e %e\n", p[i].x[0], p[i].x[1], p[i].x[2]);
  fclose(fp);
  printf("particle.txt written\n");
  */


  /*
  fp= fopen("density.txt", "w");
  for(size_t ix = 0; ix < local_nx; ix++)
   for(size_t iy = 0; iy < nc; iy++)
    for(size_t iz = 0; iz < nc; iz++)
      fprintf(fp, "%e\n", density[(ix*nc + iy)*nzpad + iz]);
  fclose(fp);
  printf("density.txt written\n");
  abort();
  */
  
  msg_printf(msg_verbose, "CIC density assignment finished.\n");
}

void check_total_density(float_t const * const density)
{
  // Checks <delta> = 0
  // Input: delta(x)
  double sum= 0.0;
  const size_t local_nx= fft_pm->local_nx;
  
  for(size_t ix = 0; ix < local_nx; ix++)
    for(size_t iy = 0; iy < nc; iy++)
      for(size_t iz = 0; iz < nc; iz++)
	sum += density[(ix*nc + iy)*nzpad + iz];

  double sum_global;
  MPI_Reduce(&sum, &sum_global, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if(comm_this_node() == 0) {
    double tol= FLOAT_EPS*nc*nc*nc;

    if(fabs(sum_global) > tol)
      msg_abort("Error: total CIC density error is large: %le > %le\n", 
		sum_global, tol);

    msg_printf(msg_debug, 
	      "Total CIC density OK within machine precision: %lf (< %.2lf).\n",
	       sum_global, tol);
  }
}


void compute_delta_k(void)
{
  // Fourier transform delta(x) -> delta(k) and copy it to delta_k
  //  Input:  delta(x) in fft_pm->fx
  //  Output: delta(k) in delta_k

  msg_printf(msg_verbose, "delta(x) -> delta(k)\n");
  fft_execute_forward(fft_pm);

  // Copy density(k) in fft_pm to density_k
  // because FFT requires twice larger RAM for working memory
  const size_t nckz= nc/2 + 1;
  const size_t local_nky= fft_pm->local_nky;

  complex_t* pm_k= fft_pm->fk;
  
#ifdef _OPENMP
  #pragma omp parallel for default(shared)
#endif
  for(size_t iy=0; iy<local_nky; iy++) {
    for(size_t ix=0; ix<nc; ix++) {
      for(size_t iz=0; iz<nckz; iz++){
	size_t index= (nc*iy + ix)*nckz + iz;
	delta_k[index][0]= pm_k[index][0];
	delta_k[index][1]= pm_k[index][1];
	
	//printf("%e %e\n", pm_k[index][0], pm_k[index][1]);
      }
    }
  }
}

void compute_force_mesh(const int axis)
{
  // Calculate one component of force mesh from precalculated density(k)
  //   Input:   delta(k)   mesh delta_k
  //   Output:  force_i(k) mesh fft_pm->fx

  complex_t* const fk= fft_pm->fk;
  
  //k=0 zero mode force is zero
  fk[0][0]= 0;
  fk[0][1]= 0;

  const float_t f1= -1.0/pow(nc, 3.0)/(2.0*M_PI/boxsize);
  const size_t nckz=nc/2+1;
  const size_t local_nky= fft_pm->local_nky;
  const size_t local_iky0= fft_pm->local_iky0;


#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
  for(size_t iy_local=0; iy_local<local_nky; iy_local++) {
    int iy= iy_local + local_iky0;
    int iy0= iy <= (nc/2) ? iy : iy - nc;

    float_t k[3];
    k[1]= (float_t) iy0;

    for(size_t ix=0; ix<nc; ix++) {
      int ix0= ix <= (nc/2) ? ix : ix - nc;
      k[0]= (float_t) ix0;

      int kzmin= (ix==0 && iy==0); // skip (0,0,0) to avoid zero division

      for(size_t iz=kzmin; iz<nckz; iz++){
	k[2]= (float_t) iz;

	float f2= f1/(k[0]*k[0] + k[1]*k[1] + k[2]*k[2])*k[axis];

	size_t index= (nc*iy_local + ix)*nckz + iz;
	fk[index][0]= -f2*delta_k[index][1];
	fk[index][1]=  f2*delta_k[index][0];
      }
    }
  }


  /*
  FILE* fp= fopen("fk.txt", "w");
  for(size_t iy_local=0; iy_local<local_nky; iy_local++) {
    for(size_t ix=0; ix<nc; ix++) {
      for(size_t iz=0; iz<nckz; iz++){
	size_t index= (nc*iy_local + ix)*nckz + iz;       
	fprintf(fp, "%e %e %e %e\n",
		delta_k[index][0], delta_k[index][1],
		fk[index][0], fk[index][1]);
     }
    }
  }
  fclose(fp);
  printf("fk.txt written\n");
  abort();
  */

  fft_execute_inverse(fft_pm); // f_k -> f(x)
}

// Does 3-linear interpolation
// particles= Values of mesh at particle positions P.x
void force_at_particle_locations(Particles* const particles, const int np, 
				 const int axis)
{
  const Particle* p= particles->p;
  
  const float_t dx_inv= nc/boxsize;
  const size_t local_nx= fft_pm->local_nx;
  const size_t local_ix0= fft_pm->local_ix0;
  const float_t* fx= fft_pm->fx;
  float3* f= particles->force;
  
#ifdef _OPENMP
  #pragma omp parallel for default(shared)     
#endif
  for(size_t i=0; i<np; i++) {
    float_t x=p[i].x[0]*dx_inv;
    float_t y=p[i].x[1]*dx_inv;
    float_t z=p[i].x[2]*dx_inv;
            
    int ix0= (int) floorf(x); /// iI !!! floor for double?
    int iy0= (int) y; //J
    int iz0= (int) z;   //K
    
    float_t wx1= x - (float) ix0; //D1
    float_t wy1= y - (float) iy0; //D2
    float_t wz1= z - (float) iz0; //D3

    float_t wx0= 1 - wx1; // T1
    float_t wy0= 1 - wy1; // T2
    float_t wz0= 1 - wz1; // T3

    if(iy0 >= nc) iy0= 0;
    if(iz0 >= nc) iz0= 0;
            
    int ix1= ix0 + 1;  // I1
    int iy1= iy0 + 1; if(iy1 >= nc) iy1= 0; // J1
    int iz1= iz0 + 1; if(iz1 >= nc) iz1= 0; // K1

    ix0 -= local_ix0;
    ix1 -= local_ix0;

    f[i][axis]= 0;

    if(0 <= ix0 && ix0 < local_nx) {
      f[i][axis] += 
	grid_val(fx, ix0, iy0, iz0)*wx0*wy0*wz0 +
	grid_val(fx, ix0, iy0, iz1)*wx0*wy0*wz1 +
	grid_val(fx, ix0, iy1, iz0)*wx0*wy1*wz0 +
	grid_val(fx, ix0, iy1, iz1)*wx0*wy1*wz1;
    }
    if(0 <= ix1 && ix1 < local_nx) {
      f[i][axis] += 
	grid_val(fx, ix1, iy0, iz0)*wx1*wy0*wz0 +
	grid_val(fx, ix1, iy0, iz1)*wx1*wy0*wz1 +
	grid_val(fx, ix1, iy1, iz0)*wx1*wy1*wz0 +
	grid_val(fx, ix1, iy1, iz1)*wx1*wy1*wz1;
    }
  }
  
}

void add_buffer_forces(Particles* const particles, const size_t np)
{
  // !! Non-MPI version
  Particle* const p= particles->p;
  
  const int np_local= particles->np_local;
  float3* const force= particles->force;

  for(int j=np_local; j<np; j++) {
    size_t i= p[j].id - 1;
    assert(p[i].id == p[j].id);
    force[i][0] += force[j][0];
    force[i][1] += force[j][1];
    force[i][2] += force[j][2];
  }

  /*
  FILE* fp= fopen("force.txt", "w");
  for(size_t i = 0; i < particles->np_local; i++)
    fprintf(fp, "%e %e %e\n", force[i][0], force[i][1], force[i][2]);
  fclose(fp);
  printf("force.txt written\n");
  abort();
  */

}



