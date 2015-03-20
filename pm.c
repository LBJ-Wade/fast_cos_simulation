#include <math.h>
#include <assert.h>
#include <gsl/gsl_rng.h>
#include "msg.h"
#include "mem.h"
#include "config.h"
#include "cosmology.h"
//#include "power.h"
#include "comm.h"
#include "particle.h"
#include "fft.h"

static int pm_factor;
static size_t nc, nzpad;
static float_t boxsize;

static FFT* fft_density; 
static FFT* fft_force;

static inline void grid_assign(float_t * const d, 
	    const size_t ix, const size_t iy, const size_t iz, const float_t f)
{
#ifdef _OPENMP
  #pragma omp atomic
#endif
  d[(ix*nc + iy)*nzpad + iz] += f;
}

static size_t send_buffer_positions(Particles* const particles);
static void pm_assign_cic_density(Particles* particles, size_t np);
static void check_total_density(float_t const * const density);

void pm_init(const int nc_pm, const int pm_factor_,
	     Mem* const mem_density, Mem* const mem_force,
	     const float_t boxsize_)
{
  nc= nc_pm;
  pm_factor= pm_factor_;
  nzpad= 2*(nc/2 + 1);
  boxsize= boxsize_;
  
  mem_use_from_zero(mem_density, 0);
  mem_use_from_zero(mem_force, 0);
  
  fft_density= fft_alloc("PM-density", nc, mem_density, 1);
  fft_force= fft_alloc("PM-force", nc, mem_force, 1);

  
}

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
  float_t zero= 0;

  // Periodic wrap up and make a buffer copy near x-edges
  for(size_t i=0; i<np; i++) {
    if(p[i].x[1] < zero) p[i].x[1] += boxsize;
    else if(p[i].x[1] >= boxsize) p[i].x[1] -= boxsize;
    
    if(p[i].x[2] < zero) p[i].x[2] += boxsize;
    else if(p[i].x[2] >= boxsize) p[i].x[2] -= boxsize;
    
    if(p[i].x[0] < zero) p[i].x[0] += boxsize;
    else if(p[i].x[0] >= boxsize) p[i].x[0] -= boxsize;
    
    if(p[i].x[0] < x_left) {
      if(ibuf >= nbuf)
	msg_abort("Error: not enough space for buffer particles. "
		  "%lu %lu\n", ibuf, nbuf);
      
      p[ibuf].x[0]= p[i].x[0] + boxsize;
      p[ibuf].x[1]= p[i].x[1];
      p[ibuf].x[2]= p[i].x[2];
      p[ibuf].id= 0;

      ibuf++;
    }
    else if(p[i].x[0] > x_right) {
      if(ibuf >= nbuf)
	msg_abort("Error: not enough space for buffer particles. "
		  "%ud %ud\n", ibuf, nbuf);

      p[ibuf].x[0]= p[i].x[0] - boxsize;
      p[ibuf].x[1]= p[i].x[1];
      p[ibuf].x[2]= p[i].x[2];
      p[ibuf].id= 0;
      ibuf++;
    }


#ifdef CHECK
    assert(p[i].x[0] >= 0 && p[i].x[0] <= boxsize);
    assert(p[i].x[1] >= 0 && p[i].x[1] <= boxsize);
    assert(p[i].x[2] >= 0 && p[i].x[2] <= boxsize);
#endif

  }

  return ibuf;
}


  
void pm_assign_cic_density(Particles* particles, size_t np) 
{
  // particles are assumed to be periodiclly wraped up in y,z direction
  // and buffer particles are
  // using fft_density
  
  float_t* const density= (float*) fft_density->fx;
  Particle* p= particles->p;
  const size_t local_nx= fft_density->local_nx;
  const size_t local_ix0= fft_density->local_ix0;
 
  msg_printf(msg_verbose, "particle position -> density mesh\n");

  //const float_t dx= Ngrid/boxsize
  //const float scaleBox=((float) Ngrid)/((float) BoxSize);
  //const float WPAR= pow(PM_factor, 3);
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

  msg_printf(msg_verbose, "CIC density assignment finished.\n");
}

void check_total_density(float_t const * const density)
{
  double sum= 0.0;
  const size_t local_nx= fft_density->local_nx;
  
  for(int ix = 0; ix < local_nx; ix++)
    for(int iy = 0; iy < nc; iy++)
      for(int iz = 0; iz < nc; iz++)
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


void compute_density_k(void)
{
  // FFT density(x) mesh -> density(k)
  fft_execute_forward(fft_density);
 
  // copy density(k) in fftdata to density_k
  // why copying???
  const size_t nckz= nc/2 + 1;
  const size_t local_ny= fft_density->local_nx;

  assert(local_ny == fft_force->local_nx);
  complex_t* density_k= fft_density->fk;
  complex_t* force_k= fft_force->fk;
  
#ifdef _OPENMP
  #pragma omp parallel for default(shared)
#endif
  for(size_t iy=0; iy<local_ny; iy++) {
    for(size_t ix=0; ix<nc; ix++) {
      for(size_t iz=0; iz<nckz; iz++){
	size_t index= iz + nckz*(ix + nc*iy);
	force_k[index][0]= density_k[index][0];
	force_k[index][1]= density_k[index][1];
      }
    }
  }
}

void pm_compute_forces(Particles* particles)
{
  size_t np_plus_buffer= send_buffer_positions(particles);
  pm_assign_cic_density(particles, np_plus_buffer);
  check_total_density(fft_density->fx);
}
