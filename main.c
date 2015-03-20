#include <stdio.h>
#include <assert.h>

#include "config.h"
#include "particle.h"
#include "util.h"
#include "comm.h"
#include "msg.h"
#include "power.h"
#include "mem.h"
#include "fft.h"
#include "cosmology.h"
#include "lpt.h"
#include "cola.h"
#include "pm.h"

Particles* alloc_particles(const int nc);

int main(int argc, char* argv[])
{
  // Setup MPI Init  as comm_mpi_init?
  //
  comm_mpi_init(&argc, &argv);
  msg_set_loglevel(msg_debug);
  
  PowerSpectrum* ps= power_alloc("camb_matterpower.dat", 0.812);

  // Parameters
  int nc= 64;
  float boxsize= 1000.0f;
  const unsigned long seed= 100;
  const double omega_m= 0.273;

  const int nstep= 10;

  const int pm_factor= 3;
  const int nc_pm= pm_factor*nc;
  const float a_init= 1.0/nstep;

  // Memory management
  Mem* mem1= mem_init("mem1"); // mainly for density
  mem_reserve(mem1, 9*fft_mem_size(nc), "LPT");
  mem_reserve(mem1, fft_mem_size(nc_pm), "PM-density");
  mem_alloc_reserved(mem1);

  Mem* mem2= mem_init("mem2");
  mem_reserve(mem2, fft_mem_size(nc_pm), "PM-force");
  mem_alloc_reserved(mem2);
  
  Particles* particles= alloc_particles(nc);
  particles->omega_m= omega_m;
  
  // 2LPT initial condition / displacement

  cosmology_init(omega_m);
  lpt_init(nc, boxsize, mem1);
  pm_init(nc_pm, pm_factor, mem1, mem2, boxsize);

  lpt_set_displacements(seed, ps, a_init, particles);
  particles->a_v= 1.0/nstep; // origial a_v

  for(int istep=1; istep<=nstep; istep++) {
    float_t a_vel= (istep + 0.5)/nstep;
    float_t a_pos= (istep + 1.0)/nstep;

    pm_compute_forces(particles);
    cola_kick(particles, a_vel);
    cola_drift(particles, a_pos);
  }
     
  msg_printf(msg_info, "Hello World\n");

  comm_mpi_finalise();
}

Particles* alloc_particles(const int nc)
{
  Particles* particles= calloc(sizeof(Particles), 1); assert(particles);
  size_t np_alloc= (size_t)((1.25*nc/comm_n_nodes())*nc*nc);
  particles->p= malloc(np_alloc*sizeof(Particle)); assert(particles->p);
  particles->force= calloc(3*np_alloc, sizeof(float)); assert(particles->force);


  particles->np_allocated= np_alloc;

  msg_printf(msg_verbose, "%lu Mbytes allocated for %lu particles\n",
	     mbytes(np_alloc), np_alloc);

  return particles;
}
