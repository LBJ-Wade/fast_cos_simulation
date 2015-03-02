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
  //int nc_pm= 3*nc;

  const float a_init= 1.0;

  // Memory management
  Mem* mem1= mem_init("mem1"); // mainly for density
  mem_reserve(mem1, 9*fft_mem_size(nc), "LPT");
  mem_alloc_reserved(mem1);
  
  Particles* particles= alloc_particles(nc);


  // 2LPT initial condition / displacement
  const double omega_m= 0.273;
  cosmology_init(omega_m);
  lpt_init(nc, boxsize, mem1);
  lpt_set_displacements(seed, ps, a_init, particles);

     
  msg_printf(msg_info, "Hello World\n");

  comm_mpi_finalise();
}

Particles* alloc_particles(const int nc)
{
  Particles* particles= calloc(sizeof(Particles), 1); assert(particles);
  size_t np_alloc= (size_t)((1.25*nc/comm_n_nodes())*nc*nc);
  particles->p= malloc(np_alloc*sizeof(Particle)); assert(particles->p);

  particles->np_allocated= np_alloc;

  msg_printf(msg_verbose, "%lu Mbytes allocated for particles\n",
	     mbytes(np_alloc));

  return particles;
}
