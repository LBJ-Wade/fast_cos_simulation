#include <assert.h>
#include <stdbool.h>
#include <R.h>
#include <Rinternals.h>
#include "msg.h"
#include "comm.h"
#include "power.h"
#include "cosmology.h"
#include "lpt.h"

static bool initialised= false;
static PowerSpectrum* ps= 0;
static Particles particles;

void rfs_init(void)
{
  comm_mpi_init(0,0);
  msg_set_loglevel(msg_debug);

  const double omega_m0= 0.273;
  cosmology_init(omega_m0);
}
  
  
void rfs_read_powerspectrum(char** filename)
{
  printf("-%s-\n", *filename);

  if(ps) power_free(ps);
  ps= power_alloc(*filename, 0);
}

void rfs_lpt_init(int* nc_, double *boxsize)
{
  assert(ps);
  int nc= *nc_;
  
  lpt_init(nc, *boxsize, 0);

  // Allocates memory for particle
  size_t np_alloc= (int)(1.25*nc*nc*nc);
  particles.p= malloc(sizeof(Particle)*np_alloc);
  particles.np_allocated= np_alloc;

  // ToDo arguments
  int seed= 1;
  int a=1;
  lpt_set_displacements(seed, ps, a, &particles);
  
}

