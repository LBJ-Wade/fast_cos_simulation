#include <assert.h>
#include <stdbool.h>
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include "msg.h"
#include "comm.h"
#include "power.h"
#include "cosmology.h"
#include "lpt.h"

static bool initialised= false;
static PowerSpectrum* ps= 0;
static Particles* particles= NULL;

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
  if(particles == NULL) {
    particles= malloc(sizeof(Particles));
    particles->p= NULL;
  }

  size_t np_alloc= (int)(1.25*nc*nc*nc);


  particles->np_allocated= np_alloc;
  
  if(particles->p && particles->np_allocated != np_alloc)
    free(particles->p);
  else
    particles->p= malloc(sizeof(Particle)*np_alloc);


  // ToDo arguments
  int seed= 1;
  int a=1;
  lpt_set_displacements(seed, ps, a, particles);
  
}

SEXP rfs_particles_with_id(SEXP indices)
{
  // Retrun a data.frame of particles for given index
  // ToDo: change to id not index

  const int n= length(indices);
  int const * const index= INTEGER(indices);

  const int ncol= 4; // id, x, y, z
  const int nrow= particles->np_local < n ? particles->np_local : n;
  printf("nrow= %d\n", nrow);

  SEXP list, col_names, row_names;
  PROTECT(list = allocVector(VECSXP, ncol));
  PROTECT(col_names = allocVector(STRSXP, ncol));
  PROTECT(row_names = allocVector(STRSXP, nrow));
  
  // column names
  SET_STRING_ELT(col_names, 0, mkChar("id"));
  SET_STRING_ELT(col_names, 1, mkChar("x"));
  SET_STRING_ELT(col_names, 2, mkChar("y"));
  SET_STRING_ELT(col_names, 3, mkChar("z"));
  setAttrib(list, R_NamesSymbol, col_names);

  // row names 1,2,3, ...
  char rname[20];
  for(int i=0; i<nrow; i++) {
    sprintf(rname, "%d", i+1);
    SET_STRING_ELT(row_names, i, mkChar(rname));
  }
  setAttrib(list, R_RowNamesSymbol, row_names);
    
  // class=data.frame
  setAttrib(list, R_ClassSymbol, ScalarString(mkChar("data.frame")));
  

  SEXP x,y,z, id;
  PROTECT(x = allocVector(REALSXP, nrow));
  PROTECT(y = allocVector(REALSXP, nrow));
  PROTECT(z = allocVector(REALSXP, nrow));
  PROTECT(id = allocVector(INTSXP, nrow));

  double* const rx= REAL(x);
  double* const ry= REAL(y);
  double* const rz= REAL(z);
  int* const ii= INTEGER(id);

  Particle const * const p= particles->p;
  for(int i=0; i<nrow; i++) {
    rx[i]= p[index[i]].x[0];
    ry[i]= p[index[i]].x[1];
    rz[i]= p[index[i]].x[2];
    ii[i]= (int) p[index[i]].id;
  }

  // set the vectors in list$x <- x ...
  SET_VECTOR_ELT(list, 0, id);
  SET_VECTOR_ELT(list, 1, x);
  SET_VECTOR_ELT(list, 2, y);
  SET_VECTOR_ELT(list, 3, z);

  UNPROTECT(7);
  return list;
}
