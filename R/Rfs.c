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
#include "util.h"

static bool initialised= false;
static PowerSpectrum* ps= 0;
static Particles* particles= NULL;
static int nc;
static float_t boxsize;

void rfs_init(void)
{
  comm_mpi_init(0,0);
  msg_set_loglevel(msg_debug);

  const double omega_m0= 0.273;
  cosmology_init(omega_m0);
}
  
  
void rfs_read_powerspectrum(char** filename)
{
  //printf("-%s-\n", *filename);

  if(ps) power_free(ps);
  ps= power_alloc(*filename, 0);
}

void rfs_lpt_init(int* nc_, double* boxsize_)
{
  assert(ps);
  nc= *nc_;
  boxsize= *boxsize_;
  
  lpt_init(nc, boxsize, 0);

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
  int a=0;
  lpt_set_displacements(seed, ps, a, particles);
  
}

void rfs_set_LPT(int* seed, double* a)
{
  lpt_set_displacements(*seed, ps, *a, particles);
}

SEXP rfs_particles_with_id(SEXP indices)
{
  // Return a data.frame of particles for given index
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

SEXP rfs_particles_in_cuboid(SEXP cuboid)
{
  // Return a data.frame of particles in the given cuboid
  // cuboid= NULL      => everything in the box
  // cuboid= c(z0, z1) : cuboid is a slice for if cuboid has only 2 indeces
  // cuboid= c(x0, x1, y0, y1, z0, z1)
  float_t left[]= {0,0,0};
  float_t right[]= {boxsize, boxsize, boxsize};
  
  const int n= length(cuboid);
  double const * const arg= REAL(cuboid);
  if(n == 0)
    ;
  else if(n == 2) {
    left[2]= arg[0];
    right[2]= arg[1];
  }
  else if(n == 6) {
    for(int i=0; i<3; i++) {
      left[i]= arg[2*i];
      right[i]= arg[2*i + 1];
    }
  }
  else {
    error("Error: number of cuboid argument must be 2 or 6\n");
    return R_NilValue;
  }

  msg_printf(msg_verbose, "particles in cuboid %.1f %.1f %.1f %.1f %.1f %.1f\n",
	 left[0], right[0], left[1], right[1], left[2], right[2]);

  // count number of particles in cuboid
  size_t count= 0;
  const size_t np= particles->np_local;
  Particle* p= particles->p;
  for(size_t i=0; i<np; i++) {
    periodic_wrapup_p(p + i, boxsize);
    if(left[0] <= p[i].x[0] && p[i].x[0] < right[0] &&
       left[1] <= p[i].x[1] && p[i].x[1] < right[1] &&
       left[2] <= p[i].x[2] && p[i].x[2] < right[2])
      //printf("%e %e %e\n", p[i].x[0], p[i].x[1], p[i].x[2]);
      count++;
  }

  const int ncol= 4; // id, x, y, z
  const int nrow= count;
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

  //Particle const * const p= particles->p;
  size_t index= 0;
  for(int i=0; i<np; i++) {
    if(left[0] <= p[i].x[0] && p[i].x[0] < right[0] &&
       left[1] <= p[i].x[1] && p[i].x[1] < right[1] &&
       left[2] <= p[i].x[2] && p[i].x[2] < right[2]) {

      rx[index]= p[i].x[0];
      ry[index]= p[i].x[1];
      rz[index]= p[i].x[2];
      ii[index]= (int) p[i].id;

      index++;
    }
  }

  assert(index == count);
  
  // set the vectors in list$x <- x ...
  SET_VECTOR_ELT(list, 0, id);
  SET_VECTOR_ELT(list, 1, x);
  SET_VECTOR_ELT(list, 2, y);
  SET_VECTOR_ELT(list, 3, z);

  UNPROTECT(7);
  return list;
}
