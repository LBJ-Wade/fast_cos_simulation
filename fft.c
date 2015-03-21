///
/// \file  fft.c
/// \brief interface for FFTW
///

#include <stdlib.h>
#include <assert.h>
#include <fftw3-mpi.h>
#include "config.h"
#include "mem.h"
#include "msg.h"
#include "util.h"
#include "fft.h"

// FFTW() adds fftw_ or fftwf_ prefix depending on DOUBLEPRECISION

#ifdef DOUBLEPRECISION
#define FFTW(f) fftw_ ## f
#else
#define FFTW(f) fftwf_ ## f
#endif


#ifdef MPI

FFT* fft_alloc(const char name[], const int nc, Mem* mem, const int transposed)
{
  // Allocates memory for FFT real and Fourier space and initilise fftw_plans

  //msg_printf(msg_debug, "fft_alloc(%s, nc=%d)\n", name, nc);
  //msg_printf(msg_debug, "malloc(%lu)\n", sizeof(FFT));
  FFT* const fft= malloc(sizeof(FFT)); assert(fft);
  fft->nc= nc;

  //msg_printf(msg_debug, "fft_alloc nc= %d\n", nc);

  ptrdiff_t ncomplex= 0;
  if(transposed) {
    ncomplex= FFTW(mpi_local_size_3d_transposed)(nc, nc, nc, MPI_COMM_WORLD,
	                 &fft->local_nx, &fft->local_ix0,
			 &fft->local_nky, &fft->local_iky0);
  }
  else {
    ncomplex= FFTW(mpi_local_size_3d)(nc, nc, nc, MPI_COMM_WORLD,
			    &fft->local_nx, &fft->local_ix0);
    fft->local_nky= fft->local_iky0= 0;
  }

  size_t size= sizeof(complex_t)*ncomplex;
  assert(fft->local_nx >= 0); assert(fft->local_ix0 >= 0);
  
    
  if(mem == 0)
    mem= mem_alloc(name, size);
  
  void* buf= mem_use_remaining(mem, size);
    // Call mem_use_from_zero(mem, 0) before this to use mem from the beginning.

  fft->ncomplex= ncomplex;

  fft->fx= buf; fft->fk= buf;

  unsigned flag= 0;
  if(transposed) flag= FFTW_MPI_TRANSPOSED_OUT;
  fft->forward_plan= FFTW(mpi_plan_dft_r2c_3d)(nc, nc, nc, fft->fx, fft->fk,
				       MPI_COMM_WORLD, FFTW_MEASURE | flag);

  unsigned flag_inv= 0;
  if(transposed) {
    flag_inv= FFTW_MPI_TRANSPOSED_IN;
    msg_printf(msg_debug, "FFTW transposed in/out\n");
  }
  
  fft->inverse_plan= FFTW(mpi_plan_dft_c2r_3d)(nc, nc, nc, fft->fk,fft->fx,
                                     MPI_COMM_WORLD, FFTW_MEASURE | flag_inv);

  // ToDo: FFTW_MPI_TRANSPOSED_IN/FFTW_MPI_TRANSPOSED_OUT would be faster

  return fft;
}

size_t fft_mem_size_working(const int nc, const int transposed)
{
  // return the memory size necessary for the 3D FFT
  ptrdiff_t local_nx, local_ix0, local_nky, local_iky0;

  ptrdiff_t ncomplex= 0;
  if(transposed)
    ncomplex= FFTW(mpi_local_size_3d_transposed)(nc, nc, nc, MPI_COMM_WORLD,
	           &local_nx, &local_ix0, &local_nky, &local_iky0);
  else
    ncomplex= FFTW(mpi_local_size_3d)(nc, nc, nc, MPI_COMM_WORLD,
				      &local_nx, &local_ix0);


  //printf("ncomplex %d %ld %lu %lu\n", nc, ncomplex, sizeof(complex_t), size_align(sizeof(complex_t)*ncomplex));
  
  return size_align(sizeof(complex_t)*ncomplex);
}

size_t fft_mem_size_fk(const int nc, const int transposed)
{
  // return the memory size necessary for the 3D FFT data in k space
  ptrdiff_t local_nx, local_ix0, local_nky, local_iky0;

  if(transposed)
    FFTW(mpi_local_size_3d_transposed)(nc, nc, nc, MPI_COMM_WORLD,
	           &local_nx, &local_ix0, &local_nky, &local_iky0);
  else {
    FFTW(mpi_local_size_3d)(nc, nc, nc, MPI_COMM_WORLD,
				      &local_nx, &local_ix0);
    local_nky= local_nx;
  }
  msg_printf(msg_verbose,
	     "FFT_mem_size_fk %lu %ld\n", sizeof(complex_t)*nc*nc*local_nky,
	     local_nky);
  
  const size_t nckz= nc/2 + 1;
  
  return size_align(sizeof(complex_t)*nc*local_nky*nckz);
}

size_t fft_local_nx(const int nc)
{
  ptrdiff_t local_nx, local_ix0; 
  FFTW(mpi_local_size_3d)(nc, nc, nc, MPI_COMM_WORLD, &local_nx, &local_ix0);

  return local_nx;
}


void fft_finalize(void)
{
  FFTW(mpi_cleanup)();
}

void fft_execute_forward(FFT* const fft)
{
  FFTW(mpi_execute_dft_r2c)(fft->forward_plan, fft->fx, fft->fk);
}

void fft_execute_inverse(FFT* const fft)
{
  FFTW(mpi_execute_dft_c2r)(fft->inverse_plan, fft->fk, fft->fx);
}


#else
// Serial version

FFT* fft_alloc(const char name[], const int nc, Mem* mem, unsigned flags)
{
  FFT* const fft= malloc(sizeof(FFT)); assert(fft);
  fft->nc= nc;
  fft->local_nx= nc;
  fft->local_ix0= 0;

  const size_t nckz= nc/2 + 1;
  ptrdiff_t ncomplex= nc*nc*nckz;
    
  size_t size= sizeof(complex_t)*ncomplex;
    
  if(mem == 0)
    mem= mem_alloc(name, size);

  fft->ncomplex= ncomplex;
  void* buf= mem_use_remaining(mem, size);
  fft->fx= buf; fft->fk= buf;

  unsigned flag0= FFTW_ESTIMATE; // can use FFTW_MEASURE for many realizations
  // serial version are mainly used for interactive jobs
  // small overhead with FFTW_ESTIMATE is probablly better than FFTW_MEASURE
  
  fft->forward_plan= FFTW(plan_dft_r2c_3d)(nc, nc, nc, fft->fx, fft->fk,
					   flag0  | flags);
  fft->inverse_plan= FFTW(plan_dft_c2r_3d)(nc, nc, nc, fft->fk, fft->fx,
					   flag0 | flags);

  return fft;
}

void fft_finalize(void)
{
  FFTW(cleanup)();
}

void fft_execute_forward(FFT* const fft)
{
  FFTW(execute_dft_r2c)(fft->forward_plan, fft->fx, fft->fk);
}

void fft_execute_inverse(FFT* const fft)
{
  FFTW(execute_dft_c2r)(fft->inverse_plan, fft->fk, fft->fx);
}

#endif

// No change whehter with or without MPI

void fft_free(FFT* const fft)
{
  FFTW(destroy_plan)(fft->forward_plan);
  FFTW(destroy_plan)(fft->inverse_plan);
  
  if(fft->allocated == true) {
    free(fft->fx);
    fft->fx= NULL; fft->fk= NULL;
  }
}

void* fft_malloc(size_t size)
{
  return FFTW(malloc)(size);
}

// Quote
// "it is probably better for you to simply create multiple plans
//  (creating a new plan is quick once one exists for a given size)
// -- FFTW3 Manual for version 3.3.3 section 4.6

// "To prevent memory leaks, you must still call fftw_destroy_plan
//  before executing fftw_cleanup."
