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

#ifdef MPI

FFT* fft_alloc(const char name[], const int nc, Mem* mem, unsigned flags)
{
  // Allocates memory for FFT real and Fourier space and initilise fftw_plans

  msg_printf(msg_debug, "fft_alloc(%s, nc=%d)", name, nc);
  msg_printf(msg_debug, "malloc(%lu)\n", sizeof(FFT));
  FFT* const fft= malloc(sizeof(FFT)); assert(fft);
  fft->nc= nc;

  msg_printf(msg_debug, "fft_alloc nc= %d\n", nc);

  ptrdiff_t ncomplex=
#ifdef DOUBLEPRECISION
    fftw_mpi_local_size_3d(nc, nc, nc, MPI_COMM_WORLD,
			   &fft->local_nx, &fft->local_ix0);    
#else
    fftwf_mpi_local_size_3d(nc, nc, nc, MPI_COMM_WORLD,
			    &fft->local_nx, &fft->local_ix0);
#endif

  size_t size= sizeof(complex_t)*ncomplex;
  assert(fft->local_nx >= 0); assert(fft->local_ix0 >= 0);
  
    
  if(mem == 0)
    mem= mem_alloc(name, size);
  else
    mem_use_remaining(mem, size);
    // Call mem_use_from_zero(mem, 0) before this to use mem from the beginning.

  fft->ncomplex= ncomplex;
  fft->fx= mem->buf; fft->fk= mem->buf;

#ifdef DOUBLEPRECISION  
  fft->forward_plan= fftw_mpi_plan_dft_r2c_3d(nc, nc, nc, fft->fx, fft->fk,
				       MPI_COMM_WORLD, FFTW_MEASURE | flags);
  fft->inverse_plan= fftw_mpi_plan_dft_c2r_3d(nc, nc, nc, fft->fk, fft->fx,
				       MPI_COMM_WORLD, FFTW_MEASURE | flags);
#else
  fft->forward_plan= fftwf_mpi_plan_dft_r2c_3d(nc, nc, nc, fft->fx, fft->fk,
				       MPI_COMM_WORLD, FFTW_MEASURE | flags);
  fft->inverse_plan= fftwf_mpi_plan_dft_c2r_3d(nc, nc, nc, fft->fk, fft->fx,
				       MPI_COMM_WORLD, FFTW_MEASURE | flags);
#endif  
  
  // ToDo: FFTW_MPI_TRANSPOSED_IN/FFTW_MPI_TRANSPOSED_OUT would be faster

  return fft;
}

void fft_finalize(void)
{
#ifdef DOUBLEPRECISION
  fftw_mpi_cleanup();
#else
  fftwf_mpi_cleanup();
#endif
}

void fft_execute_forward(FFT* const fft)
{
#ifdef DOUBLEPRECISION
  fftw_mpi_execute_dft_r2c(fft->forward_plan, fft->fx, fft->fk);
#else
  fftwf_mpi_execute_dft_r2c(fft->forward_plan, fft->fx, fft->fk);
#endif
}

void fft_execute_inverse(FFT* const fft)
{
#ifdef DOUBLEPRECISION
  fftw_mpi_execute_dft_c2r(fft->inverse_plan, fft->fk, fft->fx);
#else
  fftwf_mpi_execute_dft_c2r(fft->inverse_plan, fft->fk, fft->fx);
#endif
}


#else
// No MPI version

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
  else
    mem_use_remaining(mem, size);

  fft->ncomplex= ncomplex;
  fft->fx= mem->buf; fft->fk= mem->buf;

  unsigned flag0= FFTW_ESTIMATE; // can use FFTW_MEASURE for many realizations
  // serial version are mainly used for interactive jobs
  // small overhead with FFTW_ESTIMATE is probablly better than FFTW_MEASURE
  
#ifdef DOUBLEPRECISION  
  fft->forward_plan= fftw_plan_dft_r2c_3d(nc, nc, nc, fft->fx, fft->fk,
				          flag0  | flags);
  fft->inverse_plan= fftw_plan_dft_c2r_3d(nc, nc, nc, fft->fk, fft->fx,
					      flag0 | flags);
#else
  fft->forward_plan= fftwf_plan_dft_r2c_3d(nc, nc, nc, fft->fx, fft->fk,
					       flag0 | flags);
  fft->inverse_plan= fftwf_plan_dft_c2r_3d(nc, nc, nc, fft->fk, fft->fx,
					       flag0 | flags);
#endif  

  return fft;
}

void fft_finalize(void)
{
#ifdef DOUBLEPRECISION
  fftw_cleanup();
#else
  fftwf_cleanup();
#endif
}

void fft_execute_forward(FFT* const fft)
{
#ifdef DOUBLEPRECISION
  fftw_execute_dft_r2c(fft->forward_plan, fft->fx, fft->fk);
#else
  fftwf_execute_dft_r2c(fft->forward_plan, fft->fx, fft->fk);
#endif
}

void fft_execute_inverse(FFT* const fft)
{
#ifdef DOUBLEPRECISION
  fftw_execute_dft_c2r(fft->inverse_plan, fft->fk, fft->fx);
#else
  fftwf_execute_dft_c2r(fft->inverse_plan, fft->fk, fft->fx);
#endif
}

#endif

// No change whehter with or without MPI

void fft_free(FFT* const fft)
{
#ifdef DOUBLEPRECISION
  fftw_destroy_plan(fft->forward_plan);
  fftw_destroy_plan(fft->inverse_plan);
#else
  fftwf_destroy_plan(fft->forward_plan);
  fftwf_destroy_plan(fft->inverse_plan);
#endif
  
  if(fft->allocated == true) {
    free(fft->fx);
    fft->fx= NULL; fft->fk= NULL;
  }
}

// Quote
// "it is probably better for you to simply create multiple plans
//  (creating a new plan is quick once one exists for a given size)
// -- FFTW3 Manual for version 3.3.3 section 4.6

// "To prevent memory leaks, you must still call fftw_destroy_plan
//  before executing fftw_cleanup."
