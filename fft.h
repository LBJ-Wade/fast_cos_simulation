#ifndef FFT_H
#define FFT_H 1

#include <stdbool.h>
#include "config.h"
#include "mem.h"

typedef struct {
  char*       name;
  int         nc;
  float_t*    fx;
  complex_t*  fk;
  ptrdiff_t   local_nx, local_ix0;
  ptrdiff_t   local_nky, local_iky0;
  fftwf_plan  forward_plan, inverse_plan;
  ptrdiff_t   ncomplex;
  bool        allocated;
} FFT;

size_t fft_mem_size_working(const int nc, const int transposed);
size_t fft_mem_size_fk(const int nc, const int transposed);
size_t fft_local_nx(const int nc);
  
FFT* fft_alloc(const char name[], const int nc, Mem* mem, const int transposed);
void fft_execute_forward(FFT* const fft);
void fft_execute_inverse(FFT* const fft);
void fft_free(FFT* const fft);
void fft_finalize(void);
void* fft_malloc(size_t size);

#endif
