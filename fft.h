#ifndef FFT_H
#define FFT_H 1

#include <stdbool.h>

typedef struct {
  char*       name;
  int         nc;
  float_t*    fx;
  complex_t*  fk;
  ptrdiff_t   local_nx, local_ix0;
  fftwf_plan  forward_plan, inverse_plan;
  ptrdiff_t   ncomplex;
  bool        allocated;
} FFT;

FFT* fft_alloc(const char name[], const int nc, Mem* mem, unsigned flags);
void fft_execute_forward(FFT* const fft);
void fft_execute_inverse(FFT* const fft);
void fft_free(FFT* const fft);
void fft_finalize(void);

#endif
