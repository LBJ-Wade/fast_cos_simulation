#ifndef CONFIG_H
#define CONFIG_H 1

#include <math.h>
#include <fftw3-mpi.h>

#ifdef DOUBLEPRECISION
typedef double float_t;
typedef fftw_complex complex_t;
#else
typedef fftwf_complex complex_t;
typedef float float_t;
#endif

#ifndef M_PI
#define M_PI        3.14159265358979323846264338327950288
#endif

// Memory allignment for SIMD instruction (see fftw3 source kernel/align.c)
#define ALGN 16

#endif
