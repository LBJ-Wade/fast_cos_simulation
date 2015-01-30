///
/// \file  config.h
/// \brief definitions of basic variables
///

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

// Memory alignment for SIDM instructions, see
// Section 3.1 SIMD alignment and fftw_malloc in FFTW3 manulal, and
// FFTW3 kernel/align.c source code
#define ALGN 16

size_t size_align(size_t size);

#endif
