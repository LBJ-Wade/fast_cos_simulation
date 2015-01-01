#ifndef CONFIG_H
#define CONFIG_H 1

#include <math.h>

#ifndef M_PI
#define M_PI        3.14159265358979323846264338327950288
#endif

#ifdef DOUBLEPRECISION
typedef double float_t;
#else
typedef float float_t;
#endif

#endif
