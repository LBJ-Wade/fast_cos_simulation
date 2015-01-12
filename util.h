#ifndef UTIL_H
#define UTIL_H 1

#include <string.h>
#include <assert.h>
#include "particle.h"

char* util_new_str(char const * const);

size_t mbytes(size_t bytes);

static inline void periodic_wrapup_p(Particle* const p, const float_t boxsize)
{
  for(int k=0; k<3; k++) {
    while(p->x[k] < 0) p->x[k] += boxsize;
    while(p->x[k] >= boxsize) p->x[k] -= boxsize;
    assert(0 <= p->x[k] && p->x[k] < boxsize); // !!! heavy assert !!!
  }
}


#endif
