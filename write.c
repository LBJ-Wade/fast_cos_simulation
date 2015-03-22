#include <stdio.h>
#include <assert.h>
#include "particle.h"

void write_particles_txt(char filename[], Particles* particles)
{
  FILE* fp= fopen(filename, "w"); assert(fp);

  Particle const * const p= particles->p;
  for(size_t i = 0; i < particles->np_local; i++) {
    float_t x[3];
    const float_t boxsize= particles->boxsize;

    for(int k=0; k<3; k++) {
      x[k]= p[i].x[k];
      if(x[k] < 0) x[k] += boxsize;
      if(x[k] >= boxsize) x[k] -= boxsize;
    }

    fprintf(fp, "%e %e %e %e %e %e\n",
	    x[0], x[1], x[2],
	    p[i].v[0], p[i].v[1], p[i].v[2]);
  }
  
  fclose(fp);
  printf("%s written\n", filename);
}

