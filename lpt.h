#ifndef LPT_H
#define LPT_H 1

void lpt_init(const int nc, const double boxsize, Mem* mem);
void lpt_set_displacements(const unsigned long seed, PowerSpectrum* const ps,
			   const double a, Particles* particles);
  
#endif
