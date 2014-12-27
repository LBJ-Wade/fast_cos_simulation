#ifndef POWER_H
#define POWER_H 1

typedef struct {
  int n;
  double* log_k;
  double* log_P;
  gsl_interp *interp;
  gsl_interp_accel *acc;
  
} PowerSpectrum;

PowerSpectrum* power_init(const char filename[]);
void power_free();

double power_spectrum(PowerSpectrum* const ps, const double k);

#endif
