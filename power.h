#ifndef POWER_H
#define POWER_H 1

typedef struct {
  int n;
  double* log_k;
  double* log_P;           // logP= log_k + n
  gsl_interp *interp;
  gsl_interp_accel *acc;
} PowerSpectrum;

PowerSpectrum* power_init(const char filename[], const double sigma8_check);
void power_free();

double power_spectrum(PowerSpectrum* const ps, const double k);

#endif
