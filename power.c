#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <gsl/gsl_spline.h>
#include "comm.h"
#include "msg.h"
#include "power.h"

static void read_power_spectrum_file(const char filename[], PowerSpectrum* const ps);

PowerSpectrum* power_init(const char filename[])
{
  PowerSpectrum* const ps= malloc(sizeof(PowerSpectrum)); assert(ps);

  if(comm_this_node() == 0)
    read_power_spectrum_file(filename, ps);

  comm_bcast_int(&ps->n, 1);

  if(comm_this_node() != 0) {
    ps->log_k= malloc(sizeof(double)*2*ps->n);
    ps->log_P= ps->log_k + ps->n;
  }

  comm_bcast_double(ps->log_k, 2*ps->n);

  
  ps->interp= gsl_interp_alloc(gsl_interp_cspline, ps->n);
  ps->acc= gsl_interp_accel_alloc();

  const int n_required= (int) gsl_interp_min_size(ps->interp);
  if(ps->n < n_required)
    msg_abort("Error: Not enough power spectrum data points for cubic spline; %d data points < %d required\n", ps->n, n_required);
    
  return ps;
}

void power_free(PowerSpectrum* ps)
{
  gsl_interp_accel_free(ps->acc);
  gsl_interp_free(ps->interp);
}

double power_spectrum(PowerSpectrum* const ps, const double k)
{
  double log_P= gsl_interp_eval(ps->interp, ps->log_k, ps->log_P, k, ps->acc);
  return exp(log_P);
}
  

 void read_power_spectrum_file(const char filename[], PowerSpectrum* const ps)
{
  FILE* fp= fopen(filename, "r");
  if(fp == 0)
    msg_abort("Error: Unable to read input power spectrum: %s\n", filename);

  int nalloc= 1000;
    // Initial guess of number of data lines. Best if it's the number of the
    // power spectrum number of data entries, but works for any number.
  double* buf= malloc(sizeof(double)*2*nalloc);

  char line[128];
  int nlines= 0;
  double k, P;
  double k_prev= 0;

  while(fgets(line, 127, fp)) {
    if(nlines == nalloc) {
      nalloc *= 2;
      buf= realloc(buf, sizeof(double)*2*nalloc); assert(buf);
    }
    
    if(line[0] == '#')
      continue;
    else if(sscanf(line, "%lg %lg", &k, &P) == 2) {
      if(k < k_prev)
	msg_abort("Error: wavenumber k in the power spectrum file must be sorted in increasing order. %dth data k=%e > previous k= %e\n", nlines, k_prev, k);
      k_prev= k;
      
      buf[2*nlines    ]= k;
      buf[2*nlines + 1]= P;
      nlines++;
    }
    else {
      msg_printf(warn, "Warning: Unable to understand a line in the power spectrum file; following data are ignored: %s", line);
      break;
    }
  }

  int ret= fclose(fp); assert(ret == 0);
  
  msg_printf(verbose, "Found %d pairs of values in the input power spectrum\n", nlines);

  double* const v_logk= malloc(2*nlines*sizeof(double)); assert(v_logk);
  double* const v_logP= v_logk + nlines;

  for(int j=0; j<nlines; j++) {
    v_logk[j]= log(buf[2*j    ]);
    v_logP[j]= log(buf[2*j + 1]);
  }
  free(buf);
  
  ps->log_k= v_logk;
  ps->log_P= v_logP;
  ps->n= nlines;
}
