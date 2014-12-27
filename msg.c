#include <stdio.h>
#include <stdarg.h>
#include <math.h>
#include <mpi.h>
#include "comm.h"
#include "msg.h"

static enum LogLevel log_level;

void msg_set_loglevel(const enum LogLevel lv)
{
  log_level= lv;
}


void msg_printf(const enum LogLevel msg_level, const char *fmt, ...)
{
  if(comm_this_node() == 0 && msg_level >= log_level) {
    va_list argp;

    va_start(argp, fmt);
    vfprintf(stdout, fmt, argp);
    fflush(stdout);
    va_end(argp);
  }
}

void msg_abort(const char *fmt, ...)
{
  va_list argp;

  va_start(argp, fmt);
  vfprintf(stderr, fmt, argp);
  va_end(argp);

  MPI_Abort(MPI_COMM_WORLD, 1);
}

void msg_assert_double(const char file[], const unsigned int line, const double x, const double x_expected, const double eps)
{
  double rel_error= fabs(x-x_expected)/x_expected;
  if(rel_error >= eps) {
    msg_abort("Assesion error at %s %u: %e != %e, %e error > %e required\n",
	  file, line, x, x_expected, rel_error, eps);
  }
}
