///
/// \file  msg.c
/// \brief Prints message to standard output
///

#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include "comm.h"
#include "msg.h"

static enum LogLevel log_level;
static char prefix[8]= "";

void msg_set_loglevel(const enum LogLevel lv)
{
  log_level= lv;
}

void msg_set_prefix(const char prefix_[])
{
  // Start messages with given prefix such as '#'
  strncpy(prefix, prefix_, 7);
}


void msg_printf(const enum LogLevel msg_level, char const * const fmt, ...)
{
  if(comm_this_node() == 0 && msg_level >= log_level) {
    va_list argp;

    va_start(argp, fmt);

    fprintf(stdout, "%s", prefix);
    vfprintf(stdout, fmt, argp);
    fflush(stdout);

    va_end(argp);
  }
}

void msg_abort(char const * const fmt, ...)
{
  va_list argp;

  if(log_level <= msg_fatal) {
    va_start(argp, fmt);
    vfprintf(stdout, fmt, argp);
    va_end(argp);
  }

  comm_abort();
}

void msg_assert_double(const char file[], const unsigned int line, const double x, const double x_expected, const double eps)
{
  double rel_error= fabs(x-x_expected)/x_expected;
  if(rel_error >= eps) {
    msg_abort("Assesion error at %s %u: %e != %e, %e error > %e required\n",
	  file, line, x, x_expected, rel_error, eps);
  }
}
