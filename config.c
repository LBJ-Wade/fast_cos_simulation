#include <assert.h>
#include "config.h"
#include "msg.h"


void config_assert(void)
{
  assert(ALGN % sizeof(float_t) == 0);
  assert(ALGN % sizeof(double) == 0);
}

size_t size_align(size_t size)
{
  if(size % ALGN != 0)
    size += ALGN - (size % ALGN);

  return size;
}


