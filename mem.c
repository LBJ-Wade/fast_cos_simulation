#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <fftw3-mpi.h>
#include "msg.h"
#include "util.h"
#include "mem.h"

Mem* mem_init(const char name[])
{
  // return Mem* with zero memory
  Mem* mem= malloc(sizeof(Mem)); assert(mem);
  mem->size= 0;
  mem->buf= 0;

  mem->name= util_new_str(name);

  return mem;
}

void mem_free(Mem* const mem)
{
  if(mem->buf) free(mem);
  mem->buf= 0;
}

void mem_reserve(Mem* const mem, const size_t size, char const * const msg)
{
  if(size > mem->size)
    mem->size= size;

  if(msg)
    msg_printf(info, "%s requested %lu MB for %s\n",
	       msg, size/(1024*1024), mem->name);
}

void mem_alloc_reserved(Mem* const mem)
{
  if(mem->buf) free(mem->buf);

#ifdef DOUBLEPRECISION
  mem->buf= fftw_malloc(mem->size);
#else
  mem->buf= fftwf_malloc(mem->size);
#endif
  
  if(mem->buf == 0)
    msg_abort("Error: Unable to allocate %lu MB for %s\n",
	      mem->size/(1024*1024), mem->name);
  else
    msg_printf(info, "%lu MB allocated for mem %s\n",
	       mem->size/(1024*1024), mem->name);
}

Mem* mem_alloc(const char name[], const size_t size)
{
  // return Mem* with 'size' bytes
  Mem* mem= mem_init(name);
  mem_reserve(mem, size, NULL);
  mem_alloc_reserved(mem);

  return mem;
}
