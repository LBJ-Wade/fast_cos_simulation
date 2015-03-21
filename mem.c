///
/// \file  mem.c
/// \brief Memory (RAM) management
///

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <fftw3-mpi.h>
#include "config.h"
#include "msg.h"
#include "util.h"
#include "mem.h"
#include "fft.h"

// step 1: Mem* mem= mem_init("name");
// step 2: mem_reserve(mem, size1, "usage1");
//       : mem_reserve(mem, size2, "usage2");
// step 3: mem_alloc_reserved(mem);
//

Mem* mem_init(const char name[])
{
  // return Mem* with zero memory
  Mem* mem= malloc(sizeof(Mem)); assert(mem);
  mem->size_alloc= mem->size_using= 0;
  mem->buf= NULL;

  mem->name= util_new_str(name);

  msg_printf(msg_verbose, "Memory %s initilised.\n", name);
  
  return mem;
}

void mem_free(Mem* const mem)
{
  if(mem->buf) free(mem);
  mem->buf= 0;
}

void mem_reserve(Mem* const mem, size_t size, char const * const msg)
{
  if(size % ALGN != 0)
    size += ALGN - (size % ALGN);
  assert(size % ALGN == 0);
  
  if(size > mem->size_using)
    mem->size_using= size;   // this is the amount going to be allocated

  if(msg)
    msg_printf(msg_info, "%s requested %lu MB for %s\n",
	       msg, size/(1024*1024), mem->name);
}

void mem_alloc_reserved(Mem* const mem)
{
  if(mem->buf) free(mem->buf);

  mem->buf= fft_malloc(mem->size_using);

  if(mem->buf == 0)
    msg_abort("Error: Unable to allocate %lu MB for %s\n",
	      mem->size_using/(1024*1024), mem->name);
  else
    msg_printf(msg_info, "%lu MB allocated for mem %s\n",
	       mem->size_using/(1024*1024), mem->name);

  mem->size_alloc= mem->size_using;
  mem->size_using= 0;
}

Mem* mem_alloc(const char name[], const size_t size)
{
  // return Mem* with 'size' bytes
  Mem* mem= mem_init(name);
  mem_reserve(mem, size, NULL);
  mem_alloc_reserved(mem);

  return mem;
}

void* mem_use_from_zero(Mem* const mem, size_t size)
{
  size= size_align(size); 

  if(size > mem->size_alloc)
    msg_abort("Error: Unable to use %lu MB in Mem %s (only %lu MB allocated)\n",
	      mbytes(size), mem->name, mbytes(mem->size_alloc));

  mem->size_using= size;
  return mem->buf;
}

void* mem_use_remaining(Mem* const mem, size_t size)
{
  size= size_align(size);

  //printf("debug %lu %lu\n", mem->size_using, mem->size_alloc);
  
  if(size + mem->size_using > mem->size_alloc)
    msg_abort("Error: Unable to use %lu MB in Mem %s; %lu MB allocated, "
	      "%lu remaining.\n",
	      mbytes(size), mem->name, mbytes(mem->size_alloc),
	      mbytes(mem->size_alloc - mem->size_using));
  
  assert(mem->size_using % sizeof(float) == 0);
  float* p= mem->buf;
  size_t n= mem->size_using / sizeof(float); //printf("n= %lu\n", n);

  mem->size_using += size;
  msg_printf(msg_verbose, "Using %lu of %lu in memory %s\n",
	     mem->size_using, mem->size_alloc, mem->name);


  return p+n;
}

  
