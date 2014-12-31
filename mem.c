#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "msg.h"
#include "mem.h"

Mem* mem_init(const char name_[])
{ 
  Mem* mem= malloc(sizeof(Mem)); assert(mem);
  mem->size= 0;
  mem->buf= 0;

  const int n= strlen(name_);
  mem->name= malloc(n+1);
  strcpy(mem->name, name_);

  return mem;
}

void mem_free(Mem* const mem)
{
  if(mem->buf) free(mem);
  mem->buf= 0;
}

void mem_request(Mem* const mem, const size_t size, const char msg[])
{
  if(size > mem->size)
    mem->size= size;

  if(msg)
    msg_printf(normal, "%s requested %lu MB for %s\n",
	       msg, size/(1024*1024), mem->name);
}

void mem_alloc(Mem* const mem)
{
  if(mem->buf) free(mem->buf);

  mem->buf= malloc(mem->size);

  if(mem->buf == 0) {
    msg_abort("Error: Unable to allocate %lu MB for %s\n",
	      mem->size/(1024*1024), mem->name);
  }
}
