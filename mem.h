#ifndef MEM_H
#define MEM_H 1

typedef struct {
  char* name;
  void* buf;
  size_t size;
} Mem;

Mem* mem_init(const char name[]);
void mem_reserve(Mem* const mem, const size_t size, const char msg[]);
void mem_alloc_reserved(Mem* const mem);

Mem* mem_alloc(const char name[], const size_t size);
  
void mem_free(Mem* const mem);

#endif
