#ifndef MEM_H
#define MEM_H 1

typedef struct {
  char* name;
  void* buf;
  size_t size;
} Mem;

Mem* mem_init(const char name[]);
void mem_request(Mem* const mem, const size_t size, const char msg[]);
void mem_alloc(Mem* const mem);
void mem_free(Mem* const mem);

#endif
