#include <stdlib.h>
#include "util.h"


char* util_new_str(char const * const src)
{
  if(src == NULL) {
    char* new_str= malloc(1);
    new_str[0]= '\n';
    return new_str;
  }
  
  const int n= strlen(src);
  char* new_str= malloc(n+1);
  strcpy(new_str, src);

  return new_str;
}

size_t mbytes(size_t bytes)
{
  return bytes/(1024*1024);
}

