#ifndef _GRID_
#define _GRID_

#include "myvector.h"
class grid  
{
 public:
  size_t index_vertex[2][2]; 
  grid();
  grid(const size_t (&index_vertex)[2][2]);
  ~grid();
};

#endif
