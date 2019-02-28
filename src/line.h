#ifndef _LINE_
#define _LINE_

#include "head.h"

class line
{
 public:
  size_t index_vertex[2];
  line(){}
  ~line(){}
  line(const size_t (&index_vertex)[2]);
  line(const size_t index_vertex_a,const size_t index_vertex_b);
};
#endif
