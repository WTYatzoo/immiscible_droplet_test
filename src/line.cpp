#include "line.h"

line::line(const size_t (&index_vertex)[2])
{
  size_t i;
  for(i=0;i<2;++i)
    {
      this->index_vertex[i]=index_vertex[i];
    }
}

line::line(const size_t index_vertex_a,const size_t index_vertex_b)
{
  index_vertex[0]=index_vertex_a;
  index_vertex[1]=index_vertex_b;
}
