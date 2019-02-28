#ifndef _TRIANGLE_
#define _TRIANGLE_
#include "myvector.h"
#include "vertex.h"

class triangle
{
 public:
  size_t index_vertex[3];
  double area;
  myvector normal;
  triangle(){}
  ~triangle(){}
  triangle(const size_t (&index_vertex)[3]);
  triangle(const size_t index_vertex_a,const size_t index_vertex_b,const size_t index_vertex_c);

  int cal_area(myvector &location_cm,std::vector<vertex > &myvertexs);
  int cal_normal(myvector &location_cm,std::vector<vertex > &myvertexs);
};
#endif

