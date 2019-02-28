#include "triangle.h"

triangle::triangle(const size_t (&index_vertex)[3])
{
  size_t i;
  for(i=0;i<3;++i)
    {
      this->index_vertex[i]=index_vertex[i];
    }
}
triangle::triangle(const size_t index_vertex_a,const size_t index_vertex_b,const size_t index_vertex_c)
 {
   index_vertex[0]=index_vertex_a;
   index_vertex[1]=index_vertex_b;
   index_vertex[2]=index_vertex_c;
 }

int triangle::cal_area(myvector &location_cm,std::vector<vertex > &myvertexs)
{
  myvector x=myvertexs[index_vertex[1]].location-myvertexs[index_vertex[0]].location;
  myvector y=myvertexs[index_vertex[2]].location-myvertexs[index_vertex[0]].location;
  myvector help=x.cross(y);
  this->area=fabs(help.len()*0.5);
  return 0;
}
int triangle::cal_normal(myvector &location_cm,std::vector<vertex > &myvertexs)
{
  myvector x=myvertexs[index_vertex[1]].location-myvertexs[index_vertex[0]].location;
  myvector y=myvertexs[index_vertex[2]].location-myvertexs[index_vertex[0]].location;
  myvector help=x.cross(y);

  myvector test=myvertexs[index_vertex[0]].location-location_cm;
  double mark=help.dot(test);

  double EPS=1e-6;
  help.normalize();
  if(mark>=EPS)
    {      
      this->normal=help;
    }
  else
    {
      this->normal=-1*help;
    }
  return 0;
}
