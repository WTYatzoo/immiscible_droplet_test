#ifndef _VERTEX_
#define _VERTEX_

#include "head.h"
#include "myvector.h"

class vertex 
{
 public:
  myvector location;
  myvector velocity; //for shallow water surface's vertex, we only capture the x,y component, let z component be zero

  myvector location_ori;
  
  myvector velocity_smooth_save;

  myvector location_temp;
  myvector velocity_temp;
  vertex();
  vertex(const myvector &location);
  ~vertex();
};

#endif
