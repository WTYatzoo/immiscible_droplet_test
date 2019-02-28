#ifndef _DROPLET_
#define _DROPLET_

#include "head.h"
#include "myvector.h"
#include "vertex.h"
#include "line.h"

class droplet
{
 public:
  myvector location_cm;
  std::vector<vertex > myvertexs;
  std::vector<line > mylines;

  size_t num_vertex;// 一圈有多少个点
  double area;
  double surface_tension_coefficient;
  size_t* index_for_vertex; // need it! because we need some local query for norm
  droplet();
  ~droplet();
  droplet(myvector location_cm, double long_axis,double short_axis,double dmetric);

  int getLine(double dmetric,double long_axis,double short_axis);

  int cal_location_as_MCF(double dt);
  int cal_intermediate_velocity(double dt);
  int cal_intermediate_location(double dt);
  int update_location_velocity(double dt);
  int smooth_velocity_backward();
  int smooth_velocity_forward();
  int smooth_velocity_jacobi_style();
};

#endif 
