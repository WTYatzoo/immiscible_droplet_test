#ifndef _SWE_SOLVER_
#define _SWE_SOLVER_

#include "vertex.h"
#include "triangle.h"
#include "rigid_body.h"
class SWE_solver
{
 public:
  size_t n_x,n_y;
  size_t frame;
  double length,width;
  double dmetric;
  double dt;
  double gravity; 
  double density; // fluid density
  // visualization of the free surface of the water 
  std::vector<vertex > myvertexs; // need to save in vtk
  std::vector<triangle > mytri_face; //need  to save in vtk
  rigid_body mySphere;
  double height_max_ori,height_min_ori;
  size_t num_vertex;
  size_t num_tri;
  
  double** height;
  double*** displaced_height;
  int sign;
  size_t** num_computed_vertex; // for computing the displaced height, we use an average of the displaced_height of the center of the triangle in this cell  
  double** v_x;
  double** v_y;
  double** dhdx;
  double** dhdy;
  //存放当前advected value，保证不污染原始的数据 
  double** height_temp;
  double** v_x_temp;
  double** v_y_temp;
  size_t** index_for_vertex;
  
  std::string out_dir;
  SWE_solver();
  ~SWE_solver();
  SWE_solver(boost::property_tree::ptree &para_tree);
  int init();
  int solve();
  int cal_loc();
  int cal_velocity();
  int advect_height();
  int advect_v_x();
  int advect_v_y();
  int set_temp_to_ori();
  int update_height();
  int update_velocity();
  int boundary_condition_process();
  int triangulate();

  int coupling_F_to_R();
  int coupling_R_to_F();
};

#endif

