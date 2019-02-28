#ifndef _RIGID_BODY_
#define _RIGID_BODY_

#include "head.h"
#include "myvector.h"
#include "vertex.h"
#include "triangle.h"

class rigid_body
{
 public:
  myvector location_cm; //location of the center of mass
  myvector location_cm_ori; //oringinal location of the center of the mass
  Eigen:: Quaterniond q; //初始时的对应的rotation matrix 是单位阵

  myvector velocity_cm; //velocity of the center of mass
  myvector angular_velocity;

  myvector torque;
  myvector force_external;
  
  double length_side; //边长 or 半径
  double mass;
  double gravity;
  Eigen::Matrix<double,3,3> inertia_tensor_local; //in body's cordinate, this will not be changed  

  size_t num_vertex;
  size_t num_tri;
  std::vector<vertex > myvertexs;
  std::vector<triangle > mytri_face;

  rigid_body();
  ~rigid_body();
  rigid_body(std::string type,myvector location_cm,double length_side,double dmetric,double density,double gravity);//cube or sphere
  int getMesh(std::string type,double dmetric);
  int getBodySpaceInertiaMatrix();

  int cal_velocity_angularVelocity_location_orientation(double dt); // update velocity & angular velocity & location & orientation
  int cal_location_velocity_forMesh();
  int cal_force_torque();
  int cal_normal();
  int cal_area();
};

#endif  
