#include "rigid_body.h"

using namespace std;
using namespace Eigen;

rigid_body::rigid_body()
{
  ;
}

rigid_body::~rigid_body()
{
  ;
}

rigid_body::rigid_body(std::string type,myvector location_cm,double length_side,double dmetric,double density,double gravity)
{
  // initial state 
  this->location_cm=this->location_cm_ori=location_cm; 
  this->velocity_cm=myvector(0,0,0); //with no linear velocity
  this->torque=myvector(0,0,0); 
  this->angular_velocity=myvector(0,0,0);
  this->force_external=myvector(0,0,0);
  this->length_side=length_side;
  this->gravity=gravity;
  {
    Matrix3d A=MatrixXd::Identity(3,3);
    this->q=Quaterniond(A);
    q.normalize();
  }
  if(type=="cube")
    {
      //length_side 表示边长
      this->mass=pow(length_side,3)*density;
      getMesh(type,dmetric);
    }
  else if(type=="sphere")
    {
      //length_side 表示半径
      this->mass=4.0/3.0*pi*pow(length_side,3)*density;
      getMesh(type,dmetric);
    }

  getBodySpaceInertiaMatrix();
  cal_area();
}

int rigid_body::getBodySpaceInertiaMatrix()
{
  size_t i=0;
  inertia_tensor_local.fill(0);

  double mass_average=this->mass/num_vertex;
  for(i=0;i<num_vertex;i++)
    {
      myvector dis=myvertexs[i].location-location_cm;
      inertia_tensor_local(0,0)+=(dis(1)*dis(1)+dis(2)*dis(2));
      inertia_tensor_local(1,1)+=(dis(0)*dis(0)+dis(2)*dis(2));
      inertia_tensor_local(2,2)+=(dis(0)*dis(0)+dis(1)*dis(1));

      inertia_tensor_local(0,1)-=(dis(0)*dis(1)); 
      inertia_tensor_local(0,2)-=(dis(0)*dis(2));
      inertia_tensor_local(1,2)-=(dis(1)*dis(2));
    }
  inertia_tensor_local(1,0)=inertia_tensor_local(0,1);
  inertia_tensor_local(2,0)=inertia_tensor_local(0,2);
  inertia_tensor_local(2,1)=inertia_tensor_local(1,2);

  inertia_tensor_local*=mass_average;  
  return 0;
}

int rigid_body::cal_velocity_angularVelocity_location_orientation(double dt)
{
  Matrix3d R=this->q.matrix();
  MatrixXd IM=MatrixXd::Identity(3,3);
  MatrixXd Mass_inverse=(1.0/this->mass)*IM;

  //check the inertia_tensor_global 's computing method , the equation below is wrong 
  //  Matrix3d Inertia_Tensor_global=/*R**/inertia_tensor_local/**R.transpose()*/;

  Matrix3d Inertia_Tensor_global=R*inertia_tensor_local*R.transpose();
  
  Matrix3d Inertia_Tensor_global_inserse=Inertia_Tensor_global.inverse();
  Vector3d velocity_cm_v=velocity_cm.get_v(); Vector3d force_external_v=force_external.get_v();
  Vector3d velocity_cm_v_new=Mass_inverse*force_external_v*dt+velocity_cm_v;

  Vector3d angular_velocity_v=angular_velocity.get_v(); Vector3d torque_v=torque.get_v();
  Vector3d help_3d=torque_v-angular_velocity_v.cross(Inertia_Tensor_global*angular_velocity_v);
  Vector3d angular_velocity_v_new=Inertia_Tensor_global_inserse*help_3d*dt+angular_velocity_v;

  
  /*
    myvector location_cm; //location of the center of mass
    Eigen:: Quaterniond q; //初始时的对应的rotation matrix 是单位阵

    myvector velocity_cm; //velocity of the center of mass
    myvector angular_velocity;

    myvector torque;
    myvector force_external;
  */

  Matrix<double,4,3> Q;
  {
    Q(0,0)=-1*q.x(); Q(0,1)=-1*q.y(); Q(0,2)=-1*q.z();
    Q(1,0)=q.w(); Q(1,1)=q.z(); Q(1,2)=-1*q.y();
    Q(2,0)=-1*q.z(); Q(2,1)=q.w(); Q(2,2)=q.x();
    Q(3,0)=q.y(); Q(3,1)=-1*q.x(); Q(3,2)=q.w();
  }

  //cout<<"q:"<<q.w()<<" "<<q.x()<<" "<<q.y()<<" "<<q.z()<<endl;
  // cout<<"Q:"<<endl<<Q<<endl;
  
  Vector4d help_4d=Q*angular_velocity_v;
  Quaterniond dq=Quaterniond(help_4d(0)*dt,help_4d(1)*dt,help_4d(2)*dt,help_4d(3)*dt);
  Quaterniond q_new=Quaterniond(q.w()+dq.w(),q.x()+dq.x(),q.y()+dq.y(),q.z()+dq.z());
  q_new.normalize();
  Vector3d location_cm_v=location_cm.get_v();
  Vector3d location_cm_v_new=velocity_cm_v*dt+location_cm_v;

  location_cm.set_v(location_cm_v_new);
  q=Quaterniond(q_new);
  velocity_cm.set_v(velocity_cm_v_new);
  angular_velocity.set_v(angular_velocity_v_new);
  
  return 0;
}

int rigid_body::cal_location_velocity_forMesh()
{
  size_t i,j;
  Matrix3d R=q.matrix();
  //cout<<"rotation matrix:"<<endl<<R;
  Vector3d location_cm_v=location_cm.get_v();
  for(i=0;i<num_vertex;i++)
    {
      Vector3d location_new_now=R*((myvertexs[i].location_ori-location_cm_ori).get_v())+location_cm_v;
      myvertexs[i].location.set_v(location_new_now);

      myvertexs[i].velocity=velocity_cm+angular_velocity.cross(myvertexs[i].location-location_cm);
    }
  return 0;
}

int rigid_body::cal_area()
{
  size_t i;
  for(i=0;i<num_tri;i++)
    {
      mytri_face[i].cal_area(location_cm,myvertexs);
    }
  return 0;
}
int rigid_body::cal_normal()
{
  size_t i;
  for(i=0;i<num_tri;i++)
    {
      mytri_face[i].cal_normal(location_cm,myvertexs);
    }
  return 0;
}


int rigid_body::cal_force_torque()
{ 
  size_t i,j;
  torque=myvector(0,0,0); force_external=myvector(0,0,0);
  //add gravity force
  force_external+=(myvector(0,0,-1*gravity)*mass);
  
  //add simple contact force test:default the ground that may contact is z=-2with normal vector is z=1
  myvector normal_vector=myvector(0,0,1);
  myvector ground_point=myvector(0,0,-120); //test the contact

  double EPS=1e-20;

  double stiff_coefficient=100000;
  for(i=0;i<num_vertex;i++)
    {
      double signal_dis=(myvertexs[i].location-ground_point).dot(normal_vector);
      if(signal_dis>=EPS)
	{
	  ; // Under this condition, no penetration happens.
	}
      else
	{
	  myvector force_here=myvector(0,0,stiff_coefficient*fabs(signal_dis));
	  force_external+=force_here;

	  printf("contact\n");
	  // do not add torque because of the contact, because it leads unstable and visual unpleasing rotation 
	  //  torque+=((myvertexs[i].location-location_cm).cross(force_here));
	}
    }

  // velocity damping
  double damping_coefficient=0;
  force_external+=(-1*damping_coefficient*velocity_cm);
  
  return 0;
}
int rigid_body::getMesh(std::string type,double dmetric)
{
  while(!myvertexs.empty())
    {
      myvertexs.pop_back();
    }
  while(!mytri_face.empty())
    {
      mytri_face.pop_back();
    }
  if(type=="cube")
    {
      ;
    }
  else if(type=="sphere")
    {
      size_t n=2;

      myvector loc_1=myvector(length_side*sin(0)*cos(0),length_side*sin(0)*sin(0),length_side*cos(0));
      while(1)
	{
	  double d_theta=pi/n;
	  double d_phi=2*pi/n;
	  myvector loc_2=myvector(length_side*sin(d_theta)*cos(d_phi),length_side*sin(d_theta)*sin(d_phi),length_side*cos(d_theta));
	 
	  if((loc_2-loc_1).len()<dmetric)
	    {
	      break;
	    }
	  else
	    {
	      n++;
	    }
	}
      if(n<30)
	{
	  n=30; //最小划分数
	}

      size_t i,j,a,b;
      size_t** index_for_vertex=(size_t**)new size_t*[n+1];
      for(i=0;i<n+1;i++)
	{
	  index_for_vertex[i]=new size_t[n+1];
	}
      double theta_now,phi_now;
     
      double EPS=1e-10;
      double dmetric_x=pi/(double)n;
      double dmetric_y=2*pi/(double)n;

      printf("n: %d\n",n);
      size_t index_now_vertex=0;

      for(j=0,phi_now=0;phi_now<2*pi-dmetric_y+EPS;phi_now+=dmetric_y,++j)
	{
	  for(i=0,theta_now=0;theta_now<pi+EPS;theta_now+=dmetric_x,++i)
	    {	 
	      index_for_vertex[i][j]=index_now_vertex;
	      myvertexs.push_back(vertex(myvector(length_side*sin(theta_now)*cos(phi_now),length_side*sin(theta_now)*sin(phi_now),length_side*cos(theta_now))+location_cm));
	      index_now_vertex++;
	    }
	}

      for(i=0,theta_now=0;theta_now<pi+EPS;theta_now+=dmetric_x,++i)
	{
	  index_for_vertex[i][n]=index_for_vertex[i][0];
	}
      
      size_t index_vertex_now[2][2];
      for(i=0;i<n;i++)
	{
	  for(j=0;j<n;j++)
	    {
	      for(a=0;a<2;a++)
		{
		  for(b=0;b<2;b++)
		    {
		      index_vertex_now[a][b]=index_for_vertex[i+a][j+b];
		    }
		}
	      mytri_face.push_back(triangle(index_vertex_now[0][0],index_vertex_now[1][1],index_vertex_now[0][1])) ;
	      mytri_face.push_back(triangle(index_vertex_now[0][0],index_vertex_now[1][0],index_vertex_now[1][1])) ;
	    }
	}

      if(index_for_vertex!=NULL)
	{
	  for(i=0;i<n+1;++i)
	    {
	      delete[] index_for_vertex[i];
	    }
	  delete[] index_for_vertex;
	}
      num_vertex=myvertexs.size();
      num_tri=mytri_face.size();
    }
  
  return 0;
}

