#include "myvector.h"
#include "SWE_solver.h"
#include "rigid_body.h"
#include "droplet.h"
#include "io.h"

using namespace std;

const static int offset[4][2]={
  {
    0,0
  },
  {
    -1,0
  },
  {
    -1,-1
  },
  {
    0,-1
  }
};
SWE_solver::SWE_solver()
{
  ;
}
SWE_solver::~SWE_solver()
{
  //decoupled delete 
  size_t i,j;
  if(index_for_vertex!=NULL)
    {
      for(i=0;i<n_x+1;++i)
	{
	  delete[] index_for_vertex[i];
	}
      delete[] index_for_vertex;
    }
  if(height!=NULL)
    {
      for(i=0;i<n_x+1;++i)
	{
	  delete[] height[i];
	}
      delete[] height;
    }

  if(num_computed_vertex!=NULL)
    {
      for(i=0;i<n_x+1;++i)
	{
	  delete[] num_computed_vertex[i];
	}
      delete[] num_computed_vertex;
    }
  
  if(displaced_height!=NULL)
    {
      for(i=0;i<2;++i)
	{
	  for(j=0;j<n_x+1;++j)
	    {
	      delete[] displaced_height[i][j];
	    }
	  delete[] displaced_height[i];
	}
      delete[] displaced_height;
    }
  if(v_x!=NULL)
    {
      for(i=0;i<n_x+1;++i)
	{
	  delete[] v_x[i];
	}
      delete[] v_x;
    }
  if(v_y!=NULL)
    {
      for(i=0;i<n_x+1;++i)
	{
	  delete[] v_y[i];
	}
      delete[] v_y;
    }
  if(height_temp!=NULL)
    {
      for(i=0;i<n_x+1;++i)
	{
	  delete[] height_temp[i];
	}
      delete[] height_temp;
    }
  
  if(v_x_temp!=NULL)
    {
      for(i=0;i<n_x+1;++i)
	{
	  delete[] v_x_temp[i];
	}
      delete[] v_x_temp;
    }
  if(v_y_temp!=NULL)
    {
      for(i=0;i<n_x+1;++i)
	{
	  delete[] v_y_temp[i];
	}
      delete[] v_y_temp;
    }
   if(dhdx!=NULL)
    {
      for(i=0;i<n_x+1;++i)
	{
	  delete[] dhdx[i];
	}
      delete[] dhdx;
    }
    if(dhdy!=NULL)
    {
      for(i=0;i<n_x+1;++i)
	{
	  delete[] dhdy[i];
	}
      delete[] dhdy;
    }
}

SWE_solver::SWE_solver(boost::property_tree::ptree &para_tree)
{
  out_dir=para_tree.get<string>("out_dir.value"); // 此时给到最深的一层
  dt=para_tree.get<double>("dt.value");
  dmetric=para_tree.get<double>("dmetric.value"); 
  n_x=para_tree.get<int>("n_x.value");
  n_y=para_tree.get<int>("n_y.value");
  frame=para_tree.get<int>("frame.value");
  height_max_ori=para_tree.get<double>("height_max_ori.value");
  height_min_ori=para_tree.get<double>("height_min_ori.value");
  gravity=para_tree.get<double >("gravity.value");
  density=para_tree.get<double >("fluid_density.value");

  double rigid_body_density=para_tree.get<double >("rigid_body_density.value");
  string type="sphere";
  mySphere=rigid_body(type,myvector(30,20,18),5,2*dmetric,rigid_body_density,gravity);
  init();
  solve();
}

//right
int SWE_solver::init()
{
  /*
    double** height;  // x: 0 n_x-1 y: 0 n_y-1
    double** v_x; // x: 0 n_x y: 0 n_y-1
    double** v_y; // x: 0 n_x-1 y: 0 n_y
    size_t** index_for_vertex; // x: 0 n_x y:0 n_y
  */ 
  while(!myvertexs.empty())
    {
      myvertexs.pop_back();
    }
  while(!mytri_face.empty())
    {
      mytri_face.pop_back();
    }

  this->height=NULL;
  
  this->displaced_height=NULL;
  this->num_computed_vertex=NULL;
  
  this->v_x=NULL;
  this->v_y=NULL;

  this->dhdy=NULL;
  this->dhdx=NULL;
  
  this->index_for_vertex=NULL;
  this->height_temp=NULL;
  this->v_x_temp=NULL;
  this->v_y_temp=NULL;
  this->length=dmetric*n_x; this->width=dmetric*n_y;

  this->sign=0;
  double length_now,width_now;
  size_t i,j;
  size_t x,y;
  index_for_vertex=(size_t**)new size_t*[n_x+1];
  height=(double**)new double*[n_x+1];
  v_x=(double**)new double*[n_x+1];
  v_y=(double**)new double*[n_x+1];

  //for rigid body and fluid coupling 
  dhdx=(double**)new double*[n_x+1];
  dhdy=(double**)new double*[n_x+1];

  height_temp=(double**)new double*[n_x+1];
  v_x_temp=(double**)new double*[n_x+1];
  v_y_temp=(double**)new double*[n_x+1];

  num_computed_vertex=(size_t**)new size_t*[n_x+1];
  // displaced_height buffer 0 & 1, use sign to exchange the new and old buffer 
  {
    displaced_height=(double***)new double**[2];
    for(i=0;i<2;i++)
      {
	displaced_height[i]=new double*[n_x+1];
	for(j=0;j<n_x+1;j++)
	  {
	    displaced_height[i][j]=new double[n_y+1];
	  }
      }
  }
  
  for(i=0;i<n_x+1;i++)
    {
      num_computed_vertex[i]=new size_t[n_y+1];
      
      index_for_vertex[i]=new size_t[n_y+1];
      height[i]=new double[n_y+1];
      
      v_x[i]=new double[n_y+1];
      v_y[i]=new double[n_y+1];

      dhdx[i]=new double[n_y+1];
      dhdy[i]=new double[n_y+1];

      height_temp[i]=new double[n_y+1];
      v_x_temp[i]=new double[n_y+1];
      v_y_temp[i]=new double[n_y+1];
    }

  size_t index_now_vertex=0;

  double EPS=1e-6;
  for(i=0,length_now=0;length_now<length+EPS;length_now+=dmetric,++i)
    {
      for(j=0,width_now=0;width_now<width+EPS;width_now+=dmetric,++j)
	{	
	  index_for_vertex[i][j]=index_now_vertex;
	  myvertexs.push_back(vertex(myvector(length_now,width_now,0))); //初始化时高度未知
	  index_now_vertex++;	    
	}
    }

  for(i=0;i<=n_x;i++)
    {
      for(j=0;j<=n_y;j++)
	{
	  v_x[i][j]=v_y[i][j]=0;
	  displaced_height[0][i][j]=displaced_height[1][i][j]=height[i][j]=0; //让最外层height失效
	}
    }

  for(i=0;i<n_x;i++)
    {
      for(j=0;j<n_y;j++)
	{
	  height[i][j]=height_min_ori; //single dam test
	}
    }

  size_t scale=1;
  for(i=70*scale;i<=80*scale;i++) // the base scale is the 151*151
    {
      for(j=70*scale;j<=80*scale;j++)
	{
	  height[i][j]=height_max_ori; //single dam test
	}
    }
  num_vertex=myvertexs.size();

  printf("num_vertex: %d\n",num_vertex);
  return 0;
}


int SWE_solver::solve()
{
  int i,j;
  io myio=io();
  
  // test rigid body simulation, no coupling with fluid
  // {
  //   string type="sphere";
  //   rigid_body mySphere=rigid_body(type,myvector(0,0,0),0.1,0.01,100);
  //   string path_name_sphere=out_dir+"/sphere.vtk";
  //   myio.saveTriAsVTK(mySphere.myvertexs,mySphere.mytri_face,path_name_sphere);

  //   int inner_loop_rigid_body=500;
  //   for(i=0;i<1000;i++)
  //     {
  // 	for(j=0;j<inner_loop_rigid_body;j++)
  // 	  {
  // 	    mySphere.cal_force_torque();
  // 	    mySphere.cal_velocity_angularVelocity_location_orientation(dt);
  // 	    mySphere.cal_location_forMesh();	      
  // 	    if(j==0)
  // 	      {
  // 		stringstream ss;string frame_str;ss<<i; ss>>frame_str;
  // 		string path_name=out_dir+"/sphere_"+frame_str+".vtk";
  // 		myio.saveTriAsVTK(mySphere.myvertexs,mySphere.mytri_face,path_name);
  // 		printf("[INFO]:: saveAsVTK succeed\n");
  // 	      }
  // 	  }
  //     }
  // }  

  //test droplet with 2017 hyperbolic mean curvature flow simulation 
  // {
  //   droplet myDroplet=droplet(myvector(0,0,1),0.01,0.0008,0.00004);
    
  //   string path_name_droplet=out_dir+"/droplet_bad.vtk";
  //   myio.saveLineAsVTK(myDroplet.myvertexs,myDroplet.mylines,path_name_droplet);
  //   int inner_loop_droplet=50;
  //   for(i=0;i<0;i++)
  //     {
  // 	for(j=0;j<inner_loop_droplet;j++)
  // 	  {
  // 	    // myDroplet.cal_location_as_MCF(dt);
	   
  // 	      myDroplet.cal_intermediate_velocity(dt);
  // 	      myDroplet.cal_intermediate_location(dt);
  // 	      myDroplet.update_location_velocity(dt);
	      
  // 	      myDroplet.smooth_velocity_jacobi_style();

  // 	      int smooth_step_num=0;
  // 	      while(smooth_step_num<0)
  // 		{
  // 		  ///	  myDroplet.smooth_velocity_jacobi_style();

  // 	       		  myDroplet.smooth_velocity_forward();
  // 		  	  myDroplet.smooth_velocity_backward();
  // 		  smooth_step_num++;
  // 		}
	      
  // 	    if(j==0)
  // 	    {
  // 	      stringstream ss;string frame_str;ss<<i; ss>>frame_str;
  // 	      string path_name=out_dir+"/droplet_bad_"+frame_str+".vtk";
  // 	      myio.saveLineAsVTK(myDroplet.myvertexs,myDroplet.mylines,path_name);
  // 	      printf("[INFO]:: saveAsVTK succeed\n");
  // 	    }
  // 	  }
  //     }
  // }

  //test shallow water equation simulation,  no coupling with rigid body 
  // {
  //   int inner_loop_SWE=50;
  //   for(i=0;i<1000;i++)
  //     {
  // 	for(j=0;j<inner_loop_SWE;j++)
  // 	  {
  // 	    cal_loc();
  // 	    cal_velocity();
  // 	    if(j==0)
  // 	      {
  // 		triangulate();
  // 		stringstream ss;string frame_str;ss<<i; ss>>frame_str;
  // 		string path_name=out_dir+"/free_surface_"+frame_str+".vtk";
  // 		myio.saveTriAsVTK(myvertexs,mytri_face,path_name);
  // 		printf("[INFO]:: saveAsVTK succeed\n");
  // 	      }
	    
  // 	    //for next frame's simulation
  // 	    advect_height();
  // 	    advect_v_x();
  // 	    advect_v_y();
  // 	    set_temp_to_ori();
  // 	    update_height();
  // 	    update_velocity();
  // 	    boundary_condition_process();
  // 	  }      
  //     }
  // }


  //test rigid body simulation, coupling with fluid
  {
    int inner_loop=500;
    for(i=0;i<300;i++)
      {
	for(j=0;j<inner_loop;j++)
	  {	    	    
	    //fluid part
	    cal_loc();
  	    cal_velocity();  	    	    
  	    //for next frame's simulation
  	    advect_height();
  	    advect_v_x();
  	    advect_v_y();
  	    set_temp_to_ori();
  	    update_height();
  	    update_velocity();
  	    boundary_condition_process();

	    /*
	    //rigid body part with fluid's influence
	    //    mySphere.velocity_cm(1)=4;
	    mySphere.cal_force_torque();
	    mySphere.cal_normal();
	    coupling_F_to_R();
	    mySphere.cal_velocity_angularVelocity_location_orientation(dt);
	    mySphere.cal_location_velocity_forMesh();

	    //resolve fluid part with rigid body's influence 
	    coupling_R_to_F(); 
	    if(j==0)
	      {
		stringstream ss;string frame_str;ss<<i; ss>>frame_str;
		string path_name=out_dir+"/sphere_"+frame_str+".vtk";
		myio.saveTriAsVTK(mySphere.myvertexs,mySphere.mytri_face,path_name);
		printf("[INFO]:: Rigid_Body saveAsVTK succeed\n");
		} */
	    if(j==0)
  	      {
  		triangulate();
  		stringstream ss;string frame_str;ss<<i; ss>>frame_str;
  		string path_name=out_dir+"/free_surface_"+frame_str+".vtk";
  		myio.saveTriAsVTK(myvertexs,mytri_face,path_name);
  		printf("[INFO]:: Fluid saveAsVTK succeed\n");
  	      }
	  }
      }
  }
  return 0;
}

//compute the influence of rigid body to fluid
int SWE_solver::coupling_R_to_F()
{
  size_t i,j,k;
  int sign_new=sign; int sign_old=(sign+1)%2;

  for(i=0;i<n_x+1;i++)
    {
      for(j=0;j<n_y+1;j++)
	{
	  num_computed_vertex[i][j]=0;
	  displaced_height[sign_new][i][j]=0;
	}
    }

  size_t index_vertex_here[3];
  for(i=0;i<mySphere.num_tri;++i)
    {
      for(j=0;j<3;j++)
	{
	  index_vertex_here[j]=mySphere.mytri_face[i].index_vertex[j];
	}
      myvector location_cm_tri=myvector(0,0,0);
      for(j=0;j<3;j++)
	{
	  location_cm_tri+=mySphere.myvertexs[index_vertex_here[j]].location;
	}
      location_cm_tri/=3.0;
      size_t which_x,which_y;
      which_x=size_t((location_cm_tri(0)-0.5*dmetric)/dmetric); which_y=size_t((location_cm_tri(1)-0.5*dmetric)/dmetric);
      if(height[which_x][which_y]>location_cm_tri(2))
	{
	  displaced_height[sign_new][which_x][which_y]+=(height[which_x][which_y]-location_cm_tri(2));
	  num_computed_vertex[which_x][which_y]++;
	}
    }
  for(i=0;i<n_x+1;i++)
    {
      for(j=0;j<n_y+1;j++)
	{
	  if(num_computed_vertex[i][j]>0)
	    {
	      displaced_height[sign_new][i][j]=displaced_height[sign_new][i][j]/num_computed_vertex[i][j];
	    }
	}
    }

  double scale=0.7;
  for(i=1;i<n_x;i++)
    {
      for(j=1;j<n_y;j++)
	{
	  double height_add=0.25*scale*(displaced_height[sign_new][i][j]-displaced_height[sign_old][i][j]);
	  height[i-1][j]+=height_add; height[i+1][j]+=height_add; height[i][j-1]+=height_add; height[i][j+1]+=height_add;
	  height[i][j]-=(height_add*4);
	}
    }
  
  sign=(sign+1)%2;
  return 0;
}

//compute the influence of fluid to rigid body 
int SWE_solver::coupling_F_to_R()
{
  size_t i,j,k;
  size_t index_vertex_here[3];

  myvector normal_ground=myvector(0,0,1);
  for(i=0;i<mySphere.num_tri;++i)
    {
      for(j=0;j<3;j++)
	{
	  index_vertex_here[j]=mySphere.mytri_face[i].index_vertex[j];
	}
      myvector location_cm_tri=myvector(0,0,0);
      myvector velocity_cm_tri=myvector(0,0,0);
      for(j=0;j<3;j++)
	{
	  location_cm_tri+=mySphere.myvertexs[index_vertex_here[j]].location;
	  velocity_cm_tri+=mySphere.myvertexs[index_vertex_here[j]].velocity;
	}
      location_cm_tri/=3.0;
      velocity_cm_tri/=3.0;

      size_t which_x,which_y;
      double mod_x,mod_y;  
      which_x=size_t((location_cm_tri(0)-0.5*dmetric)/dmetric); which_y=size_t((location_cm_tri(1)-0.5*dmetric)/dmetric); 	 
      mod_x=(location_cm_tri(0)-0.5*dmetric-which_x*dmetric)/dmetric;
      mod_y=(location_cm_tri(1)-0.5*dmetric-which_y*dmetric)/dmetric;	 	  
	  
      double height_cm_tri=height[which_x][which_y]*(1-mod_x)*(1-mod_y)+height[which_x][which_y+1]*(1-mod_x)*mod_y+height[which_x+1][which_y]*(1-mod_y)*mod_x+height[which_x+1][which_y+1]*mod_x*mod_y;

      if(height_cm_tri<location_cm_tri(2))
	{
	  ;
	}
      else
	{
	  //  printf("density : %lf gravity : %lf area: %lf normal: %lf\n",density,gravity,mySphere.mytri_face[i].area,mySphere.mytri_face[i].normal(2));
	  myvector bouyancy_force=-1*density*gravity*mySphere.mytri_face[i].area*(height_cm_tri-location_cm_tri(2))*mySphere.mytri_face[i].normal(2)*normal_ground;
	  mySphere.force_external+=bouyancy_force;
	  // printf("bouyancy_force: %lf %lf %lf\n",bouyancy_force(0),bouyancy_force(1),bouyancy_force(2));
	}


      // get dh/dx from bilinear interpolation 
      which_x=int(location_cm_tri(0)/dmetric); mod_x=(location_cm_tri(0)-which_x*dmetric)/dmetric;
      which_y=int((location_cm_tri(1)-dmetric*0.5)/dmetric); mod_y=(location_cm_tri(1)-dmetric*0.5-which_y*dmetric)/dmetric;
      double v_x_now=v_x[which_x][which_y]*(1-mod_x)*(1-mod_y)+v_x[which_x+1][which_y]*(1-mod_y)*mod_x+v_x[which_x][which_y+1]*(1-mod_x)*mod_y+v_x[which_x+1][which_y+1]*mod_y*mod_x;
      double dhdx_now=dhdx[which_x][which_y]*(1-mod_x)*(1-mod_y)+dhdx[which_x+1][which_y]*(1-mod_y)*mod_x+dhdx[which_x][which_y+1]*(1-mod_x)*mod_y+dhdx[which_x+1][which_y+1]*mod_y*mod_x;

      //get dh/dy from bilinear interpolation 
      which_x=int((location_cm_tri(0)-0.5*dmetric)/dmetric); mod_x=(location_cm_tri(0)-0.5*dmetric-which_x*dmetric)/dmetric;
      which_y=int(location_cm_tri(1)/dmetric); mod_y=(location_cm_tri(1)-which_y*dmetric)/dmetric;
	  
      double v_y_now=v_y[which_x][which_y]*(1-mod_x)*(1-mod_y)+v_y[which_x+1][which_y]*(1-mod_y)*mod_x+v_y[which_x][which_y+1]*(1-mod_x)*mod_y+v_y[which_x+1][which_y+1]*mod_y*mod_x;
      double dhdy_now=dhdy[which_x][which_y]*(1-mod_x)*(1-mod_y)+dhdy[which_x+1][which_y]*(1-mod_y)*mod_x+dhdy[which_x][which_y+1]*(1-mod_x)*mod_y+dhdy[which_x+1][which_y+1]*mod_y*mod_x;

      myvector v_fluid_now=myvector(v_x_now,v_y_now,v_x_now*dhdx_now+v_y_now*dhdy_now);
      myvector v_relative=velocity_cm_tri-v_fluid_now;

      //  printf("v_relative: %lf %lf %lf\n",v_relative(0),v_relative(1),v_relative(2));
      double area_eff;
      double omega=0;
      double coe_drag=2000;
      double coe_lift=1000;
      if(height_cm_tri<location_cm_tri(2)||v_relative.dot(mySphere.mytri_face[i].normal)<0)
	{
	  ;
	}
      else
	{
	  area_eff=(mySphere.mytri_face[i].normal.dot(v_relative)/v_relative.len()*omega+(1-omega))*mySphere.mytri_face[i].area; // from 07 sig wave particles 
	  myvector f_drag=-0.5*coe_drag*area_eff*v_relative.len()*v_relative;
	  myvector n_v=mySphere.mytri_face[i].normal.cross(v_relative); n_v.normalize();
	  //	  printf("f_drag: %lf %lf %lf area_eff: %lf area: %lf\n",f_drag(0),f_drag(1),f_drag(2),area_eff,mySphere.mytri_face[i].area);
	  myvector f_lift=-0.5*coe_lift*area_eff*v_relative.len()*(v_relative.cross(n_v)); 
	  mySphere.force_external+=f_drag;
	  mySphere.force_external+=f_lift;
	}
      
    }
  return 0;
}
//right
int SWE_solver::boundary_condition_process()
{
  size_t i,j;
  
  //down
  for(i=0;i<n_x;i++)
    {
      height[i][0]=height[i][1];
    }
  for(i=0;i<n_x;i++)
    {
      v_y[i][1]=0;
    }
  for(i=0;i<=n_x;i++)
    {
      v_x[i][0]=0;
    }
  
  //left
  for(j=0;j<n_y;j++)
    {
      height[0][j]=height[1][j];
    }
  for(j=0;j<n_y;j++)
    {
      v_x[1][j]=0;
    }
  for(j=0;j<=n_y;j++)
    {
      v_y[0][j]=0;
    }

  
  //right
  for(j=0;j<n_y;j++)
    {
      height[n_x-1][j]=height[n_x-2][j];
    }
  for(j=0;j<n_y;j++)
    {
      v_x[n_x-1][j]=0;
    }
  for(j=0;j<=n_y;j++)
    {
      v_y[n_x-1][j]=0;
    }

  //up
  for(i=0;i<n_x;i++)
    {
      height[i][n_y-1]=height[i][n_y-2];
    }
  for(i=0;i<n_x;i++)
    {
      v_y[i][n_y-1]=0;
    }
  for(i=0;i<=n_x;i++)
    {
      v_x[i][n_y-1]=0;
    }
  
  return 0;
}

//right
int SWE_solver::update_height()
{
  size_t i,j;
  for(j=1;j<n_y;j++)
    {
      for(i=1;i<n_x;i++)
	{
	  height[i][j]-=0.5*(height[i][j]*((v_x[i+1][j]-v_x[i][j])/dmetric+(v_y[i][j+1]-v_y[i][j])/dmetric)*dt);
	}
    }
  return 0;
}

// right
int SWE_solver::update_velocity()
{
  size_t i,j;
  for(j=1;j<n_y;j++)
    {
      for(i=2;i<n_x;i++)
	{
	  v_x[i][j]+=gravity*(height[i-1][j]-height[i][j])*dt/dmetric;

	  //for coupling
	  dhdx[i][j]=(height[i-1][j]-height[i][j])/dmetric;
	}
    }
  for(j=2;j<n_y;j++)
    {
      for(i=1;i<n_x;i++)
	{
	  v_y[i][j]+=gravity*(height[i][j-1]-height[i][j])*dt/dmetric;

	  //for coupling
	  dhdy[i][j]=(height[i][j-1]-height[i][j])/dmetric;
	}
    }
  return 0;
}


int SWE_solver::set_temp_to_ori()
{
  size_t i,j;
  for(i=1;i<n_x-1;i++)
    {
      for(j=1;j<n_y-1;j++)
	{
       	  height[i][j]=height_temp[i][j];
	}
    }
  for(i=1;i<n_x;i++)
    {
      for(j=1;j<n_y-1;j++)
	{
	  v_x[i][j]=v_x_temp[i][j];
	}
    }
  for(i=1;i<n_x-1;i++)
    {
      for(j=1;j<n_y;j++)
	{
	  v_y[i][j]=v_y_temp[i][j];
	}
    }
  return 0;
}

int SWE_solver::advect_v_x()
{
  size_t i,j,a,b;
  
  myvector loc_now;
  myvector velocity_now;
  myvector loc_back;

  size_t index_vertex_now[2][2];
  for(i=1;i<n_x;i++)
    {
      for(j=1;j<n_y-1;j++)
	{
	  loc_now(0)=dmetric*i;  loc_now(1)=dmetric*j+dmetric*0.5; loc_now(2)=0;
	  velocity_now(0)=v_x[i][j]; velocity_now(1)=0.25*(v_y[i][j]+v_y[i-1][j]+v_y[i][j+1]+v_y[i-1][j+1]); velocity_now(2)=0;
	  loc_back=loc_now-dt*velocity_now;
	  size_t which_x,which_y;
	  double mod_x,mod_y;

	  which_x=int(loc_back(0)/dmetric); mod_x=(loc_back(0)-which_x*dmetric)/dmetric;
	  which_y=int((loc_back(1)-dmetric*0.5)/dmetric); mod_y=(loc_back(1)-dmetric*0.5-which_y*dmetric)/dmetric;
	  v_x_temp[i][j]=v_x[which_x][which_y]*(1-mod_x)*(1-mod_y)+v_x[which_x+1][which_y]*(1-mod_y)*mod_x+v_x[which_x][which_y+1]*(1-mod_x)*mod_y+v_x[which_x+1][which_y+1]*mod_y*mod_x;	  	  
	}
    }
  return 0;
}

int SWE_solver::advect_v_y()
{
  size_t i,j,a,b;

  myvector loc_now;
  myvector velocity_now;
  myvector loc_back;

  size_t index_vertex_now[2][2];
  for(i=1;i<n_x-1;i++)
    {
      for(j=1;j<n_y;j++)
	{
	  loc_now(0)=dmetric*i+dmetric*0.5;  loc_now(1)=dmetric*j; loc_now(2)=0;
	  velocity_now(1)=v_y[i][j]; velocity_now(0)=0.25*(v_x[i][j]+v_x[i+1][j]+v_x[i][j-1]+v_x[i+1][j-1]); velocity_now(2)=0;
	  loc_back=loc_now-dt*velocity_now;

	  size_t which_x,which_y;
	  double mod_x,mod_y;

	  which_x=int((loc_back(0)-0.5*dmetric)/dmetric); mod_x=(loc_back(0)-0.5*dmetric-which_x*dmetric)/dmetric;
	  which_y=int(loc_back(1)/dmetric); mod_y=(loc_back(1)-which_y*dmetric)/dmetric;
	  
	  v_y_temp[i][j]=v_y[which_x][which_y]*(1-mod_x)*(1-mod_y)+v_y[which_x+1][which_y]*(1-mod_y)*mod_x+v_y[which_x][which_y+1]*(1-mod_x)*mod_y+v_y[which_x+1][which_y+1]*mod_y*mod_x;

	  
	}
    }
  return 0;
}

int SWE_solver::advect_height()
{
  size_t i,j,a,b;

  myvector loc_now;
  myvector velocity_now;
  myvector loc_back;

  size_t index_vertex_now[2][2];
  for(i=1;i<n_x-1;i++)
    {
      for(j=1;j<n_y-1;j++)
	{
	  loc_now(0)=dmetric*i+dmetric*0.5;  loc_now(1)=dmetric*j+dmetric*0.5; loc_now(2)=0;
	  velocity_now(0)=0.5*(v_x[i][j]+v_x[i+1][j]); velocity_now(1)=0.5*(v_y[i][j]+v_y[i][j+1]); velocity_now(2)=0;

	  size_t which_x,which_y;
	  double mod_x,mod_y;	  
	  
	  loc_back=loc_now-dt*velocity_now;	  
	  which_x=size_t((loc_back(0)-0.5*dmetric)/dmetric); which_y=size_t((loc_back(1)-0.5*dmetric)/dmetric); 	 
	  mod_x=(loc_back(0)-0.5*dmetric-which_x*dmetric)/dmetric;
	  mod_y=(loc_back(1)-0.5*dmetric-which_y*dmetric)/dmetric;	 	  
	  
	  height_temp[i][j]=height[which_x][which_y]*(1-mod_x)*(1-mod_y)+height[which_x][which_y+1]*(1-mod_x)*mod_y+height[which_x+1][which_y]*(1-mod_y)*mod_x+height[which_x+1][which_y+1]*mod_x*mod_y;
	    	        	  
	}
    }
  return 0;
}

//right
int SWE_solver::cal_loc()
{
  size_t i,j,k;
  int i_now,j_now;
  size_t num_related;
  double height_sum_related;

  double EPS=1e-6;
  for(i=0;i<=n_x;i++)
    {
      for(j=0;j<=n_y;j++)
	{
	  num_related=0; height_sum_related=0;
	  for(k=0;k<4;k++)
	    {
	      i_now=i+offset[k][0]; j_now=j+offset[k][1];
	      if(i_now>=0&&i_now<n_x&&j_now>=0&&j_now<n_y&&fabs(height[i_now][j_now])>=EPS)
		{
		  num_related++;
		  height_sum_related+=height[i_now][j_now];
		}
	    }
	  int index_now_vertex=index_for_vertex[i][j];
	  myvertexs[index_now_vertex].location.z=height_sum_related/num_related;
	}
    }
  return 0;
}

int SWE_solver::cal_velocity()
{
  size_t i,j,k;
  int i_now,j_now;
  size_t num_v_x_related;
  double v_x_sum_related;
  
  size_t num_v_y_related;
  double v_y_sum_related;

  double EPS=1e-6;
  for(i=0;i<=n_x;i++)
    {
      for(j=0;j<=n_y;j++)
	{
	  num_v_x_related=0; v_x_sum_related=0;
	  num_v_y_related=0; v_y_sum_related=0;
	  
	  {
	    if(i-1>=0&&i-1<n_x)
	      {
		num_v_y_related++;
		v_y_sum_related+=v_y[i-1][j];
	      }
	    if(i>=0&&i<n_x)
	      {
		num_v_y_related++;
		v_y_sum_related+=v_y[i][j];
	      }
	    
	  }
	  	  
	  {
	    if(j>=0&&j<n_y)
	      {
		num_v_x_related++;
		v_x_sum_related+=v_x[i][j];
	      }
	    if(j-1>=0&&j-1<n_y)
	      {
		num_v_x_related++;
		v_x_sum_related+=v_x[i][j-1];
	      }
	  }	
	  int index_now_vertex=index_for_vertex[i][j];
	  myvertexs[index_now_vertex].velocity(0)=v_x_sum_related/num_v_x_related;
	  myvertexs[index_now_vertex].velocity(1)=v_y_sum_related/num_v_y_related;
	}
    }
  return 0;
}

//right
int SWE_solver::triangulate()
{
  size_t index_vertex_now[2][2];
  size_t i,j,a,b;
  while(!mytri_face.empty())
    {
      mytri_face.pop_back();
    }
  for(i=0;i<n_x;i++)
    {
      for(j=0;j<n_y;j++)
	{
	  for(a=0;a<2;a++)
	    {
	      for(b=0;b<2;b++)
		{
		  index_vertex_now[a][b]=index_for_vertex[i+a][j+b];
		}
	    }
	  if(myvertexs[index_vertex_now[0][0]].location(2)+myvertexs[index_vertex_now[1][1]].location(2)>myvertexs[index_vertex_now[0][1]].location(2)+myvertexs[index_vertex_now[1][0]].location(2))
	    {
	      mytri_face.push_back(triangle(index_vertex_now[0][0],index_vertex_now[1][1],index_vertex_now[0][1])) ;
	      mytri_face.push_back(triangle(index_vertex_now[0][0],index_vertex_now[1][0],index_vertex_now[1][1])) ;

	    }
	  else
	    {
	      mytri_face.push_back(triangle(index_vertex_now[0][0],index_vertex_now[1][0],index_vertex_now[0][1])) ;
	      mytri_face.push_back(triangle(index_vertex_now[1][0],index_vertex_now[1][1],index_vertex_now[0][1])) ;
	    }
	}
    }
  return 0;
}

