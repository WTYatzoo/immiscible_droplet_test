#include "droplet.h"

using namespace Eigen;
using namespace std;

droplet::droplet()
{
  ;
}

droplet::~droplet()
{
  delete[] index_for_vertex;
}

droplet::droplet(myvector location_cm, double long_axis,double short_axis,double dmetric)
{
  this->surface_tension_coefficient=0.1;
  this->location_cm=location_cm;
  getLine(dmetric,long_axis,short_axis);
  //myvertexs[0].velocity=myvector(0.01,0,0);
}

int droplet::getLine(double dmetric,double long_axis,double short_axis)
{
  while(!myvertexs.empty())
    {
      myvertexs.pop_back();
    }
  while(!mylines.empty())
    {
      mylines.pop_back();
    }
  
  {
    size_t n=2;
    myvector loc_1=myvector(long_axis*cos(0),short_axis*sin(0),0);
    while(1)
      {
	double d_phi=2*pi/n;
	myvector loc_2=myvector(long_axis*cos(d_phi),short_axis*sin(d_phi),0);

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
	n=30;
      }

    this->num_vertex=n; //!!!!!!
    size_t i,a ;
    this->index_for_vertex=(size_t*)new size_t[n+1];
    double phi_now;

    double EPS=1e-10;
    double dmetric_x=2*pi/(double)n;

    printf("n: %d\n",n);
    size_t index_now_vertex=0;

    for(i=0,phi_now=0;phi_now<2*pi-dmetric_x+EPS;phi_now+=dmetric_x,++i)
      {
	index_for_vertex[i]=index_now_vertex;
	myvertexs.push_back(vertex(myvector(long_axis*cos(phi_now),short_axis*sin(phi_now),0)+location_cm));
	index_now_vertex++;
      }
    index_for_vertex[n]=index_for_vertex[0];

    size_t index_vertex_now[2];

    this->area=0;
    for(i=0;i<n;i++)
      {
	for(a=0;a<2;a++)
	  {
	    index_vertex_now[a]=index_for_vertex[i+a];
	  }
	area+=myvertexs[index_vertex_now[0]].location(0)*myvertexs[index_vertex_now[1]].location(1)-myvertexs[index_vertex_now[0]].location(1)*myvertexs[index_vertex_now[1]].location(0);
	mylines.push_back(line(index_vertex_now));
      }
    area=fabs(0.5*area);
    
  }
  return 0;
}

int droplet::cal_location_as_MCF(double dt)
{
  size_t index_adjoin[2];
  double len_adjoin;
  myvector help;
  size_t i,j;

  MatrixXd L=MatrixXd::Random(num_vertex*2,num_vertex*2); 
  VectorXd loc_old(num_vertex*2);
  VectorXd loc_new(num_vertex*2);
  L.fill(0); 
  for(i=0;i<num_vertex;i++)
    {
      loc_old(i*2)=myvertexs[i].location(0); loc_old(i*2+1)=myvertexs[i].location(1);
      index_adjoin[0]=(i+1)%num_vertex;
      index_adjoin[1]=((i-1)+num_vertex)%num_vertex;
      len_adjoin=0;
      help=myvector(0,0,0);
      for(j=0;j<2;j++)
	{
	  len_adjoin+=(myvertexs[index_adjoin[j]].location-myvertexs[i].location).len();

	  help+=(myvertexs[index_adjoin[j]].location-myvertexs[i].location); //方向可能向内也可能向外!!!!!
	}

      len_adjoin*=0.5; double len_adjoin_inverse=1.0/len_adjoin;
      L(i*2,index_adjoin[0]*2)+=len_adjoin_inverse; L(i*2,i*2)+=-2*len_adjoin_inverse; L(i*2,index_adjoin[1]*2)+=len_adjoin_inverse;

      L(i*2+1,index_adjoin[0]*2+1)+=len_adjoin_inverse; L(i*2+1,i*2+1)+=-2*len_adjoin_inverse; L(i*2+1,index_adjoin[1]*2+1)+=len_adjoin_inverse;												       
    }

  MatrixXd A_dense=MatrixXd::Identity(num_vertex*2,num_vertex*2)-dt*2*surface_tension_coefficient*L;
  SparseMatrix<double > A(2*num_vertex,2*num_vertex);
  vector< Triplet<double > > tripletsForA;

  //共轭梯度法求解线性方程组从原先的适用于正定matrix拓展到适用于ambitrary matrix
  ConjugateGradient<SparseMatrix<double> > cg;
  cg.setMaxIterations(150);

  double EPS=1e-10; 

  for(i=0;i<num_vertex*2;++i)
    {
      for(j=0;j<num_vertex*2;++j)
	{
	  if(fabs(A_dense(i,j))>= EPS)
	    {
	      tripletsForA.emplace_back(i,j,A_dense(i,j));
	    }
	}
    }
  
  A.setFromTriplets(tripletsForA.begin(),tripletsForA.end());
  A.makeCompressed();

  {
    cg.compute(A);
    loc_new=cg.solve(loc_old);
  }

  for(i=0;i<num_vertex;i++)
    {
      myvector location_new=myvector(loc_new(2*i),loc_new(2*i+1),0);
      myvertexs[i].velocity_temp=(location_new-myvertexs[i].location)/dt;
      myvertexs[i].location_temp=location_new;

      // myvertexs[i].velocity=myvertexs[i].velocity_temp;
      // myvertexs[i].location=myvertexs[i].location_temp;
    }
  return 0;
}

int droplet::cal_intermediate_velocity(double dt)
{
  size_t index_adjoin[2];
  double len_adjoin;
  myvector help;
  size_t i,j;
  for(i=0;i<num_vertex;i++)
    {
      index_adjoin[0]=(i+1)%num_vertex;
      index_adjoin[1]=((i-1)+num_vertex)%num_vertex;
      len_adjoin=0;
      help=myvector(0,0,0);
      for(j=0;j<2;j++)
	{
	  len_adjoin+=(myvertexs[index_adjoin[j]].location-myvertexs[i].location).len();
	  help+=(myvertexs[index_adjoin[j]].location-myvertexs[i].location); //方向可能向内也可能向外!!!!!
	}
      len_adjoin*=0.5;
      
      myvertexs[i].velocity_temp=myvertexs[i].velocity+dt*2*surface_tension_coefficient*help/len_adjoin;
    }
  return 0;
}

int droplet::cal_intermediate_location(double dt)
{
  size_t i;
  for(i=0;i<num_vertex;i++)
    {
      myvertexs[i].location_temp=myvertexs[i].location+dt*myvertexs[i].velocity_temp;
    }
  return 0;
}

int droplet::smooth_velocity_jacobi_style()
{
  size_t index_adjoin[2];
  double len_adjoin;
  myvector help;

  size_t i,j;
  for(i=0;i<num_vertex;i++) 
    {
      index_adjoin[0]=((i+1)+num_vertex)%num_vertex;
      index_adjoin[1]=((i-1)+num_vertex)%num_vertex;
      len_adjoin=0;
      help=myvector(0,0,0);
      for(j=0;j<2;j++)
	{
	  len_adjoin+=(myvertexs[index_adjoin[j]].velocity-myvertexs[i].velocity).len();
	  help+=(myvertexs[index_adjoin[j]].velocity-myvertexs[i].velocity); //方向可能向内也可能向外!!!!!
	}
      len_adjoin*=0.5;
      myvertexs[i].velocity_temp=(1-1e-4)*(myvertexs[i].velocity+1e-6*help/len_adjoin);

      //   myvertexs[i].velocity_temp=0.5/3*myvertexs[i].velocity_temp+2.5/3*myvertexs[i].velocity;
    }

  for(i=0;i<num_vertex;i++)
    {
      myvertexs[i].velocity=myvertexs[i].velocity_temp;
    }
  return 0;
}
int droplet::smooth_velocity_backward() //gauss-seidel style explicit mean curvature flow smooth from backward
{
  size_t index_adjoin[2];
  double len_adjoin;
  myvector help;

  int i,j;
  for(i=num_vertex-1;i>=0;i--) 
    {
      index_adjoin[0]=((i+1)+num_vertex)%num_vertex;
      index_adjoin[1]=((i-1)+num_vertex)%num_vertex;
      len_adjoin=0;
      help=myvector(0,0,0);
      for(j=0;j<2;j++)
	{
	  len_adjoin+=(myvertexs[index_adjoin[j]].velocity-myvertexs[i].velocity).len();
	  help+=(myvertexs[index_adjoin[j]].velocity-myvertexs[i].velocity); //方向可能向内也可能向外!!!!!
	}
      len_adjoin*=0.5;
      myvertexs[i].velocity=myvertexs[i].velocity+1e-6*help/len_adjoin;
    }
  return 0;
}


int droplet::smooth_velocity_forward() //gauss-seidel style explicit mean curvature flow smooth from forward
{
  size_t index_adjoin[2];
  double len_adjoin;
  myvector help;
  size_t i,j;
  for(i=0;i<num_vertex;i++)
    {
      index_adjoin[0]=(i+1)%num_vertex;
      index_adjoin[1]=((i-1)+num_vertex)%num_vertex;
      len_adjoin=0;
      help=myvector(0,0,0);
      for(j=0;j<2;j++)
	{
	  len_adjoin+=(myvertexs[index_adjoin[j]].velocity-myvertexs[i].velocity).len();
	  help+=(myvertexs[index_adjoin[j]].velocity-myvertexs[i].velocity); //方向可能向内也可能向外!!!!!
	}
      len_adjoin*=0.5;
      myvertexs[i].velocity=(1-1e-4)*myvertexs[i].velocity+1e-6*help/len_adjoin;
    }
  return 0;
}

int droplet::update_location_velocity(double dt)
{
  size_t i,j,a;
  size_t index_vertex_now[2];
  double area_temp=0;
  double len_boundary_temp=0;
  for(i=0;i<num_vertex;i++)
    {
      for(a=0;a<2;a++)
	{
	  index_vertex_now[a]=index_for_vertex[i+a];
	}
      len_boundary_temp+=(myvertexs[index_vertex_now[0]].location_temp-myvertexs[index_vertex_now[1]].location_temp).len();
      area_temp+=myvertexs[index_vertex_now[0]].location_temp(0)*myvertexs[index_vertex_now[1]].location_temp(1)-myvertexs[index_vertex_now[0]].location_temp(1)*myvertexs[index_vertex_now[1]].location_temp(0);
    }
  area_temp=fabs(0.5*area_temp);
    
  double scale=(area-area_temp)/len_boundary_temp;

  size_t index_adjoin[2];
  for(i=0;i<num_vertex;i++)
    {
      index_adjoin[0]=(i+1)%num_vertex;
      index_adjoin[1]=((i-1)+num_vertex)%num_vertex;
      myvector norm=myvector(0,0,0);

      myvector norm_seg_1=(myvertexs[index_adjoin[0]].location_temp-myvertexs[i].location_temp).cross(myvector(0,0,1)); norm_seg_1.normalize();
      norm_seg_1*=(myvertexs[index_adjoin[0]].location_temp-myvertexs[i].location_temp).len();

      myvector norm_seg_2=(myvertexs[i].location_temp-myvertexs[index_adjoin[1]].location_temp).cross(myvector(0,0,1)); norm_seg_2.normalize();
      norm_seg_2*=(myvertexs[index_adjoin[1]].location_temp-myvertexs[i].location_temp).len();

      norm=norm_seg_1+norm_seg_2; norm.normalize();
      
      myvertexs[i].location=myvertexs[i].location_temp+scale*norm;
      myvertexs[i].velocity=myvertexs[i].velocity_temp+scale/dt*norm;
      // myvertexs[i].location=myvertexs[i].location_temp;
      // myvertexs[i].velocity=myvertexs[i].velocity_temp;
    }
  return 0;
}
