#include "io.h"
using namespace std;

int io::saveTriAsVTK(std::vector<vertex> &myvertexs,std::vector<triangle > &mytri_face,const std::string name)
{
  FILE *fp;
  size_t i,j;
  printf("%s \n",name.c_str());
  fp=fopen(name.c_str(),"w");
  fprintf(fp,"# vtk DataFile Version 2.0\n");
  fprintf(fp,"triangle\n");
  fprintf(fp,"ASCII\n\n");
  fprintf(fp,"DATASET UNSTRUCTURED_GRID\n");
  fprintf(fp,"POINTS %d double\n",myvertexs.size());
  for(i=0;i<myvertexs.size();++i)
    {
      fprintf(fp,"%lf %lf %lf\n",myvertexs[i].location.x,myvertexs[i].location.y,myvertexs[i].location.z);
    }
  fprintf(fp,"CELLS %d %d\n",mytri_face.size(),mytri_face.size()*4);
  for(i=0;i<mytri_face.size();++i)
    {
      fprintf(fp,"3");
      for(j=0;j<3;++j)
	{
	  fprintf(fp," %d",mytri_face[i].index_vertex[j]);
	}
      fprintf(fp,"\n");

    }
  fprintf(fp,"CELL_TYPES %d\n",mytri_face.size());
  for(i=0;i<mytri_face.size();++i)
    {
      fprintf(fp,"5\n");
    }

  fprintf(fp,"POINT_DATA %d\n",myvertexs.size());
  fprintf(fp,"VECTORS v double\n");
  //  fprintf(fp,"LOOKUP_TABLE default\n"); //vectors 不要这句  scalars 要这句
  for(i=0;i<myvertexs.size();++i)
    {
      fprintf(fp,"%lf %lf 0\n",myvertexs[i].velocity(0),myvertexs[i].velocity(1));
    }
  /*
  fprintf(fp,"POINT_DATA %d\n",myvertexs.size());
  fprintf(fp,"SCALARS v_y double\n");
  fprintf(fp,"LOOKUP_TABLE default\n");
  for(i=0;i<myvertexs.size();++i)
    {
      fprintf(fp,"%lf\n",myvertexs[i].velocity(1));
      }*/
  fclose(fp);
  return 0;
}

int io::saveLineAsVTK(std::vector<vertex> &myvertexs,std::vector<line > &mylines,const std::string name)
{
  FILE *fp;
  size_t i,j;
  printf("%s \n",name.c_str());
  fp=fopen(name.c_str(),"w");
  
  fprintf(fp,"# vtk DataFile Version 2.0\n");
  fprintf(fp,"line\n");
  fprintf(fp,"ASCII\n\n");
  fprintf(fp,"DATASET UNSTRUCTURED_GRID\n");
  fprintf(fp,"POINTS %d double\n",myvertexs.size());
  
  for(i=0;i<myvertexs.size();++i)
    {
      fprintf(fp,"%lf %lf %lf\n",myvertexs[i].location.x,myvertexs[i].location.y,myvertexs[i].location.z);
    }
  fprintf(fp,"CELLS %d %d\n",mylines.size(),mylines.size()*3);
  for(i=0;i<mylines.size();++i)
    {
      fprintf(fp,"2");
      for(j=0;j<2;++j)
	{
	  fprintf(fp," %d",mylines[i].index_vertex[j]);
	}
      fprintf(fp,"\n");

    }
  fprintf(fp,"CELL_TYPES %d\n",mylines.size());
  for(i=0;i<mylines.size();++i)
    {
      fprintf(fp,"3\n");
    }

  fprintf(fp,"POINT_DATA %d\n",myvertexs.size());
  fprintf(fp,"VECTORS v double\n");
  //  fprintf(fp,"LOOKUP_TABLE default\n"); //vectors 不要这句  scalars 要这句
  for(i=0;i<myvertexs.size();++i)
    {
      fprintf(fp,"%lf %lf 0\n",myvertexs[i].velocity(0),myvertexs[i].velocity(1));
      }
  /*
  fprintf(fp,"POINT_DATA %d\n",myvertexs.size());
  fprintf(fp,"SCALARS v_y double\n");
  fprintf(fp,"LOOKUP_TABLE default\n");
  for(i=0;i<myvertexs.size();++i)
    {
      fprintf(fp,"%lf\n",myvertexs[i].velocity(1));
      }*/
   fclose(fp);
   return 0;
}

