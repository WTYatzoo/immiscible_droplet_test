#include "head.h"
#include "grid.h"

using namespace std;

grid::grid()
{
  ;
}

grid::~grid()
{
  ;
}

grid::grid(const size_t (&index_vertex)[2][2])
{
  size_t i,j;
  for(i=0;i<2;++i)
    {
      for(j=0;j<2;++j)
	{
	  this->index_vertex[i][j]=index_vertex[i][j];
	}
    }
}

