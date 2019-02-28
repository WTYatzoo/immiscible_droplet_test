#include "vertex.h"

vertex::vertex()
{
  ;
}

vertex::vertex(const myvector &location)
{
  this->location_ori=this->location=this->location_temp=location;
  this->velocity=this->velocity_temp=myvector(0,0,0);
}

vertex::~vertex()
{
  ;
}
