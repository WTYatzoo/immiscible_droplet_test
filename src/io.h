#ifndef _IO_
#define _IO_
#include "head.h"
#include "vertex.h"
#include "triangle.h"
#include "line.h"
class io
{
 public:
  io(){}
  ~io(){}
  //name 不用引用，因为可能实参直接是"xxx" 
  int saveTriAsVTK(std::vector<vertex> &myvertexs,std::vector<triangle > &mytri_face,const std::string name);
  int saveLineAsVTK(std::vector<vertex> &myvertexs,std::vector<line > &mylines,const std::string name);
};
#endif
