#include "head.h"
#include "SWE_solver.h"
using namespace std;

void readcmdline(int argc, char* argv[],boost::property_tree::ptree &para_tree)
{
  size_t i;
  for(i=1;i<argc;++i)
    {
      string para_here=argv[i];
      size_t pos=para_here.find("=");
      if(pos!= string::npos)
	{
	  string key=para_here.substr(0,pos);
	  string value=para_here.substr(pos+1);
	  para_tree.put(key+".value",value);
	  printf("--[cmdline para] %s %s \n",key.c_str(),value.c_str());
	}
    }
  return;
}
int main(int argc, char *argv[])
{
  // test access point
  boost::property_tree::ptree para_tree;
  readcmdline(argc,argv,para_tree);
  string prog=para_tree.get<string>("prog.value");
  if(prog=="SWE_solver")
    {
      SWE_solver* mySWE_solver=new SWE_solver(para_tree);
      delete mySWE_solver;
    }
  return 0;
}
