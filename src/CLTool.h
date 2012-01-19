
#ifndef __PLUMED_CLTool_h
#define __PLUMED_CLTool_h
#include <cstdio>
#include <vector>
#include <string>

namespace PLMD{

class PlumedCommunicator;

class CLToolOptions{
  friend class CLToolRegister;
  std::vector<std::string> line;
public:
  CLToolOptions(const std::string &name):
    line(1,name) { }
};

/**
Interface to all the command-line tools.

This class just define an interface, and does not implement anything.
Inherits from this class to create a new command-line tool.
*/
class CLTool
{
public:
/// virtual function mapping to the specific main for each tool
  virtual int main(int argc, char **argv,FILE*in,FILE*out,PlumedCommunicator&pc)=0;
/// virtual function returning a one-line descriptor for the tool
  virtual std::string description()const{return "(no description available)";};
/// virtual destructor to allow inheritance
  virtual ~CLTool(){};
};

}


#endif
