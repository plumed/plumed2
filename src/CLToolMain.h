
#ifndef __PLUMED_CLToolMain_h
#define __PLUMED_CLToolMain_h
#include <cstdio>
#include <vector>
#include <string>
#include "WithCmd.h"


namespace PLMD{

class PlumedCommunicator;

/**
Interface to all the command-line tools.

This class just define an interface, and does not implement anything.
Inherits from this class to create a new command-line tool.
*/
class CLToolMain:
public WithCmd
{
/// arguments for command-line mode:
  int argc;
/// arguments for command-line mode:
  std::vector<std::string> argv;
  FILE*in;
  FILE*out;
  PlumedCommunicator&comm;
  static int run(int argc, char **argv,FILE*in,FILE*out,PlumedCommunicator&pc);
public:
  CLToolMain();
  ~CLToolMain();
/// virtual function returning a one-line descriptor for the tool
  void cmd(const std::string& key,void*val=NULL);
};

}


#endif
