
#ifndef __PLUMED_CLToolMain_h
#define __PLUMED_CLToolMain_h
#include <cstdio>
#include <vector>
#include <string>
#include "WithCmd.h"


namespace PLMD{

class PlumedCommunicator;

/**
Class providing access to command line tools.

This class provides an interface using the "cmd()" syntax to all the
command-line tools. In this manner, it allows all the registered
tools to be called directly from a PlumedMain object using proper
commands. It is only accessed via the cmd() function.
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
/**
Send messages to the CLToolMain.

Messages are used to set the MPI communicator, input, output etc. See the CLToolMain.cpp
file for details. A typical usage is:
\verbatim
#include "CLToolMain.h"
int main(int argc,char**argv){
  PLMD::CLToolMain cltoolMain;
  cltoolMain.cmd("setArgc",&argc);
  cltoolMain.cmd("setArgv",argv);
  int ret;
  cltoolMain.cmd("run",&ret);
  return ret;
}
\endverbatim
This will run the tool registered with name argv[1] with options argv[2]...argv[argc-1].

*/
  void cmd(const std::string& key,void*val=NULL);
};

}


#endif
