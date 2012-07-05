
#ifndef __PLUMED_CLToolMain_h
#define __PLUMED_CLToolMain_h
#include <cstdio>
#include <vector>
#include <string>
#include "WithCmd.h"


namespace PLMD{

class PlumedCommunicator;

/**
Class providing cmd() access to command line tools.

This class provides an interface using the "cmd()" syntax to all the
command-line tools.
It is only accessed via the cmd() function, which can
be used to set the arguments, communicators and IO descriptors and 
to run the tool.
It can run all the tools registered via the PLUMED_REGISTER_CLTOOL macro,
or the scripts which are located in PLUMED_ROOT/scripts.

A typical usage is:
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

This class is also used in the \ref PlumedMain class to provide
the same functionalities through the external plumed interface, which
is available also for C and FORTRAN. Thus, the preferred approach is to do something like
\verbatim
#include "Plumed.h"
int main(int argc,char**argv){
  PLMD::Plumed p;
  p.cmd("CLTool setArgc",&argc);
  p.cmd("CLTool setArgv",argv);
  int ret;
  p.cmd("CLTool run",&ret);
  return ret;
}
\endverbatim

See the file \ref main.cpp for a similar example.

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
*/
  void cmd(const std::string& key,void*val=NULL);
};

}


#endif
