#include "PlumedException.h"
#include "CLTool.h"
#include "CLToolRegister.h"
#include <string>
#include <cstdlib>
#include <cstdio>
#include <iostream>

using namespace std;
using namespace PLMD;

/**
This is the entry point to the command line tools
included in the plumed library.

\todo
Features implemented here should be merged with the tools/plumed script.
I mean: the final user should see a single "plumed" executable
which should be usable as 
plumed driver
plumed mklib
...
Probably it is better to transfer all the parsing into
the routine down here, and then to call extra tools which
need bash/patch/other commands using system()
In this manner, all the features which do not require system()
will be available also e.g. on Windows (I am thinking
about VMD plugin etc)
*/

int CLTool::globalMain(int argc, char **argv,FILE*in,FILE*out){
  int i;
  for(i=1;i<argc;i++){
    string arg(argv[i]);
    if(arg.length()>0){
      if(arg[0]=='-'){
        plumed_merror("Unknown option " + arg);
      } else break;
    }
  }
  if(i==argc) {
    fprintf(out,"%s","nothing to do\n");
    exit(0);
  }
  CLTool *cl;
  string command(argv[i]);
  cl=cltoolRegister().create(command);
  if(!cl) plumed_merror("unknown command " + command);

  int ret=cl->main(argc-i,&argv[i],in,out);

  delete cl;

  return ret;
}
