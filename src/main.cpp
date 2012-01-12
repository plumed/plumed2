#include "CLTool.h"

using namespace PLMD;

/// This main just wraps the CLTool::globalMain static
/// function. The object file generated from this .cpp
/// is the only part of the plumed library that should
/// not be linked with external MD codes, so as 
/// to avoid linker error.
int main(int argc,char**argv){
  return CLTool::globalMain(argc,argv);
}
