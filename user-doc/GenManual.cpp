#include <iostream>
#include "../src/Plumed.h"

int main ( int argc, char* argv[] ) {

  if ( argc!=2 ){
     std::cout<<"Wrong number of arguments to GenManual routine "<<argc<<std::endl;
     for(unsigned i=0;i<argc;++i) std::cout<<argv[i]<<" ";
     std::cout<<std::endl;
     abort();
  }
  std::string actionName=argv[1];
  std::cout<<"Creating manual for argument "<<actionName<<std::endl;
  std::string fname=actionName + ".man";
  FILE* ofile; ofile=fopen(fname.c_str(),"w+");

  plumed plumedmain;
  plumedmain=plumed_create();
  plumed_cmd(plumedmain,"setLog",ofile); 
  char* aname=const_cast<char*>( actionName.c_str() );
  plumed_cmd(plumedmain,"GenerateManual",aname);
  plumed_finalize(plumedmain); fclose(ofile);
  return 1;
}
