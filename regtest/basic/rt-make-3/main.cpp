#include "plumed/tools/File.h"
#include <fstream>
#include <iostream>

using namespace PLMD;

int main(){
  OFile out;
  out.open("test");
  out<<"prova1\n";
  out.close();

  out.enforceRestart();
  out.open("test");
  out<<"prova2\n";
  out.close();

  IFile in;
  std::string line;

  in.open("test");
  std::ofstream oo("output");
  while(in.getline(line)) oo<<line<<"\n";
  in.close();

  return 0;
}
