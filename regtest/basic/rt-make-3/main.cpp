#include "plumed/tools/File.h"
#include "plumed/tools/Tools.h"
#include "plumed/tools/Exception.h"
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

  int i=0;
  plumed_assert(Tools::convert("2*3",i));
  plumed_assert(i==6);

  i=0;
  plumed_assert(Tools::convert("(10+3)/13",i));
  plumed_assert(i==1);

  i=0;
  plumed_assert(Tools::convert("exp(log(7))",i));
  plumed_assert(i==7);

  i=0;
  plumed_assert(!Tools::convert("10.1",i));

  long int l=0;
  plumed_assert(Tools::convert("2397083434877565865",l));
  plumed_assert(l==2397083434877565865L);

  l=0;
  // this version is using lepton conversions and should fail:
  plumed_assert(!Tools::convert("1*2397083434877565865",l));

  return 0;
}
