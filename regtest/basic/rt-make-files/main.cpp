#include "plumed/tools/File.h"
#include "plumed/tools/Tools.h"
#include "plumed/tools/Exception.h"
#include <fstream>
#include <iostream>

using namespace PLMD;

int main(){

//write
  int           out_int=-1;
  long int      out_longint=-2397083434877565865L;
  unsigned      out_unsigned=1;
  long unsigned out_longunsigned=2397083434877565865L;
  double        out_double=3.14;

  OFile out;
  out.open("test");
  out.printField("int",out_int);
  out.printField("long_int",out_longint);
  out.printField("unsigned",out_unsigned);
  out.printField("long_unsigned",out_longunsigned);
  out.printField("double",out_double);
  out.printField();
  out.close();

//read
  int           in_int=0;
  long int      in_longint=0;
  unsigned      in_unsigned=0;
  long unsigned in_longunsigned=0;
  double        in_double=0.0;

  IFile in;
  in.open("test");
  in.scanField("int",in_int);
  in.scanField("long_int",in_longint);
  in.scanField("unsigned",in_unsigned);
  in.scanField("long_unsigned",in_longunsigned);
  in.scanField("double",in_double);
  in.close();

//check
  plumed_assert(out_int==in_int);
  plumed_assert(out_longint==in_longint);
  plumed_assert(out_unsigned==in_unsigned);
  plumed_assert(out_longunsigned==in_longunsigned);
  plumed_assert(out_double==in_double);

  return 0;
}
