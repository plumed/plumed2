#include "ActionAtomistic.h"
#include "ActionRegister.h"
#include <cstdio>

using namespace PLMD;
using namespace std;

namespace PLMD
{

//+PLUMEDOC GENERIC DUMPATOMS
/**

The atoms listed are written out on a file (currently only xyz format).
The positions are those stored in plumed exactly at the point where
the directive is put in the input file. This is useful for debugging directives that
edit atom positions (e.g. \ref WHOLEMOLECULES).

\par Example
This will print out the positions of atoms 1-100 on the file file.xyz every 10 steps
\verbatim
DUMPATOMS ATOMS=1-100 FILE=file.xyz STRIDE=10
\endverbatim

*/
//+ENDPLUMEDOC

class GenericDumpAtoms:
  public ActionAtomistic
{
  FILE*fp;
public:
  GenericDumpAtoms(const ActionOptions&);
  ~GenericDumpAtoms();
  void calculate();
  void apply(){};
};

PLUMED_REGISTER_ACTION(GenericDumpAtoms,"DUMPATOMS")

GenericDumpAtoms::GenericDumpAtoms(const ActionOptions&ao):
  ActionAtomistic(ao)
{
  registerKeyword(2, "ATOMS", "the atom indices whose positions you would like to print out");
  registerKeyword(1, "FILE", "file on which to output coordinates");
  setNeighbourListStyle("none");

  // Read everything in ActionAtomistic
  readActionAtomistic(); parseAtomList("ATOMS", 0 ); printAllAtoms("printing atoms"); 

  std::vector<double> domain(2,0.0);
  readActionWithExternalArguments( 0, domain );

  std::string file; parse("FILE",file);
  if( file.length()==0 ) error("specified input file makes no sense");
  fp=fopen(file.c_str(),"w");
  checkRead();
}

void GenericDumpAtoms::calculate(){
  fprintf(fp,"%d\n",getNumberOfAtoms());
  const Tensor & t(getBox()); Pbc tpbc; tpbc.setBox(t);
  if(tpbc.isOrthorombic()){
    fprintf(fp," %f %f %f\n",t(0,0),t(1,1),t(2,2));
  }else{
    fprintf(fp," %f %f %f %f %f %f %f %f %f\n",
                 t(0,0),t(0,1),t(0,2),
                 t(1,0),t(1,1),t(1,2),
                 t(2,0),t(2,1),t(2,2)
           );
  }
  for(unsigned i=0;i<getNumberOfAtoms();++i){
    fprintf(fp,"X %f %f %f\n",getPositions(i)(0),getPositions(i)(1),getPositions(i)(2));
  }
}

GenericDumpAtoms::~GenericDumpAtoms(){
  fclose(fp);
}
  

}
