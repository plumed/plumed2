#include "ActionRegister.h"
#include "GenericGroup.h"

using namespace std;

namespace PLMD {

//+PLUMEDOC GENERIC GROUP
/**
Define a group of atoms

\par Syntax
\verbatim
GROUP LABEL=label ATOMS=x,y,z,...
\endverbatim
The label is associated to a group of atoms which is then automatically
expanded when used in multi-atoms options

*/
//+ENDPLUMEDOC

PLUMED_REGISTER_ACTION(GenericGroup,"GROUP")

GenericGroup::GenericGroup(const ActionOptions&ao):
  ActionAtomistic(ao)
{
  allowKeyword("ATOMS");
  int natoms=-1; unsigned ngrp=1;
  readActionAtomistic( natoms, ngrp );
  std::vector<double> domain(2,0.0);
  readActionWithExternalArguments( 0, domain );
  checkRead();
}

GenericGroup::~GenericGroup(){
  plumed.getAtoms().removeGroup(getLabel());
}

void GenericGroup::interpretAtomsKeyword( const std::vector<std::vector<unsigned> >& flist ){
  if( flist.size()!=1 ) error("cannot create multiple groups in a single line");

  log.printf("  created from the following list of atoms : " );
  for(unsigned j=0;j<flist[0].size();++j) log.printf("%s ", plumed.getAtoms().interpretIndex( flist[0][j] ).c_str() );
  log.printf("\n");
}

}
