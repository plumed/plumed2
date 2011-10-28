#include "ActionRegister.h"
#include "GenericGroup.h"

using namespace std;

namespace PLMD {

//+PLUMEDOC GROUP STATIC_GROUP
/**
Define a group of atoms

\par Example
The following contains a static group containing atoms 1-20.  Wherever the label
of the group appears after the GROUP keyword the specified list of atom will be used
to calculate the colvar.  
\verbatim
GROUP LABEL=label ATOMS=1-20
\endverbatim

*/
//+ENDPLUMEDOC

PLUMED_REGISTER_ACTION(GenericGroup,"STATIC_GROUP")

GenericGroup::GenericGroup(const ActionOptions&ao):
  ActionAtomistic(ao)
{
  allowKeyword("ATOMS"); forbidKeyword("UPDATE"); forbidKeyword("NL_CUTOFF");
  int natoms=-1; unsigned ngrp=1;
  readActionAtomistic( natoms, ngrp );
  std::vector<double> domain(2,0.0);
  readActionWithExternalArguments( 0, domain );
  checkRead();
}

GenericGroup::~GenericGroup(){
  plumed.getAtoms().removeGroup(getLabel());
}

void GenericGroup::getGroupDerivatives( std::vector<Vector>& derivatives ) const {
  return;
}

void GenericGroup::interpretAtomsKeyword( const std::vector<std::vector<unsigned> >& flist ){
  if( flist.size()!=1 ) error("cannot create multiple groups in a single line");

  log.printf("  created from the following list of atoms : " );
  for(unsigned j=0;j<flist[0].size();++j) log.printf("%s ", plumed.getAtoms().interpretIndex( flist[0][j] ).c_str() );
  log.printf("\n");
}

}
