#include "ActionAtomistic.h"
#include "ActionRegister.h"
#include <cstdio>

using namespace PLMD;
using namespace std;

namespace PLMD
{

//+PLUMEDOC GENERIC DUMPATOMS
/**
  Dump atom positions

\par syntax
\verbatim
DUMPATOMS [STRIDE=s] FILE=file.xyz ATOMS=list
\endverbatim

Listed atoms are written on file.xyz (currently only xyz format).
The positions are those stored in plumed exactly at the point where
the directive is put in the input file. This is relevant for directives
editing atom positions (e.g. \ref WHOLEMOLECULES).
*/
//+ENDPLUMEDOC

class GenericDumpAtoms:
  public ActionAtomistic
{
  FILE*fp;
public:
  GenericDumpAtoms(const ActionOptions&);
  ~GenericDumpAtoms();
  void interpretGroupsKeyword( const unsigned& natoms, const std::string& atomGroupName, const std::vector<std::vector<unsigned> >& groups );
  void interpretAtomsKeyword( const std::vector<std::vector<unsigned> >& flist );
  void updateNeighbourList( const double& cutoff, std::vector<bool>& skips ){ assert(false); }
  void calculate();
  void apply(){};
};

PLUMED_REGISTER_ACTION(GenericDumpAtoms,"DUMPATOMS")

GenericDumpAtoms::GenericDumpAtoms(const ActionOptions&ao):
  ActionAtomistic(ao)
{
  allowKeyword("ATOMS"); allowKeyword("GROUP" ); forbidKeyword("NL_CUTOFF");
  registerKeyword(1, "FILE", "file on which to output coordinates");
  int natoms=-1; unsigned ngrp=1; readActionAtomistic( natoms, ngrp );
  std::vector<double> domain(2,0.0);
  readActionWithExternalArguments( 0, domain );

  std::string file; parse("FILE",file);
  if( file.length()==0 ) error("specified input file makes no sense");
  fp=fopen(file.c_str(),"w");
  checkRead();
}

void GenericDumpAtoms::interpretGroupsKeyword( const unsigned& natoms, const std::string& atomGroupName, const std::vector<std::vector<unsigned> >& groups ){
  if( groups.size()!=1 ) error("cannot print atoms from multiple groups");

  if( atomGroupName!=getLabel() ){
     log.printf("  printing atoms specified in group %s\n", atomGroupName.c_str() );
  } else {
     log.printf("  printing the following list of atoms : " );
     for(unsigned j=0;j<groups[0].size();++j) log.printf("%s ", plumed.getAtoms().interpretIndex( groups[0][j] ).c_str() );
     log.printf("\n");
  }
} 

void GenericDumpAtoms::interpretAtomsKeyword( const std::vector<std::vector<unsigned> >& flist ){
  if( flist.size()!=1 ) error("cannot create multiple virtual atoms in a single line");

  log.printf("  printing the following list of atoms : " );
  for(unsigned j=0;j<flist[0].size();++j) log.printf("%s ", plumed.getAtoms().interpretIndex( flist[0][j] ).c_str() );
  log.printf("\n");
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
