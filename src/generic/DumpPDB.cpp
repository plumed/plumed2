/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2023 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed.org for more information.

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "core/ActionWithArguments.h"
#include "core/ActionWithValue.h"
#include "core/ActionAtomistic.h"
#include "core/ActionPilot.h"
#include "core/ActionRegister.h"
#include "core/PlumedMain.h"
#include "tools/OFile.h"
#include "tools/h36.h"

namespace PLMD {
namespace generic {

//+PLUMEDOC PRINTANALYSIS DUMPPDB
/*
Output PDB file.


\par Examples

*/
//+ENDPLUMEDOC

class DumpPDB :
  public ActionWithArguments,
  public ActionPilot
{
  std::string fmt;
  std::string file;
  std::string description;
  std::vector<std::string> argnames;
  std::vector<AtomNumber> pdb_atom_indices;
  void buildArgnames();
  void printAtom( OFile& opdbf, const unsigned& anum, const Vector& pos, const double& m, const double& q ) const ;
public:
  static void registerKeywords( Keywords& keys );
  explicit DumpPDB(const ActionOptions&);
  void calculate() override {}
  void apply() override {}
  void update() override ;
};

PLUMED_REGISTER_ACTION(DumpPDB,"DUMPPDB")

void DumpPDB::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys );
  ActionPilot::registerKeywords( keys );
  ActionWithArguments::registerKeywords( keys ); keys.use("ARG");
  keys.add("optional","ATOMS","value containing positions of atoms that should be output");
  keys.add("compulsory","STRIDE","0","the frequency with which the atoms should be output");
  keys.add("compulsory","FILE","the name of the file on which to output these quantities");
  keys.add("compulsory","FMT","%f","the format that should be used to output real numbers");
  keys.add("optional","DESCRIPTION","the title to use for your PDB output");
  keys.add("optional","ATOM_INDICES","the indices of the atoms in your PDB output");
  keys.add("optional","ARG_NAMES","the names of the arguments that are being output");
}

DumpPDB::DumpPDB(const ActionOptions&ao):
  Action(ao),
  ActionWithArguments(ao),
  ActionPilot(ao)
{
  parse("FILE",file);
  if(file.length()==0) error("name out output file was not specified");
  std::vector<std::string> atoms; parseVector("ATOMS",atoms);
  if( atoms.size()>0 ) {
    std::vector<Value*> atom_args; interpretArgumentList( atoms, plumed.getActionSet(), this, atom_args );
    if( atom_args.size()!=1 ) error("only one action should appear in input to atom args");
    std::vector<Value*> args( getArguments() ); args.push_back( atom_args[0] ); requestArguments( args );
    std::vector<std::string> indices; parseVector("ATOM_INDICES",indices); std::vector<Value*> xvals,yvals,zvals,masv,chargev;
    ActionAtomistic::getAtomValuesFromPlumedObject(plumed,xvals,yvals,zvals,masv,chargev);
    ActionAtomistic::interpretAtomList( indices, xvals, this, pdb_atom_indices );
    log.printf("  printing atoms : ");
    for(unsigned i=0; i<indices.size(); ++i) log.printf("%d ", pdb_atom_indices[i].serial() );
    log.printf("\n");
    if( pdb_atom_indices.size()!=atom_args[0]->getShape()[1]/3 ) error("mismatch between size of matrix containing positions and vector of atom indices");
  }
  log.printf("  printing configurations to PDB file to file named %s \n", file.c_str() );
  parseVector("ARG_NAMES",argnames); if( argnames.size()==0 ) buildArgnames();
  parse("FMT",fmt); fmt=" "+fmt;
  if( getStride()==0 ) { setStride(0); log.printf("  printing pdb at end of calculation \n"); }

  parse("DESCRIPTION",description);
  if( getPntrToArgument(0)->getRank()==0 || getPntrToArgument(0)->hasDerivatives() ) error("argument for printing of PDB should be vector or matrix");
  unsigned nv=getPntrToArgument(0)->getShape()[0];
  for(unsigned i=0; i<getNumberOfArguments(); ++i) {
    if( getPntrToArgument(i)->getRank()==0 || getPntrToArgument(i)->hasDerivatives() ) error("argument for printing of PDB should be vector or matrix");
    if( getPntrToArgument(i)->getShape()[0]!=nv ) error("for printing to pdb files all arguments must have same number of values");
  }
}

void DumpPDB::printAtom( OFile& opdbf, const unsigned& anum, const Vector& pos, const double& m, const double& q ) const {
  std::array<char,6> at; const char* msg = h36::hy36encode(5,anum,&at[0]);
  plumed_assert(msg==nullptr) << msg; at[5]=0; double lunits = plumed.getUnits().getLength()/0.1;
  opdbf.printf("ATOM  %s  X   RES  %4u    %8.3f%8.3f%8.3f%6.2f%6.2f\n",
               &at[0], anum, lunits*pos[0], lunits*pos[1], lunits*pos[2], m, q );
}

void DumpPDB::buildArgnames() {
  unsigned nargs = getNumberOfArguments(); if( pdb_atom_indices.size()>0 ) nargs = nargs - 1;
  if( nargs<0 ) return;

  argnames.resize(0); unsigned nvals = getPntrToArgument(0)->getShape()[0];
  if( getPntrToArgument(0)->getRank()==2 ) nvals = getPntrToArgument(0)->getShape()[0];
  for(unsigned i=0; i<nargs; ++i) {
    if( getPntrToArgument(i)->getShape()[0]!=nvals ) error("all arguments should have same number of values");
    if( getPntrToArgument(i)->getRank()==1 ) {
      argnames.push_back( getPntrToArgument(i)->getName() );
    } else if( getPntrToArgument(i)->getRank()==2 ) {
      (getPntrToArgument(i)->getPntrToAction())->getMatrixColumnTitles( argnames );
    }
    getPntrToArgument(i)->buildDataStore();
  }
}

void DumpPDB::update() {
  OFile opdbf; opdbf.link(*this);
  opdbf.setBackupString("analysis");
  opdbf.open( file );
  std::size_t psign=fmt.find("%"); plumed_assert( psign!=std::string::npos );
  std::string descr2="%s=%-" + fmt.substr(psign+1) + " ";

  unsigned totargs = 0; unsigned nvals = getPntrToArgument(0)->getShape()[0];
  for(unsigned i=0; i<getNumberOfArguments(); ++i) {
    if( getPntrToArgument(i)->getShape()[0]!=nvals ) error("all arguments should have same number of values");
    if( getPntrToArgument(i)->getRank()==1 ) totargs += 1;
    else if( getPntrToArgument(i)->getRank()==2 ) totargs += getPntrToArgument(i)->getShape()[1];
  }
  if( totargs!=argnames.size() ) buildArgnames();

  if( description.length()>0 ) opdbf.printf("# %s AT STEP %lld TIME %f \n", description.c_str(), getStep(), getTime() );
  unsigned nargs = getNumberOfArguments(), atomarg=0;
  if( pdb_atom_indices.size()>0 ) {
    if( nargs>1 ) { atomarg = nargs - 1; nargs = nargs-1; }
    else nargs = 0;
  }
  for(unsigned i=0; i<nvals; ++i) {
    unsigned n=0;
    for(unsigned j=0; j<nargs; j++) {
      Value* thisarg=getPntrToArgument(j); opdbf.printf("REMARK ");
      if( thisarg->getRank()==1 ) {
        opdbf.printf( descr2.c_str(), argnames[n].c_str(), thisarg->get(i) ); n++;
      } else if( thisarg->getRank()==2 ) {
        unsigned ncols = thisarg->getShape()[1];
        for(unsigned k=0; k<ncols; ++k) { opdbf.printf( descr2.c_str(), argnames[n].c_str(), thisarg->get(i*ncols+k) ); n++; }
      }
      opdbf.printf("\n");
    }
    if( pdb_atom_indices.size()>0 ) {
      unsigned npos = pdb_atom_indices.size(); Vector pos;
      for(unsigned k=0; k<npos; ++k) {
        pos[0]=getPntrToArgument(atomarg)->get(npos*(3*i+0) + k);
        pos[1]=getPntrToArgument(atomarg)->get(npos*(3*i+1) + k);
        pos[2]=getPntrToArgument(atomarg)->get(npos*(3*i+2) + k);
        printAtom( opdbf, pdb_atom_indices[k].serial(), pos, 1.0, 1.0 );
      }
    }
    opdbf.printf("END\n");
  }
  opdbf.close();
}

}
}
