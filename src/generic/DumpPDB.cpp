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
  std::vector<std::string> pdb_arg_names;
  std::vector<unsigned> pdb_atom_indices;
  void printAtom( OFile& opdbf, const unsigned& anum, const Vector& pos, const double& m, const double& q ) const ;
public:
  static void registerKeywords( Keywords& keys );
  explicit DumpPDB(const ActionOptions&);
  void calculate() override {}
  void apply() override {}
  void update() override ;
  void runFinalJobs() override ;
};

PLUMED_REGISTER_ACTION(DumpPDB,"DUMPPDB")

void DumpPDB::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys );
  ActionPilot::registerKeywords( keys );
  ActionWithArguments::registerKeywords( keys ); keys.use("ARG");
  keys.add("compulsory","STRIDE","0","the frequency with which the atoms should be output");
  keys.add("compulsory","FILE","the name of the file on which to output these quantities");
  keys.add("compulsory","FMT","%f","the format that should be used to output real numbers");
  keys.add("optional","DESCRIPTION","the title to use for your PDB output");
  keys.add("optional","ATOM_INDICES","the indices of the atoms in your PDB output");
  keys.add("optional","ARG_NAMES","the names to use for the arguments in the output PDB");
}

DumpPDB::DumpPDB(const ActionOptions&ao):
  Action(ao),
  ActionWithArguments(ao),
  ActionPilot(ao)
{
  parse("FILE",file);
  if(file.length()==0) error("name out output file was not specified");
  log.printf("  printing configurations to PDB file to file named %s \n", file.c_str() );
  parse("FMT",fmt); fmt=" "+fmt;
  if( getStride()==0 ) { setStride(0); log.printf("  printing pdb at end of calculation \n"); }

  parseVector("ATOM_INDICES",pdb_atom_indices); parseVector("ARG_NAMES",pdb_arg_names); parse("DESCRIPTION",description); 
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

void DumpPDB::update() {
  OFile opdbf; opdbf.link(*this);
  opdbf.setBackupString("analysis");
  opdbf.open( file ); unsigned nn=0;
  std::size_t psign=fmt.find("%"); plumed_assert( psign!=std::string::npos );  
  std::string descr2="%s=%-" + fmt.substr(psign+1) + " "; 
  // Add suitable code in here to print frames for paths here.  Gareth !!!!!
  std::vector<unsigned> argnums, posnums, matnums; bool use_real_arg_names = pdb_arg_names.size()==0;
  for(unsigned k=0;k<getNumberOfArguments();++k) { 
      if( getPntrToArgument(k)->getRank()==2 ) {
           if( matnums.size()>0 ) error("can only output one matrix at a time");
           matnums.push_back(k);
           if( pdb_atom_indices.size()==0 ) { pdb_atom_indices.resize( getPntrToArgument(k)->getShape()[1] / 3 ); for(unsigned i=0;i<pdb_atom_indices.size();++i) pdb_atom_indices[i]=i; }
           plumed_assert( pdb_atom_indices.size()==getPntrToArgument(k)->getShape()[1] / 3 );
      } else if( getPntrToArgument(k)->getName().find(".pos")!=std::string::npos ) {
           if( posnums.size()%3==0 && getPntrToArgument(k)->getName().find(".posx-")==std::string::npos ) error("x coordinate of input positions in wrong place");
           if( posnums.size()%3==1 && getPntrToArgument(k)->getName().find(".posy-")==std::string::npos ) error("y coordinate of input positions in wrong place");
           if( posnums.size()%3==2 && getPntrToArgument(k)->getName().find(".posz-")==std::string::npos ) error("z coordinate of input positions in wrong place");
           posnums.push_back(k);
      } else {
          if( use_real_arg_names ) pdb_arg_names.push_back( getPntrToArgument(k)->getName() );
          argnums.push_back( k );
      }
  }
  if( posnums.size()%3!=0 ) error("found misleading number of stored positions for output");
  if( pdb_atom_indices.size()==0 ) { pdb_atom_indices.resize( posnums.size() / 3 ); for(unsigned i=0;i<pdb_atom_indices.size();++i) pdb_atom_indices[i]=i; }

  if( description.length()>0 ) opdbf.printf("# %s AT STEP %d TIME %f \n", description.c_str(), getStep(), getTime() );
  unsigned nvals = getPntrToArgument(0)->getShape()[0]; Vector pos;
  for(unsigned j=0;j<nvals;++j) {
      if( argnums.size()>0 ) {
          opdbf.printf("REMARK ");
          for(unsigned k=0;k<argnums.size();++k) {
              Value* thisarg=getPntrToArgument(argnums[k]); opdbf.printf( descr2.c_str(), pdb_arg_names[k].c_str(), thisarg->get(j) );
          }
          opdbf.printf("\n");
      }
      if( posnums.size()==0 && matnums.size()==1 ) {
          unsigned npos = getPntrToArgument(matnums[0])->getShape()[1] / 3;
          for(unsigned k=0;k<npos;++k) {
              pos[0]=getPntrToArgument(matnums[0])->get(npos*(3*j+0) + k);
              pos[1]=getPntrToArgument(matnums[0])->get(npos*(3*j+1) + k);
              pos[2]=getPntrToArgument(matnums[0])->get(npos*(3*j+2) + k);
              printAtom( opdbf, pdb_atom_indices[k], pos, 1.0, 1.0 );
          }
      } else {
          unsigned npos = posnums.size() / 3;
          for(unsigned k=0;k<npos;++k) {
              pos[0]=getPntrToArgument(posnums[3*k+0])->get(j);
              pos[1]=getPntrToArgument(posnums[3*k+1])->get(j);
              pos[2]=getPntrToArgument(posnums[3*k+2])->get(j);
              printAtom( opdbf, pdb_atom_indices[k], pos, 1.0, 1.0 );
          }
      }
      opdbf.printf("END\n");
  }
  opdbf.close();
}

void DumpPDB::runFinalJobs() {
  if( getStride()>0 ) return ;
  setStride(1); update();
}

}
}
