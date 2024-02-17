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
#include "core/ActionWithVector.h"
#include "core/PlumedMain.h"
#include "core/ActionRegister.h"
#include "RMSDShortcut.h"
#include "tools/RMSD.h"
#include "tools/PDB.h"

namespace PLMD {
namespace colvar {

class RMSDVector : public ActionWithVector {

  bool firststep;
  bool squared;
  bool norm_weights;
  std::string type;
  std::vector<PLMD::RMSD> myrmsd;
  std::vector<double> align, displace;
  void setReferenceConfigurations();
public:
  static void registerKeywords(Keywords& keys);
  explicit RMSDVector(const ActionOptions&);
  unsigned getNumberOfDerivatives() override ;
  void performTask( const unsigned& current, MultiValue& myvals ) const override ;
  void calculate() override ;
};

PLUMED_REGISTER_ACTION(RMSDVector,"RMSD_VECTOR")

void RMSDVector::registerKeywords(Keywords& keys) {
  ActionWithVector::registerKeywords(keys); keys.use("ARG");
  keys.add("compulsory","TYPE","SIMPLE","the manner in which RMSD alignment is performed.  Should be OPTIMAL or SIMPLE.");
  keys.add("compulsory","ALIGN","1.0","the weights to use when aligning to the reference structure");
  keys.add("compulsory","DISPLACE","1.0","the weights to use when calculating the displacement from the reference structure");
  keys.addFlag("SQUARED",false," This should be set if you want mean squared displacement instead of RMSD ");
}

RMSDVector::RMSDVector(const ActionOptions&ao):
  Action(ao),
  ActionWithVector(ao),
  firststep(true)
{
  if( getPntrToArgument(0)->getRank()!=1 ) error("first argument should be vector");
  if( getPntrToArgument(1)->getRank()!=2 ) error("second argument should be matrix");
  if( getPntrToArgument(0)->getNumberOfValues()!=getPntrToArgument(1)->getShape()[1] ) error("mismatch between sizes of input vectors");
  if( getPntrToArgument(0)->getNumberOfValues()%3!=0 ) error("number of components in input arguments should be multiple of three");

  unsigned natoms = getPntrToArgument(0)->getNumberOfValues() / 3;
  type.assign("SIMPLE"); parse("TYPE",type); parseFlag("SQUARED",squared);
  std::vector<double> sqrtdisplace( natoms );
  RMSDShortcut::readAlignAndDisplace( this, norm_weights, align, displace, sqrtdisplace );

  std::vector<unsigned> shape( 1, getPntrToArgument(1)->getShape()[0] );
  addValue( shape ); setNotPeriodic(); myrmsd.resize( getPntrToArgument(1)->getShape()[0] );

  log.printf("  calculating RMSD distance between %d sets of %d atoms. Distance between vector %s of atoms and matrix of configurations in %s\n", 
                getPntrToArgument(1)->getShape()[0], natoms, getPntrToArgument(0)->getName().c_str(), getPntrToArgument(1)->getName().c_str() );
  log.printf("  method for alignment : %s \n",type.c_str() );
  if(squared)log.printf("  chosen to use SQUARED option for MSD instead of RMSD\n");
  else      log.printf("  using periodic boundary conditions\n");
}

unsigned RMSDVector::getNumberOfDerivatives() {
  return getPntrToArgument(0)->getNumberOfValues() + getPntrToArgument(1)->getNumberOfValues();
}

void RMSDVector::setReferenceConfigurations() {
  unsigned natoms = getPntrToArgument(0)->getShape()[0] / 3;
  Vector center; std::vector<Vector> pos( natoms );
  for(unsigned jconf=0; jconf<myrmsd.size(); ++jconf) {
      center.zero();
      for(unsigned i=0; i<pos.size(); ++i) {
          for(unsigned j=0; j<3; ++j) pos[i][j] = getPntrToArgument(1)->get( (3*jconf+j)*pos.size() + i );
          center+=pos[i]*align[i];
      }
      for(unsigned i=0; i<pos.size(); ++i) pos[i] -= center;
      myrmsd[jconf].clear(); myrmsd[jconf].set(align,displace,pos,type,true,norm_weights);
  }
}

// calculator
void RMSDVector::calculate() {
  if( firststep || !getPntrToArgument(1)->isConstant() ) { setReferenceConfigurations(); firststep=false; }
  runAllTasks();
}

void RMSDVector::performTask( const unsigned& current, MultiValue& myvals ) const {
  unsigned natoms = getPntrToArgument(0)->getShape()[0] / 3;

  std::vector<Vector>& pos( myvals.getFirstAtomVector() ); 
  std::vector<std::vector<Vector> > & allder( myvals.getFirstAtomDerivativeVector() ); 
  if( allder.size()!=2 ) allder.resize(2);
  std::vector<Vector> der( allder[0] );
  if( pos.size()!=natoms ) { pos.resize( natoms ); der.resize( natoms ); }
  for(unsigned i=0;i<pos.size();++i) {
      for(unsigned j=0; j<3; ++j) pos[i][j] = getPntrToArgument(0)->get( j*natoms + i );
  }
  double r = myrmsd[current].calculate( pos, der, squared ); unsigned ostrn = getConstPntrToComponent(0)->getPositionInStream();
  myvals.setValue( ostrn, r ); if( doNotCalculateDerivatives() ) return;

  for(unsigned i=0; i<natoms; i++){
      for(unsigned j=0; j<3; ++j ) { myvals.addDerivative( ostrn, j*natoms+i, der[i][j] ); myvals.updateIndex( ostrn, j*natoms+i ); }
  }
}

}
}



