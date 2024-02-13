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
#include "core/ActionWithValue.h"
#include "core/ActionWithArguments.h"
#include "core/ActionRegister.h"
#include "RMSDShortcut.h"
#include "tools/PDB.h"

namespace PLMD {
namespace colvar {

class DisplacementRMSD : 
public ActionWithValue,
public ActionWithArguments {

  bool firststep;
  bool squared;
  bool norm_weights;
  std::string type;
  PLMD::RMSD myrmsd;
  std::vector<double> align, displace, sqrtdisplace;
  std::vector<Vector> pos, direction, der;
public:
  explicit DisplacementRMSD(const ActionOptions&);
  void calculate() override;
  static void registerKeywords(Keywords& keys);
  unsigned getNumberOfDerivatives() override { return 2*getPntrToArgument(0)->getNumberOfValues(); }
  void apply() override ;
};

PLUMED_REGISTER_ACTION(DisplacementRMSD,"RMSD_DISPLACEMENT_VECTOR")

void DisplacementRMSD::registerKeywords(Keywords& keys) {
  Action::registerKeywords( keys ); ActionWithValue::registerKeywords( keys ); ActionWithArguments::registerKeywords( keys ); keys.use("ARG");
  keys.add("compulsory","ALIGN","1.0","the weights to use when aligning to the reference structure");
  keys.add("compulsory","DISPLACE","1.0","the weights to use when calculating the displacement from the reference structure");
  keys.add("compulsory","TYPE","SIMPLE","the manner in which RMSD alignment is performed.  Should be OPTIMAL or SIMPLE.");
  keys.addFlag("SQUARED",false," This should be set if you want mean squared displacement instead of RMSD ");
  keys.addFlag("UNORMALIZED",false,"by default the mean sequare deviation or root mean square deviation is calculated.  If this option is given no averaging is done");
  keys.addOutputComponent("disp","default","the vector of displacements for the atoms");
  keys.addOutputComponent("dist","default","the RMSD distance the atoms have moved");
}

DisplacementRMSD::DisplacementRMSD(const ActionOptions&ao):
  Action(ao),
  ActionWithValue(ao),
  ActionWithArguments(ao),
  firststep(true)
{
  bool unorm=false; parseFlag("UNORMALIZED",unorm); norm_weights=!unorm;
  if( getPntrToArgument(0)->getRank()!=1 || getPntrToArgument(1)->getRank()!=1 ) error("arguments should be vectors");
  if( getPntrToArgument(0)->getNumberOfValues()!=getPntrToArgument(1)->getNumberOfValues() ) error("mismatch between sizes of input vectors");
  if( getPntrToArgument(0)->getNumberOfValues()%3!=0 ) error("number of components in input arguments should be multiple of three");

  requestArguments( getArguments() );
  type.assign("SIMPLE"); parse("TYPE",type); parseFlag("SQUARED",squared);
  RMSDShortcut::readAlignAndDisplace( this, norm_weights, align, displace, sqrtdisplace );

  addComponentWithDerivatives( "dist" ); componentIsNotPeriodic("dist");
  std::vector<unsigned> shape( 1, getPntrToArgument(0)->getShape()[0] );
  addComponent( "disp", shape ); componentIsNotPeriodic("disp");
  getPntrToComponent(1)->buildDataStore();

  unsigned natoms = getPntrToArgument(0)->getNumberOfValues()/3; pos.resize(natoms); der.resize(natoms); direction.resize(natoms);
  log.printf("  calculating RMSD distance between two sets of %d atoms in vectors %s and %s\n", natoms, getPntrToArgument(1)->getName().c_str(), getPntrToArgument(0)->getName().c_str() );
  log.printf("  method for alignment : %s \n",type.c_str() );
  if(squared)log.printf("  chosen to use SQUARED option for MSD instead of RMSD\n");
  else      log.printf("  using periodic boundary conditions\n");
}


// calculator
void DisplacementRMSD::calculate() {
  if( firststep || !getPntrToArgument(1)->isConstant() ) { 
      RMSDShortcut::setReferenceConfiguration( 0, getPntrToArgument(1), align, displace, type, norm_weights, myrmsd ); 
      firststep=false;
  }
  unsigned natoms = pos.size();
  for(unsigned i=0;i<natoms;++i) {
      for(unsigned j=0; j<3; ++j) pos[i][j] = getPntrToArgument(0)->get( j*natoms + i );
  }
  double r = RMSDShortcut::calculateDisplacement( type, align, displace, sqrtdisplace, pos, myrmsd, direction, der, squared );
  // Set the total RMSD
  getPntrToComponent(0)->set(r);
  // Set the displacement vector
  Value* myval=getPntrToComponent(1);
  for(unsigned j=0; j<3; ++j) {
      for(unsigned i=0;i<natoms;++i) myval->set( j*natoms + i, direction[i][j] );
  }
  if( doNotCalculateDerivatives() ) return ;

  myval = getPntrToComponent(0);
  for(unsigned i=0; i<der.size(); ++i) {
      for(unsigned j=0; j<3; ++j) myval->setDerivative( j*natoms + i, der[i][j] );
  }
}

void DisplacementRMSD::apply() {
  if( doNotCalculateDerivatives() ) return;

  bool hasforce=false; std::vector<double> forces( getNumberOfDerivatives(), 0 );
  if( getPntrToComponent(0)->forcesWereAdded() ) { hasforce=true; getPntrToComponent(0)->applyForce( forces ); } 
  if( getPntrToComponent(1)->forcesWereAdded() ) {
      Value* mydisp = getPntrToComponent(1);
      hasforce=true; RMSDShortcut::addDisplacementForces( type, align, displace, sqrtdisplace, pos, myrmsd, direction, der, mydisp, squared );
      for(unsigned i=0; i<forces.size(); ++i) forces[i] += mydisp->getForce(i); 
  }
  if( hasforce ) { unsigned mm=0; addForcesOnArguments( 0, forces, mm, getLabel() ); }
}

}
}



