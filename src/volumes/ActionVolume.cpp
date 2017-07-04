/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2017 The plumed team
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
#include "ActionVolume.h"

namespace PLMD {
namespace volumes {

void ActionVolume::shortcutKeywords( Keywords& keys ) {
  keys.addFlag("SUM",false,"calculate the sum of all the quantities.");
}

void ActionVolume::expandShortcut( const std::string& lab, const std::vector<std::string>& words,
                                   const std::map<std::string,std::string>& keys,
                                   std::vector<std::vector<std::string> >& actions ){
  if( keys.count("SUM") ) { 
      std::vector<std::string> mc_line; mc_line.push_back( lab + "_vols:" );
      for(unsigned i=0;i<words.size();++i) mc_line.push_back(words[i]);
      actions.push_back( mc_line );
      std::vector<std::string> input; input.push_back( lab + ":" ); input.push_back("SUM");
      input.push_back("ARG=" + lab + "_vols"); actions.push_back( input );
  }
}

void ActionVolume::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys ); 
  ActionAtomistic::registerKeywords( keys );
  ActionWithValue::registerKeywords( keys );
  keys.add("atoms","ATOMS","the group of atoms that you would like to investigate");
  keys.add("compulsory","SIGMA","the width of the function to be used for kernel density estimation");
  keys.add("compulsory","KERNEL","gaussian","the type of kernel function to be used");
  keys.addFlag("OUTSIDE",false,"calculate quantities for colvars that are on atoms outside the region of interest");
}

ActionVolume::ActionVolume(const ActionOptions&ao):
  Action(ao),
  ActionAtomistic(ao),
  ActionWithValue(ao)
{
  std::vector<AtomNumber> atoms; parseAtomList("ATOMS",atoms);
  if( atoms.size()==0 ) error("no atoms were specified");
  log.printf("  examining positions of atoms ");
  for(unsigned i=0;i<atoms.size();++i){ log.printf(" %d", atoms[i].serial() );  addTaskToList( i ); }
  log.printf("\n"); ActionAtomistic::requestAtoms( atoms );

  parseFlag("OUTSIDE",not_in); sigma=0.0;
  if( keywords.exists("SIGMA") ) parse("SIGMA",sigma);
  if( keywords.exists("KERNEL") ) parse("KERNEL",kerneltype);

  if( getFullNumberOfTasks()==1 ){ ActionWithValue::addValueWithDerivatives(); }
  else { std::vector<unsigned> shape(1); shape[0]=getFullNumberOfTasks(); ActionWithValue::addValue( shape ); }
  setNotPeriodic();
}

void ActionVolume::requestAtoms( const std::vector<AtomNumber> & a ) {
  std::vector<AtomNumber> all_atoms( getAbsoluteIndexes() );
  for(unsigned i=0;i<a.size();++i) all_atoms.push_back( a[i] );
  ActionAtomistic::requestAtoms( all_atoms ); forcesToApply.resize( 3*all_atoms.size()+9 );
}

void ActionVolume::buildCurrentTaskList( std::vector<unsigned>& tflags ) const {
  tflags.assign(tflags.size(),1);  // Can surely do something more smart here
}

void ActionVolume::calculate(){
  setupRegions(); runAllTasks();
}

void ActionVolume::performTask( const unsigned& curr, MultiValue& outvals ) const {
  Vector wdf; Tensor vir; std::vector<Vector> refders( getNumberOfAtoms()-getFullNumberOfTasks() );
  double weight=calculateNumberInside( ActionAtomistic::getPosition(curr), wdf, vir, refders );
  if( not_in ) {
    weight = 1.0 - weight; wdf *= -1.; vir *=-1;
    for(unsigned i=0; i<refders.size(); ++i) refders[i]*=-1;
  }
  outvals.setValue( getPntrToOutput(0)->getPositionInStream(), weight );
  if( !doNotCalculateDerivatives() ){
      // Atom position
      outvals.putIndexInActiveArray( 3*curr+0 );  
      outvals.addDerivative( getPntrToOutput(0)->getPositionInStream(), 3*curr+0, wdf[0] );
      outvals.putIndexInActiveArray( 3*curr+1 );
      outvals.addDerivative( getPntrToOutput(0)->getPositionInStream(), 3*curr+1, wdf[1] );
      outvals.putIndexInActiveArray( 3*curr+2 );
      outvals.addDerivative( getPntrToOutput(0)->getPositionInStream(), 3*curr+2, wdf[2] );
      // Add derivatives with respect to reference positions
      unsigned vbase = 3*getFullNumberOfTasks();  
      for(unsigned i=0;i<refders.size();++i){
         outvals.addDerivative( getPntrToOutput(0)->getPositionInStream(), vbase, refders[i][0] );
         outvals.putIndexInActiveArray( vbase ); vbase++;
         outvals.addDerivative( getPntrToOutput(0)->getPositionInStream(), vbase, refders[i][1] );
         outvals.putIndexInActiveArray( vbase ); vbase++;
         outvals.addDerivative( getPntrToOutput(0)->getPositionInStream(), vbase, refders[i][2] );
         outvals.putIndexInActiveArray( vbase ); vbase++;
      }
      // Add virial 
      for(unsigned i=0;i<9;++i) outvals.putIndexInActiveArray( vbase + i );
      outvals.addDerivative( getPntrToOutput(0)->getPositionInStream(), vbase, vir(0,0) ); vbase++;
      outvals.addDerivative( getPntrToOutput(0)->getPositionInStream(), vbase, vir(0,1) ); vbase++;
      outvals.addDerivative( getPntrToOutput(0)->getPositionInStream(), vbase, vir(0,2) ); vbase++;
      outvals.addDerivative( getPntrToOutput(0)->getPositionInStream(), vbase, vir(1,0) ); vbase++;
      outvals.addDerivative( getPntrToOutput(0)->getPositionInStream(), vbase, vir(1,1) ); vbase++;
      outvals.addDerivative( getPntrToOutput(0)->getPositionInStream(), vbase, vir(1,2) ); vbase++;
      outvals.addDerivative( getPntrToOutput(0)->getPositionInStream(), vbase, vir(2,0) ); vbase++;
      outvals.addDerivative( getPntrToOutput(0)->getPositionInStream(), vbase, vir(2,1) ); vbase++;
      outvals.addDerivative( getPntrToOutput(0)->getPositionInStream(), vbase, vir(2,2) ); vbase++;
      outvals.completeUpdate();
  }
}

void ActionVolume::apply(){
  if( doNotCalculateDerivatives() ) return;
  std::fill(forcesToApply.begin(),forcesToApply.end(),0);
  if( getForcesFromValues( forcesToApply ) ) setForcesOnAtoms( forcesToApply );
}

}
}
