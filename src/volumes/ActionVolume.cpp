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
  keys.addOutputComponent("_sum","SUM","the sum of all the colvars weighted by the function that determines if we are in the region");
  keys.addFlag("MEAN",false,"calculate the average value of the colvar inside the region of interest");
  keys.addOutputComponent("_mean","MEAN","the average values of the colvar in the region of interest");
  keys.add("optional","DATA","the label of an action that calculates multicolvars.  Weighted sums based on the location of the colvars calculated by this action will be calcualted");
  keys.add("optional","LESS_THAN","calcualte the number of colvars that are inside the region of interest and that are less than a certain threshold");
  keys.addOutputComponent("_lessthan","LESS_THAN","the number of cvs in the region of interest that are less than a certain threshold");
  keys.add("optional","MORE_THAN","calcualte the number of colvars that are inside the region of interest and that are greater that a certain threshold");
  keys.addOutputComponent("_morethan","MORE_THAN","the number of cvs in the region of interest that are more than a certain threshold");
  keys.add("optional","BETWEEN","calculate the number of colvars that are inside the region of interest and that have a CV value that is between a particular set of bounds");
  keys.addOutputComponent("_between","BETWEEN","the number of cvs in the region of interest that are within a certain range");
}

void ActionVolume::expandShortcut( const std::string& lab, const std::vector<std::string>& words,
                                   const std::map<std::string,std::string>& keys,
                                   std::vector<std::vector<std::string> >& actions ){
  if( keys.count("DATA") ) {
      std::string mc_lab = keys.find("DATA")->second;
      // Create the apprpriate volume object
      std::vector<std::string> vol_input; vol_input.push_back( lab + ":" );
      for(unsigned i=0;i<words.size();++i) vol_input.push_back( words[i] );
      vol_input.push_back( "ATOMS=" + mc_lab ); actions.push_back( vol_input );
      // Now create input for sums 
      if( keys.count("SUM") || keys.count("MEAN") ){
          std::vector<std::string> me_input; me_input.push_back( lab + "_prod:" );
          me_input.push_back("MATHEVAL"); me_input.push_back("ARG1=" + mc_lab); 
          me_input.push_back("ARG2=" + lab); me_input.push_back("FUNC=x*y"); 
          me_input.push_back("PERIODIC=NO"); actions.push_back( me_input );
          std::vector<std::string> input; 
          if( keys.count("SUM") ) input.push_back( lab + "_sum:" );
          else input.push_back( lab + "_numer:"); 
          input.push_back("COMBINE");
          input.push_back("ARG=" + lab + "_prod"); input.push_back("PERIODIC=NO"); actions.push_back( input );
      }
      if( keys.count("MEAN") ){
          // Calculate denominator
          std::vector<std::string> norm_in; norm_in.push_back( lab + "_norm:" );
          norm_in.push_back("COMBINE"); norm_in.push_back("ARG=" + lab); 
          norm_in.push_back("PERIODIC=NO"); actions.push_back( norm_in );
          // And calculate final quantity which is mean of these two actions
          std::vector<std::string> me_input2; me_input2.push_back( lab + "_mean:" );
          me_input2.push_back("MATHEVAL"); 
          if( keys.count("SUM") ) me_input2.push_back("ARG1=" + lab + "_sum" );
          else me_input2.push_back("ARG1=" + lab + "_numer"); 
          me_input2.push_back("ARG2=" + lab + "_norm"); me_input2.push_back("FUNC=x/y");
          me_input2.push_back("PERIODIC=NO"); actions.push_back( me_input2 );
      }
      if( keys.count("LESS_THAN") ){
          // Calculate number less than
          std::vector<std::string> lt_inp; lt_inp.push_back( mc_lab + "_" + lab + "_lt:" );
          lt_inp.push_back("LESS_THAN"); lt_inp.push_back("ARG=" + mc_lab );
          lt_inp.push_back("SWITCH=" + keys.find("LESS_THAN")->second  ); 
          actions.push_back( lt_inp );
          // And the matheval bit
          std::vector<std::string> me_input; me_input.push_back( lab + "_lt:" ); 
          me_input.push_back("MATHEVAL"); me_input.push_back("ARG1=" + mc_lab + "_" + lab + "_lt" );
          me_input.push_back("ARG2=" + lab); me_input.push_back("FUNC=x*y");
          me_input.push_back("PERIODIC=NO"); actions.push_back( me_input );
          // And the final sum
          std::vector<std::string> input; input.push_back( lab + "_lessthan:" ); input.push_back("COMBINE");
          input.push_back("ARG=" + lab + "_lt"); input.push_back("PERIODIC=NO"); actions.push_back( input );
      }
      if( keys.count("MORE_THAN") ){
          // Calculate number less than
          std::vector<std::string> lt_inp; lt_inp.push_back( mc_lab + "_" + lab + "_mt:" );
          lt_inp.push_back("LESS_THAN"); lt_inp.push_back("ARG=" + mc_lab );
          lt_inp.push_back("SWITCH=" + keys.find("LESS_THAN")->second  ); 
          actions.push_back( lt_inp );
          // And the matheval bit
          std::vector<std::string> me_input; me_input.push_back( lab + "_mt:" ); 
          me_input.push_back("MATHEVAL"); me_input.push_back("ARG1=" + mc_lab + "_" + lab + "_mt" );
          me_input.push_back("ARG2=" + lab); me_input.push_back("FUNC=x*y");
          me_input.push_back("PERIODIC=NO"); actions.push_back( me_input );
          // And the final sum
          std::vector<std::string> input; input.push_back( lab + "_morethan:" ); input.push_back("COMBINE");
          input.push_back("ARG=" + lab + "_mt"); input.push_back("PERIODIC=NO"); actions.push_back( input );
      }
      if( keys.count("BETWEEN") ){
          // Calculate number less than
          std::vector<std::string> lt_inp; lt_inp.push_back( mc_lab + "_" + lab + "_bt:" );
          lt_inp.push_back("LESS_THAN"); lt_inp.push_back("ARG=" + mc_lab );
          lt_inp.push_back("SWITCH=" + keys.find("LESS_THAN")->second  ); 
          actions.push_back( lt_inp );
          // And the matheval bit
          std::vector<std::string> me_input; me_input.push_back( lab + "_bt:" ); 
          me_input.push_back("MATHEVAL"); me_input.push_back("ARG1=" + mc_lab + "_" + lab + "_bt" );
          me_input.push_back("ARG2=" + lab); me_input.push_back("FUNC=x*y");
          me_input.push_back("PERIODIC=NO"); actions.push_back( me_input );
          // And the final sum
          std::vector<std::string> input; input.push_back( lab + "_between:" ); input.push_back("COMBINE");
          input.push_back("ARG=" + lab + "_bt"); input.push_back("PERIODIC=NO"); actions.push_back( input );
      }
  } else if( keys.count("SUM") ) { 
      std::vector<std::string> mc_line; mc_line.push_back( lab + "_vols:" );
      for(unsigned i=0;i<words.size();++i) mc_line.push_back(words[i]);
      actions.push_back( mc_line );
      std::vector<std::string> input; input.push_back( lab + ":" ); input.push_back("COMBINE");
      input.push_back("ARG=" + lab + "_vols"); input.push_back("PERIODIC=NO"); actions.push_back( input );
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

  // Check if we have any dependencies on ActionWithValue objects
  if( getDependencies().size()>0 ){
      std::vector<ActionWithValue*> f_actions;
      for(unsigned i=0;i<getDependencies().size();++i){
          ActionWithValue* av = dynamic_cast<ActionWithValue*>( getDependencies()[i] );
          if( av ) f_actions.push_back( av );
      }
      plumed_assert( f_actions.size()<2 );
      if( f_actions.size()==1 ){ 
          std::vector<std::string> empty(1); empty[0] = f_actions[0]->getLabel();
          f_actions[0]->addActionToChain( empty, this ); 
      }
  }

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

void ActionVolume::buildCurrentTaskList( std::vector<unsigned>& tflags ) {
  tflags.assign(tflags.size(),1);  // Can surely do something more smart here
}

void ActionVolume::calculate(){
  if( actionInChain() ) return;
  setupRegions(); runAllTasks();
}

void ActionVolume::prepareForTasks(){
  if( actionInChain() ) retrieveAtoms();
  setupRegions(); ActionWithValue::prepareForTasks();
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
      outvals.addDerivative( getPntrToOutput(0)->getPositionInStream(), 3*curr+0, wdf[0] );
      outvals.updateIndex( getPntrToOutput(0)->getPositionInStream(), 3*curr+0 );
      outvals.addDerivative( getPntrToOutput(0)->getPositionInStream(), 3*curr+1, wdf[1] );
      outvals.updateIndex( getPntrToOutput(0)->getPositionInStream(), 3*curr+1 );
      outvals.addDerivative( getPntrToOutput(0)->getPositionInStream(), 3*curr+2, wdf[2] );
      outvals.updateIndex( getPntrToOutput(0)->getPositionInStream(), 3*curr+2 ); 
      // Add derivatives with respect to reference positions
      unsigned vbase = 3*getFullNumberOfTasks();  
      for(unsigned i=0;i<refders.size();++i){
         outvals.addDerivative( getPntrToOutput(0)->getPositionInStream(), vbase, refders[i][0] );
         outvals.updateIndex( getPntrToOutput(0)->getPositionInStream(), vbase ); vbase++;
         outvals.addDerivative( getPntrToOutput(0)->getPositionInStream(), vbase, refders[i][1] );
         outvals.updateIndex( getPntrToOutput(0)->getPositionInStream(), vbase ); vbase++;
         outvals.addDerivative( getPntrToOutput(0)->getPositionInStream(), vbase, refders[i][2] );
         outvals.updateIndex( getPntrToOutput(0)->getPositionInStream(), vbase ); vbase++;
      }
      // Add virial 
      outvals.addDerivative( getPntrToOutput(0)->getPositionInStream(), vbase, vir(0,0) ); vbase++;
      outvals.addDerivative( getPntrToOutput(0)->getPositionInStream(), vbase, vir(0,1) ); vbase++;
      outvals.addDerivative( getPntrToOutput(0)->getPositionInStream(), vbase, vir(0,2) ); vbase++;
      outvals.addDerivative( getPntrToOutput(0)->getPositionInStream(), vbase, vir(1,0) ); vbase++;
      outvals.addDerivative( getPntrToOutput(0)->getPositionInStream(), vbase, vir(1,1) ); vbase++;
      outvals.addDerivative( getPntrToOutput(0)->getPositionInStream(), vbase, vir(1,2) ); vbase++;
      outvals.addDerivative( getPntrToOutput(0)->getPositionInStream(), vbase, vir(2,0) ); vbase++;
      outvals.addDerivative( getPntrToOutput(0)->getPositionInStream(), vbase, vir(2,1) ); vbase++;
      outvals.addDerivative( getPntrToOutput(0)->getPositionInStream(), vbase, vir(2,2) ); vbase++;
      vbase = 3*getNumberOfAtoms();
      for(unsigned i=0;i<9;++i) outvals.updateIndex( getPntrToOutput(0)->getPositionInStream(), vbase + i );
  }
}

void ActionVolume::apply(){
  if( doNotCalculateDerivatives() ) return;
  std::fill(forcesToApply.begin(),forcesToApply.end(),0);
  if( getForcesFromValues( forcesToApply ) ) setForcesOnAtoms( forcesToApply );
}

}
}
