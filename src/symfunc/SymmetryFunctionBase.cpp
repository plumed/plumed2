/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015-2017 The plumed team
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
#include "SymmetryFunctionBase.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"
#include "multicolvar/MultiColvarBase.h"

namespace PLMD {
namespace symfunc {

void SymmetryFunctionBase::shortcutKeywords( Keywords& keys ) {
  keys.add("atoms-3","SPECIES","this keyword is used for colvars such as coordination number. In that context it specifies that plumed should calculate "
               "one coordination number for each of the atoms specified.  Each of these coordination numbers specifies how many of the "
               "other specified atoms are within a certain cutoff of the central atom.  You can specify the atoms here as another multicolvar "
               "action or using a MultiColvarFilter or ActionVolume action.  When you do so the quantity is calculated for those atoms specified "
               "in the previous multicolvar.  This is useful if you would like to calculate the Steinhardt parameter for those atoms that have a "
               "coordination number more than four for example");
  keys.add("atoms-4","SPECIESA","this keyword is used for colvars such as the coordination number.  In that context it species that plumed should calculate "
               "one coordination number for each of the atoms specified in SPECIESA.  Each of these cooordination numbers specifies how many "
               "of the atoms specifies using SPECIESB is within the specified cutoff.  As with the species keyword the input can also be specified "
               "using the label of another multicolvar");
  keys.add("atoms-4","SPECIESB","this keyword is used for colvars such as the coordination number.  It must appear with SPECIESA.  For a full explanation see "
               "the documentation for that keyword");
  keys.add("compulsory","NN","6","The n parameter of the switching function ");
  keys.add("compulsory","MM","0","The m parameter of the switching function; 0 implies 2*NN");
  keys.add("compulsory","D_0","0.0","The d_0 parameter of the switching function");
  keys.add("compulsory","R_0","The r_0 parameter of the switching function");
  keys.add("optional","SWITCH","the switching function that it used in the construction of the contact matrix");
  multicolvar::MultiColvarBase::shortcutKeywords( keys );
}

void SymmetryFunctionBase::expandMatrix( const bool& components, const std::string& lab, const std::vector<std::string>& words,
                                         const std::map<std::string,std::string>& keys,
                                         std::vector<std::vector<std::string> >& actions ){
  std::vector<std::string> matinp; matinp.push_back( lab + "_mat:" ); matinp.push_back("CONTACT_MATRIX");
  if( keys.count("SPECIES") ) {
      matinp.push_back("GROUP=" + keys.find("SPECIES")->second ); 
  } else if( keys.count("SPECIESA") ) {
      matinp.push_back("GROUPA=" + keys.find("SPECIESA")->second ); matinp.push_back("GROUPB=" + keys.find("SPECIESB")->second );
  }
  if( keys.count("SWITCH") ) { 
      matinp.push_back("SWITCH=" + keys.find("SWITCH")->second );
  } else if( keys.count("R_0") ) {
      matinp.push_back("R_0=" + keys.find("R_0")->second );
      matinp.push_back("D_0=" + keys.find("D_0")->second );
      matinp.push_back("NN=" + keys.find("NN")->second );
      matinp.push_back("MM=" + keys.find("MM")->second );
  } else {
      plumed_merror("could not interpret switching function definition");
  }
  if( components ) matinp.push_back("COMPONENTS");
  actions.push_back( matinp );
}

void SymmetryFunctionBase::registerKeywords( Keywords& keys ){
  Action::registerKeywords( keys );
  ActionWithValue::registerKeywords( keys );
  ActionWithArguments::registerKeywords( keys );
  keys.add("compulsory","WEIGHT","");
  keys.add("numbered","VECTORS","");
}

SymmetryFunctionBase::SymmetryFunctionBase(const ActionOptions&ao):
Action(ao),
ActionWithValue(ao),
ActionWithArguments(ao)
{
  std::vector<std::string> alabels(1); std::vector<Value*> wval; parseArgumentList("WEIGHT",wval);
  if( wval.size()!=1 ) error("keyword WEIGHT should be provided with the label of a single action"); 
  alabels[0]=(wval[0]->getPntrToAction())->getLabel(); (wval[0]->getPntrToAction())->addActionToChain( alabels, this );
  log.printf("  using bond weights from %s \n",wval[0]->getName().c_str() );
  nderivatives=(wval[0]->getPntrToAction())->getNumberOfDerivatives(); 

  if( keywords.exists("VECTORS") ) {
      for(unsigned i=1;i<=3;++i){
          std::vector<Value*> vecs; parseArgumentList("VECTORS",i,vecs);
          if( vecs.size()!=1 ) error("keywords VECTORS should be provided with the label of a single action");
          if( wval[0]->getRank()!=vecs[0]->getRank() ) error("rank of weights does not match rank of vector");
          if( wval[0]->getRank()==2 ){
              if( wval[0]->getShape()[0]!=vecs[0]->getShape()[0] || wval[0]->getShape()[1]!=vecs[0]->getShape()[1] ){
                  error("mismatched shapes of matrices in input"); 
              }
          } else if( wval[0]->getRank()==1 && wval[0]->getShape()[0]!=vecs[0]->getShape()[0] ) error("mismatched shapes of vectors in input");
          if( (wval[0]->getPntrToAction())->getLabel()!=(vecs[0]->getPntrToAction())->getLabel() ){
               error("found mismatched vectors and weights in input to symmetry function - current not available, please email plumed list");
          }
          alabels[0]=(vecs[0]->getPntrToAction())->getLabel(); (vecs[0]->getPntrToAction())->addActionToChain( alabels, this ); wval.push_back(vecs[0]); 
          std::string dir="x"; if( i==2 ) dir="y"; else dir="z";
          log.printf("  %s direction of bond read from %s \n",dir.c_str(),vecs[0]->getName().c_str() );
      }
  }
  requestArguments(wval,true); forcesToApply.resize( nderivatives );
  if( plumed.getAtoms().getAllGroups().count(wval[0]->getPntrToAction()->getLabel()) ){
     const auto m=plumed.getAtoms().getAllGroups().find(wval[0]->getPntrToAction()->getLabel());
     plumed.getAtoms().insertGroup( getLabel(), m->second ); 
  } 
}

void SymmetryFunctionBase::addValueWithDerivatives() {
  std::vector<unsigned> shape;
  if( getPntrToArgument(0)->getRank()==2 ){
      shape.resize(1); shape[0]=getPntrToArgument(0)->getShape()[0];
  } 
  ActionWithValue::addValue( shape ); setNotPeriodic();
}

void SymmetryFunctionBase::addComponentWithDerivatives( const std::string& name ) { 
  std::vector<unsigned> shape;
  if( getPntrToArgument(0)->getRank()==2 ){
      shape.resize(1); shape[0]=getPntrToArgument(0)->getShape()[0];
  }
  ActionWithValue::addComponent(name,shape); componentIsNotPeriodic(name); 
}

void SymmetryFunctionBase::buildCurrentTaskList( std::vector<unsigned>& tflags ) {
  plumed_assert( actionInChain() ); tflags.assign(tflags.size(),1);
}

void SymmetryFunctionBase::performTask( const unsigned& current, MultiValue& myvals ) const {
  double weight = myvals.get( getPntrToArgument(0)->getPositionInStream() );
  if( weight>epsilon ){
      Vector dir; 
      dir[0] = myvals.get( getPntrToArgument(1)->getPositionInStream() ); 
      dir[1] = myvals.get( getPntrToArgument(2)->getPositionInStream() ); 
      dir[2] = myvals.get( getPntrToArgument(3)->getPositionInStream() );
      compute( weight, dir, myvals ); 
  }
  updateDerivativeIndices( myvals );
} 

void SymmetryFunctionBase::updateDerivativeIndices( MultiValue& myvals ) const {
  if( !doNotCalculateDerivatives() && myvals.inVectorCall() ) {
      // Update derivatives for indices
      unsigned istrn = getPntrToArgument(0)->getPositionInMatrixStash();
      std::vector<unsigned>& mat_indices( myvals.getMatrixIndices( istrn ) );
      for(unsigned i=0;i<myvals.getNumberOfMatrixIndices(istrn);++i) {
          for(unsigned j=0;j<getNumberOfComponents();++j){
              unsigned ostrn = getPntrToOutput(j)->getPositionInStream();
              myvals.updateIndex( ostrn, mat_indices[i] );
          }
      }
  }
} 

void SymmetryFunctionBase::apply() {
  if( doNotCalculateDerivatives() ) return;
  std::fill(forcesToApply.begin(),forcesToApply.end(),0); unsigned mm=0;
  if( getForcesFromValues( forcesToApply ) ) setForcesOnArguments( forcesToApply, mm ); 
}

}
}

