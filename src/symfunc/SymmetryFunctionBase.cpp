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
#include "core/ActionShortcut.h"
#include "core/PlumedMain.h"
#include "core/ActionSetup.h"
#include "core/Atoms.h"
#include "multicolvar/MultiColvarBase.h"

namespace PLMD {
namespace symfunc {

void SymmetryFunctionBase::shortcutKeywords( Keywords& keys ) {
  ActionShortcut::registerKeywords( keys );
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

void SymmetryFunctionBase::expandMatrix( const bool& components, const std::string& lab, const std::string& sp_str, 
                                         const std::string& spa_str, const std::string& spb_str, ActionShortcut* action ) {
  if( sp_str.length()==0 && spa_str.length()==0 ) return;

  std::string matinp = lab  + "_mat: CONTACT_MATRIX";
  if( sp_str.length()>0 ) matinp += " GROUP=" + sp_str;
  else if( spa_str.length()>0 ) matinp += " GROUPA=" + spa_str + " GROUPB=" + spb_str;

  std::string sw_str; action->parse("SWITCH",sw_str); 
  if( sw_str.length()>0 ) {
      matinp += " SWITCH={" + sw_str + "}";
  } else {
      std::string r0; action->parse("R_0",r0); std::string d0; action->parse("D_0",d0);
      if( r0.length()==0 ) action->error("missing switching function parameters use SWITCH/R_0"); 
      std::string nn; action->parse("NN",nn); std::string mm; action->parse("MM",mm);
      matinp += " R_0=" + r0 + " D_0=" + d0 + " NN=" + nn + " MM=" + mm;
  } 
  if( components ) matinp += " COMPONENTS";
  action->readInputLine( matinp );
} 

void SymmetryFunctionBase::createSymmetryFunctionObject( const std::string& lab, const std::string& name, const bool& iscoord, const bool& norm, ActionShortcut* action ) { 
  // Read species keywords and create matrix
  std::string sp_str, specA, specB; action->parse("SPECIES",sp_str); action->parse("SPECIESA",specA); action->parse("SPECIESB",specB);
  std::map<std::string,std::string> keymap; multicolvar::MultiColvarBase::readShortcutKeywords( keymap, action );
  if( sp_str.length()==0 && specA.length()==0 ) {
     action->readInputLine( lab + ": " + name + "_MATINP " + action->convertInputLineToString() ); 
  } else {
     SymmetryFunctionBase::expandMatrix( true, lab,  sp_str, specA, specB, action );
     // Create input for symmetry function
     if( iscoord ) {
        action->readInputLine( lab + ": " + name + " WEIGHT=" + lab + "_mat.w " + action->convertInputLineToString() ); 
     } else {
        action->readInputLine( lab + ": " + name + " WEIGHT=" + lab + "_mat.w VECTORS1=" + lab + "_mat.x VECTORS2=" + lab + "_mat.y VECTORS3=" + lab + "_mat.z " + 
                       action->convertInputLineToString() );
     }
     std::string olab = lab;
     if( norm ) { 
         olab = lab + "_n"; action->readInputLine( lab + "_denom: COORDINATIONNUMBER WEIGHT=" + lab + "_mat.w");
         // Input for matheval action
         action->readInputLine( lab + "_n: MATHEVAL ARG1=" + lab + " ARG2=" + lab + "_denom FUNC=x/y PERIODIC=NO");
     }
     multicolvar::MultiColvarBase::expandFunctions( lab, olab, "", keymap, action );
  }
}

void SymmetryFunctionBase::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys );
  ActionWithValue::registerKeywords( keys );
  ActionWithArguments::registerKeywords( keys );
  keys.add("compulsory","WEIGHT","");
  keys.add("numbered","VECTORS","");
  keys.addFlag("ONESHOT",false,"This forces all the elements of the row of the matrix to be computed prior to computing the symmetry function.  "
               "It should only be ever need to be used for testing.");
  keys.addFlag("USECOLS",false,"When this flag is present the CVs are calculated by summing over the columns rather than the rows.  You are thus calculating "
               "symmetry functions for the atoms in GROUPB rather than symmetry functions for the atoms in GROUPA.  The derivatives are much "
               "more expensive when this approach is used");
}

SymmetryFunctionBase::SymmetryFunctionBase(const ActionOptions&ao):
  Action(ao),
  ActionWithValue(ao),
  ActionWithArguments(ao),
  usecols(false),
  nderivatives(0),
  done_with_matrix_comput(true)
{
  if( keywords.exists("USECOLS") ) {
    parseFlag("USECOLS",usecols);
    if( usecols ) log.printf("  calculating symmetry functions for second group \n");
  }
  std::vector<std::string> alabels(1); std::vector<Value*> wval; parseArgumentList("WEIGHT",wval);
  if( wval.size()!=1 ) error("keyword WEIGHT should be provided with the label of a single action");
  if( wval[0]->getPntrToAction() ) {
      alabels[0]=(wval[0]->getPntrToAction())->getLabel();
      ActionSetup* as = dynamic_cast<ActionSetup*>( wval[0]->getPntrToAction() );
      if( !as ) (wval[0]->getPntrToAction())->addActionToChain( alabels, this );
      nderivatives=(wval[0]->getPntrToAction())->getNumberOfDerivatives();
  }
  log.printf("  using bond weights from matrix labelled %s \n",wval[0]->getName().c_str() );

  if( keywords.exists("VECTORS") ) {
    if( !wval[0]->getPntrToAction() ) error("using weights from input matrix not available with vectors");
    for(unsigned i=1; i<=3; ++i) {
      std::vector<Value*> vecs; parseArgumentList("VECTORS",i,vecs);
      if( vecs.size()!=1 ) error("keywords VECTORS should be provided with the label of a single action");
      if( wval[0]->getRank()!=vecs[0]->getRank() ) error("rank of weights does not match rank of vector");
      if( wval[0]->getRank()==2 ) {
        if( wval[0]->getShape()[0]!=vecs[0]->getShape()[0] || wval[0]->getShape()[1]!=vecs[0]->getShape()[1] ) {
          error("mismatched shapes of matrices in input");
        }
      } else if( wval[0]->getRank()==1 && wval[0]->getShape()[0]!=vecs[0]->getShape()[0] ) error("mismatched shapes of vectors in input");
      // This checks if the weights come from a different action than the vectors
      bool found=false;
      for(unsigned j=0;j<wval.size();++j) {
          if( wval[j]->getPntrToAction()->getLabel()==vecs[0]->getPntrToAction()->getLabel() ) { found=true; break; }
      }
      if( !found ) nderivatives += (vecs[0]->getPntrToAction())->getNumberOfDerivatives();

      if( ((wval[0]->getPntrToAction())->getActionThatCalculates())->getLabel()!=((vecs[0]->getPntrToAction())->getActionThatCalculates())->getLabel() ) {
        error("found mismatched vectors and weights in input to symmetry function (2nd version) - current not available, please email plumed list");
      }
      alabels[0]=(vecs[0]->getPntrToAction())->getLabel(); ActionSetup* as2 = dynamic_cast<ActionSetup*>( vecs[0]->getPntrToAction() );
      if( as2 ) (vecs[0]->getPntrToAction())->addActionToChain( alabels, this );
      wval.push_back(vecs[0]); std::string dir="x"; if( i==2 ) dir="y"; else dir="z";
      log.printf("  %s direction of bond read from matrix labelled %s \n",dir.c_str(),vecs[0]->getName().c_str() );
    }
  }
  if( keywords.exists("ONESHOT") ) {
    bool oneshot; parseFlag("ONESHOT",oneshot);
    if( oneshot ) {
      done_with_matrix_comput=false;
      log.printf("  computing full matrix rows before computing symmetry function \n");
    }
  } else {
    done_with_matrix_comput=false;
  }
  // If we are doing the calculation as we compute matrix elements we must store all matrix elements
  // in rows.  Actually we store whole matrix because I don't want to make more complicated options
  if( !done_with_matrix_comput ) {
    // The -mat here is added to prevent this behaving like a proper stored value when updating forces
    for(unsigned i=0; i<wval.size(); ++i) wval[i]->buildDataStore( getLabel() + "-mat" );
  }
  requestArguments(wval,true); forcesToApply.resize( nderivatives );
  if( getPntrToArgument(0)->getRank()==2 ) {
    for(unsigned i=0; i<getPntrToArgument(0)->getShape()[0]; ++i) addTaskToList(i);
  }
  if( !usecols && wval[0]->getPntrToAction() ) {
    if( plumed.getAtoms().getAllGroups().count(wval[0]->getPntrToAction()->getLabel()) ) {
        const auto m=plumed.getAtoms().getAllGroups().find(wval[0]->getPntrToAction()->getLabel());
        plumed.getAtoms().insertGroup( getLabel(), m->second );
    }
  }
}

void SymmetryFunctionBase::interpretDotStar( const std::string& ulab, unsigned& nargs, std::vector<Value*>& myvals ) {
  multicolvar::MultiColvarBase::interpretDotStar( getLabel(), ulab, nargs, myvals, plumed.getActionSet() );
}

void SymmetryFunctionBase::addValueWithDerivatives() {
  std::vector<unsigned> shape;
  if( getPntrToArgument(0)->getRank()==2 ) {
    shape.resize(1);
    if( usecols ) shape[0]=getPntrToArgument(0)->getShape()[1];
    else shape[0]=getPntrToArgument(0)->getShape()[0];
    if( shape[0]==1 ) shape.resize(0);
  }
  if( shape.size()==0 ) ActionWithValue::addValueWithDerivatives( shape );
  else ActionWithValue::addValue( shape );
  setNotPeriodic();
  if( usecols ) getPntrToOutput( getNumberOfComponents()-1 )->buildColumnSums();
}

void SymmetryFunctionBase::addComponentWithDerivatives( const std::string& name ) {
  std::vector<unsigned> shape;
  if( getPntrToArgument(0)->getRank()==2 ) {
    shape.resize(1);
    if( usecols ) shape[0]=getPntrToArgument(0)->getShape()[1];
    else shape[0]=getPntrToArgument(0)->getShape()[0];
    if( shape[0]==1 ) shape.resize(0);
  }
  if( shape.size()==0 ) ActionWithValue::addComponentWithDerivatives(name,shape);
  else ActionWithValue::addComponent(name,shape);
  componentIsNotPeriodic(name);
  if( usecols ) getPntrToOutput( getNumberOfComponents()-1 )->buildColumnSums();
}

void SymmetryFunctionBase::performTask( const unsigned& current, MultiValue& myvals ) const {
  if( !myvals.inVectorCall() && done_with_matrix_comput && !myvals.inMatrixRerun() ) {
    double weight = myvals.get( getPntrToArgument(0)->getPositionInStream() );
    if( fabs(weight)>epsilon ) {
      Vector dir; dir.zero();
      if( getNumberOfArguments()==4 ) {
        dir[0] = myvals.get( getPntrToArgument(1)->getPositionInStream() );
        dir[1] = myvals.get( getPntrToArgument(2)->getPositionInStream() );
        dir[2] = myvals.get( getPntrToArgument(3)->getPositionInStream() );
      }
      compute( weight, dir, myvals );
    }
  } else if( myvals.inVectorCall() ) {
    if( !done_with_matrix_comput ) {
      // Make sure tempory derivative space is set up if required
      if( !doNotCalculateDerivatives() ) {
        std::vector<double>& tmp_w( myvals.getSymfuncTemporyDerivatives( getPntrToArgument(0)->getPositionInStream() ) );
        if( tmp_w.size()<getNumberOfComponents()*getPntrToArgument(0)->getShape()[1] ) tmp_w.resize( getNumberOfComponents()*getPntrToArgument(0)->getShape()[1], 0 );
        std::vector<double>& tmp_x( myvals.getSymfuncTemporyDerivatives( getPntrToArgument(1)->getPositionInStream() ) );
        if( tmp_x.size()<getNumberOfComponents()*getPntrToArgument(1)->getShape()[1] ) tmp_x.resize( getNumberOfComponents()*getPntrToArgument(1)->getShape()[1], 0 );
        std::vector<double>& tmp_y( myvals.getSymfuncTemporyDerivatives( getPntrToArgument(2)->getPositionInStream() ) );
        if( tmp_y.size()<getNumberOfComponents()*getPntrToArgument(2)->getShape()[1] ) tmp_y.resize( getNumberOfComponents()*getPntrToArgument(2)->getShape()[1], 0 );
        std::vector<double>& tmp_z( myvals.getSymfuncTemporyDerivatives( getPntrToArgument(3)->getPositionInStream() ) );
        if( tmp_z.size()<getNumberOfComponents()*getPntrToArgument(3)->getShape()[1] ) tmp_z.resize( getNumberOfComponents()*getPntrToArgument(3)->getShape()[1], 0 );
      }
      computeSymmetryFunction( current, myvals );
      // And now the derivatives
      if( !doNotCalculateDerivatives() ) {
        ActionWithValue* av = (getPntrToArgument(0)->getPntrToAction())->getActionThatCalculates();
        unsigned aindex_start = myvals.getNumberOfIndicesInFirstBlock();
        unsigned matind = getPntrToArgument(0)->getPositionInMatrixStash();
        unsigned my_w = getPntrToArgument(0)->getPositionInStream();
        std::vector<double>& tmp_w( myvals.getSymfuncTemporyDerivatives(my_w) );
        // Turn off matrix element storing during rerun of calculations
        myvals.setMatrixStashForRerun();
        if( getNumberOfArguments()==4 ) {
          unsigned my_x = getPntrToArgument(1)->getPositionInStream();
          std::vector<double>& tmp_x( myvals.getSymfuncTemporyDerivatives(my_x) );
          unsigned my_y = getPntrToArgument(2)->getPositionInStream();
          std::vector<double>& tmp_y( myvals.getSymfuncTemporyDerivatives(my_y) );
          unsigned my_z = getPntrToArgument(3)->getPositionInStream();
          std::vector<double>& tmp_z( myvals.getSymfuncTemporyDerivatives(my_z) );
          for(unsigned j=0; j<myvals.getNumberOfStashedMatrixElements(matind); ++j) {
            unsigned wstart=0, jind = myvals.getStashedMatrixIndex(matind,j);
            // Check for derivatives and skip recalculation if there are none
            double totder=0;
            for(unsigned i=0; i<getNumberOfComponents(); ++i) {
              totder +=tmp_w[wstart+jind] + tmp_x[wstart+jind] + tmp_y[wstart+jind] + tmp_z[wstart+jind];
              wstart++;
            }
            if( fabs(totder)<epsilon ) continue ;
            // Rerun the task required
            wstart = 0; av->runTask( av->getLabel(), myvals.getTaskIndex(), current, aindex_start + jind, myvals );
            // Now add on the derivatives
            for(unsigned i=0; i<getNumberOfComponents(); ++i) {
              unsigned ostrn = getPntrToOutput(i)->getPositionInStream();
              for(unsigned k=0; k<myvals.getNumberActive(my_w); ++k) {
                unsigned kind=myvals.getActiveIndex(my_w,k);
                myvals.addDerivative( ostrn, arg_deriv_starts[0] + kind, tmp_w[wstart+jind]*myvals.getDerivative( my_w, kind ) );
              }
              for(unsigned k=0; k<myvals.getNumberActive(my_x); ++k) {
                unsigned kind=myvals.getActiveIndex(my_x,k);
                myvals.addDerivative( ostrn, arg_deriv_starts[1] + kind, tmp_x[wstart+jind]*myvals.getDerivative( my_x, kind ) );
              }
              for(unsigned k=0; k<myvals.getNumberActive(my_y); ++k) {
                unsigned kind=myvals.getActiveIndex(my_y,k);
                myvals.addDerivative( ostrn, arg_deriv_starts[2] + kind, tmp_y[wstart+jind]*myvals.getDerivative( my_y, kind ) );
              }
              for(unsigned k=0; k<myvals.getNumberActive(my_z); ++k) {
                unsigned kind=myvals.getActiveIndex(my_z,k);
                myvals.addDerivative( ostrn, arg_deriv_starts[3] + kind, tmp_z[wstart+jind]*myvals.getDerivative( my_z, kind ) );
              }
              tmp_w[wstart+jind]=tmp_x[wstart+jind]=tmp_y[wstart+jind]=tmp_z[wstart+jind]=0;
              wstart += getPntrToArgument(0)->getShape()[1];
            }
            // Clear the matrix elements for this task
            av->clearMatrixElements( myvals );
          }
        } else {
          for(unsigned j=0; j<myvals.getNumberOfStashedMatrixElements(matind); ++j) {
            unsigned jind = myvals.getStashedMatrixIndex(matind,j);
            // Rerun the task required
            av->runTask( av->getLabel(), myvals.getTaskIndex(), current, aindex_start + jind, myvals );
            // Now add on the derivatives
            for(unsigned i=0; i<getNumberOfComponents(); ++i) {
              unsigned wstart=0, ostrn = getPntrToOutput(i)->getPositionInStream();
              for(unsigned k=0; k<myvals.getNumberActive(my_w); ++k) {
                unsigned kind=myvals.getActiveIndex(my_w,k);
                myvals.addDerivative( ostrn, arg_deriv_starts[i] + kind, tmp_w[wstart+jind]*myvals.getDerivative( my_w, kind ) );
              }
              tmp_w[wstart+jind]=0; wstart += getPntrToArgument(0)->getShape()[1];
            }
            // Clear the matrix elements for this task
            av->clearMatrixElements( myvals );
          }
        }
        // Set the myvals object to store matrix elements
        myvals.setMatrixStashForNormalRun();
      }
    }
    // Update derivatives for indices
    if( !doNotCalculateDerivatives() ) updateDerivativeIndices( myvals );
  }
}

void SymmetryFunctionBase::updateDerivativeIndices( MultiValue& myvals ) const {
  unsigned istrn = getPntrToArgument(0)->getPositionInMatrixStash();
  std::vector<unsigned>& mat_indices( myvals.getMatrixIndices( istrn ) );
  for(unsigned i=0; i<myvals.getNumberOfMatrixIndices(istrn); ++i) {
    for(unsigned j=0; j<getNumberOfComponents(); ++j) {
      unsigned ostrn = getPntrToOutput(j)->getPositionInStream();
      myvals.updateIndex( ostrn, mat_indices[i] );
    }
  }
  if( getNumberOfArguments()>1 ) {
    if( arg_deriv_starts[1]>0 && getPntrToArgument(1)->getRank()==2 ) {
       istrn = getPntrToArgument(1)->getPositionInMatrixStash();
       std::vector<unsigned>& mat_indices( myvals.getMatrixIndices( istrn ) );
       for(unsigned i=0; i<myvals.getNumberOfMatrixIndices(istrn); ++i) {
         for(unsigned j=0; j<getNumberOfComponents(); ++j) {
           unsigned ostrn = getPntrToOutput(j)->getPositionInStream();
           myvals.updateIndex( ostrn, arg_deriv_starts[1] + mat_indices[i] );
         }
       }   
    }
    // It would be easy enough to extend to allow vector components to come from different actions.  This will 
    // likely be of no use to anyone though so for the moment I just have this assert to prevent people from doing
    // wrong things.
    plumed_dbg_assert( arg_deriv_starts[2]==arg_deriv_starts[1] && arg_deriv_starts[3]==arg_deriv_starts[1] ); 
  }
}

void SymmetryFunctionBase::computeSymmetryFunction( const unsigned& current, MultiValue& myvals ) const {
  Vector dir;
  unsigned matind = getPntrToArgument(0)->getPositionInMatrixStash();
  if( getNumberOfArguments()>1 ) {
    unsigned matind_x = getPntrToArgument(1)->getPositionInMatrixStash();
    unsigned matind_y = getPntrToArgument(2)->getPositionInMatrixStash();
    unsigned matind_z = getPntrToArgument(3)->getPositionInMatrixStash();
    for(unsigned j=0; j<myvals.getNumberOfStashedMatrixElements(matind); ++j) {
      unsigned jind = myvals.getStashedMatrixIndex(matind,j);
      double weight = myvals.getStashedMatrixElement( matind, jind );
      dir[0] = myvals.getStashedMatrixElement( matind_x, jind );
      dir[1] = myvals.getStashedMatrixElement( matind_y, jind );
      dir[2] = myvals.getStashedMatrixElement( matind_z, jind );
      myvals.setSymfuncTemporyIndex(jind); compute( weight, dir, myvals );
    }
  } else {
    for(unsigned j=0; j<myvals.getNumberOfStashedMatrixElements(matind); ++j) {
      unsigned jind = myvals.getStashedMatrixIndex(matind,j);
      double weight = myvals.getStashedMatrixElement( matind, jind );
      myvals.setSymfuncTemporyIndex(jind); compute( weight, dir, myvals );
    }
  }
}

void SymmetryFunctionBase::apply() {
  if( doNotCalculateDerivatives() ) return;
  if( forcesToApply.size()!=getNumberOfDerivatives() ) forcesToApply.resize( getNumberOfDerivatives() );
  std::fill(forcesToApply.begin(),forcesToApply.end(),0); unsigned mm=0;
  if( getForcesFromValues( forcesToApply ) ) { setForcesOnArguments( 0, forcesToApply, mm ); }
}

}
}

