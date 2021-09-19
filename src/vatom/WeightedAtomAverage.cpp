/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2019 The plumed team
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
#include "WeightedAtomAverage.h"
#include "core/ActionWithArguments.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "core/Atoms.h"
#include <cmath>

using namespace std;

namespace PLMD {
namespace vatom {

void WeightedAtomAverage::registerKeywords(Keywords& keys) {
  ActionWithVirtualAtom::registerKeywords(keys);
  ActionWithValue::registerKeywords(keys); keys.remove("NUMERICAL_DERIVATIVES");
  keys.add("optional","WEIGHTS","what weights should be used when calculating the center.  If this keyword is not present the geometric center is computed. "
           "If WEIGHTS=@masses is used the center of mass is computed.  If WEIGHTS=@charges the center of charge is computed.  If "
           "the label of an action is provided PLUMED assumes that that action calculates a list of symmetry functions that can be used "
           "as weights. Lastly, an explicit list of numbers to use as weights can be provided");
  keys.addFlag("MASS",false,"calculate the center of mass");
  keys.addFlag("UNORMALIZED",false,"do not divide by the sum of the weights");
}

WeightedAtomAverage::WeightedAtomAverage(const ActionOptions&ao):
  Action(ao),
  ActionWithVirtualAtom(ao),
  myx(0), myw(0), nspace(1), bufstart(0),
  weight_mass(false),
  weight_charge(false),
  first(true),
  unorm(false),
  val_weights(NULL)
{
  vector<AtomNumber> atoms, catom; parseFlag("UNORMALIZED",unorm);
  parseAtomList("ATOMS",atoms); bool usemass = false; parseFlag("MASS",usemass);
  if(atoms.size()==0) error("at least one atom should be specified");
  if( keywords.exists("CENTER") ) {
      parseAtomList("CENTER",catom);
      if(catom.size()!=1) error("should be one central atom only");
      log.printf("  computing the gyration tensor around atom %d \n", catom[0].serial() );
  }
  std::vector<std::string> str_weights; parseVector("WEIGHTS",str_weights);
  if( usemass ) {
      if( str_weights.size()>0 ) error("USEMASS is incompatible with WEIGHTS");
      str_weights.resize(1); str_weights[0]="@masses";
  }
  std::string norm_str=""; if(unorm) norm_str = " unormalized";
  if( str_weights.size()==0 ) {
    if( catom.size()==0 ) log<<"  computing" + norm_str + " the geometric center of atoms:\n";
    else log<<"  computing" + norm_str + " gyration tensor\n";
    weights.resize( atoms.size() );
    for(unsigned i=0; i<atoms.size(); i++) weights[i] = 1.;
  } else if( str_weights.size()==1 ) {
    if( str_weights[0]=="@masses" ) {
      weight_mass=true;
      if( catom.size()==0 ) log<<"  computing" + norm_str + " the center of mass of atoms:\n";
      else log<<"  computing" + norm_str + " moment of inertia tensor\n";
    } else if( str_weights[0]=="@charges" ) {
      weight_charge=true;
      if( catom.size()==0 ) log<<"  computing" + norm_str + " the center of charge of atoms:\n";
      else log<<"  computing" + norm_str + " moment of charge tensor\n";
    } else {
      std::size_t dot=str_weights[0].find_first_of("."); unsigned nargs=0; std::vector<Value*> args;
      if( dot!=std::string::npos ) {
        ActionWithValue* action=plumed.getActionSet().selectWithLabel<ActionWithValue*>( str_weights[0].substr(0,dot) );
        if( !action ) {
          std::string str=" (hint! the actions in this ActionSet are: ";
          str+=plumed.getActionSet().getLabelList<ActionWithValue*>()+")";
          error("cannot find action named " + str_weights[0] +str);
        }
        action->copyOutput(str_weights[0].substr(dot+1))->use( this, args );
      } else {
        ActionWithValue* action=plumed.getActionSet().selectWithLabel<ActionWithValue*>( str_weights[0] );
        if( !action ) {
          std::string str=" (hint! the actions in this ActionSet are: ";
          str+=plumed.getActionSet().getLabelList<ActionWithValue*>()+")";
          error("cannot find action named " + str_weights[0] +str);
        }
        if( action->getNumberOfComponents()>1 ) error("requesting value from action " + action->getLabel() + " but action has components");
        action->copyOutput(0)->use( this, args );
      }
      if( args.size()!=1 ) error("should only have one value as input to WEIGHT");
      if( args[0]->getRank()!=1 || args[0]->getShape()[0]!=atoms.size() ) error("value input for WEIGHTS has wrong shape");
      val_weights = args[0]; std::vector<std::string> empty(1); empty[0] = (val_weights->getPntrToAction())->getLabel();
      if( catom.size()==0 && (val_weights->getPntrToAction())->valuesComputedInChain() ) (val_weights->getPntrToAction())->addActionToChain( empty, this );
      else val_weights->buildDataStore( getLabel() );
      log.printf("  atoms are weighted by values in vector labelled %s \n",val_weights->getName().c_str() );
      if( unorm ) log.printf("  final sum is not divided by sum of weights \n");
      else log<<"  final sum is divided by sum of weights \n";
    }
  } else {
    log<<" with weights:";
    if( str_weights.size()!=atoms.size() ) error("number of elements in weight vector does not match the number of atoms");
    weights.resize( atoms.size() );
    for(unsigned i=0; i<weights.size(); ++i) {
      if(i%25==0) log<<"\n";
      Tools::convert( str_weights[i], weights[i] ); log.printf(" %f",weights[i]);
    }
    log.printf("\n");
    if( unorm ) log.printf("  final sum is not divided by sum of weights \n");
    else log<<"  final sum is divided by sum of weights \n";
  }
  log.printf("  of atoms:");
  for(unsigned i=0; i<atoms.size(); ++i) {
    if(i>0 && i%25==0) log<<"\n";
    log.printf("  %d",atoms[i].serial());
  }
  log<<"\n";
  for(unsigned i=0;i<atoms.size();++i) addTaskToList(i);
  if( catom.size()>0 ) atoms.push_back( catom[0] );
  requestAtoms(atoms); if( val_weights ) addDependency( val_weights->getPntrToAction() ); 
}

unsigned WeightedAtomAverage::getNumberOfDerivatives() const {
  if( val_weights && actionInChain() ) return 3*getNumberOfAtoms() + (val_weights->getPntrToAction())->getNumberOfDerivatives(); 
  if( val_weights ) return 3*getNumberOfAtoms() + val_weights->getNumberOfValues();
  return 3*getNumberOfAtoms(); 
} 

void WeightedAtomAverage::setStashIndices( unsigned& nquants ) {
  myx = nquants; myw = nquants + getNumberOfStoredQuantities(); 
  nquants += getNumberOfStoredQuantities() + 1;
} 

void WeightedAtomAverage::getSizeOfBuffer( const unsigned& nactive_tasks, unsigned& bufsize ) {
  bufstart = bufsize; if( !doNotCalculateDerivatives() ) nspace = 1 + getNumberOfDerivatives();
  unsigned ntmp_vals = getNumberOfStoredQuantities(); bufsize += (ntmp_vals + 1)*nspace;
  if( final_vals.size()!=ntmp_vals ) {
      final_vals.resize( ntmp_vals ); weight_deriv.resize( getNumberOfDerivatives() );
      final_deriv.resize( ntmp_vals ); for(unsigned i=0;i<ntmp_vals;++i) final_deriv[i].resize( getNumberOfDerivatives() );
      if( val_weights ) {
          val_deriv.resize( 9 ); unsigned nder = val_weights->getNumberOfValues(); 
          if( actionInChain() ) nder = (val_weights->getPntrToAction())->getNumberOfDerivatives(); 
          val_forces.resize( nder ); for(unsigned i=0;i<9;++i) val_deriv[i].resize( nder );
      }
  }
  ActionWithValue::getSizeOfBuffer( nactive_tasks, bufsize );
} 

void WeightedAtomAverage::prepareForTasks( const unsigned& nactive, const std::vector<unsigned>& pTaskList ) {
  // Check that we have the mass information we need on the first step
  if( first ) {
    // Check if we have masses if we need them
    if( weight_mass ) {
        for(unsigned i=0; i<getNumberOfAtoms(); i++) {
          if(std::isnan(getMass(i))) {
            error(
              "You are trying to compute a CENTER or COM but masses are not known.\n"
              "        If you are using plumed driver, please use the --mc option"
            );
          }
        }
    // Check we have charges if we need them
    } else if ( weight_charge  ) {
        error(
            "You are trying to compute a center of charnge but chargest are not known.\n"
            "        If you are using plumed driver, please use the --mc option"
          );
    }
    first=false;
  }
  setupEntity();
}

void WeightedAtomAverage::calculate() {
  if( actionInChain() ) return;
  runAllTasks();
}

void WeightedAtomAverage::performTask( const unsigned& task_index, MultiValue& myvals ) const {
  Vector pos = getPosition( task_index ); double w;
  if( weight_mass ) {
    w = getMass(task_index);
  } else if( weight_charge ) {
    w = getCharge(task_index);
  } else if( val_weights && actionInChain() ) {
    w = myvals.get( val_weights->getPositionInStream() );
  } else if( val_weights ) {
    w = val_weights->get(task_index);
  } else {
    plumed_dbg_assert( task_index<weights.size() );
    w = weights[task_index];
  }
  unsigned ntmp_vals = getNumberOfStoredQuantities();
  myvals.addValue( myw, w ); compute( task_index, w, pos, myvals );
  if( !doNotCalculateDerivatives() && val_weights && fabs(w)>epsilon ) {
      double invw = 1 / w; unsigned base = 3*getNumberOfAtoms();
      if( actionInChain() ) {
          unsigned istrn = val_weights->getPositionInStream();
          for(unsigned k=0; k<myvals.getNumberActive(istrn); ++k) {
              unsigned kindex = myvals.getActiveIndex(istrn,k);
              double der = myvals.getDerivative( istrn, kindex );
              for(unsigned j=0; j<ntmp_vals;++j) addDerivative( j, base+kindex, der*invw*myvals.get(myx+j), myvals );
              myvals.addDerivative( myw, base+kindex, der ); myvals.updateIndex( myw, base+kindex );
          } 
      } else {
          for(unsigned j=0; j<ntmp_vals;++j) addDerivative( j, base+task_index, invw*myvals.get(myx+j), myvals );
          myvals.addDerivative( myw, base+task_index, 1.0 ); myvals.updateIndex( myw, base+task_index );
      }
  }
}

void WeightedAtomAverage::gatherForVirtualAtom( const MultiValue& myvals, std::vector<double>& buffer ) const { 
  unsigned ntmp_vals = getNumberOfStoredQuantities() + 1, bstart = bufstart;
  for(unsigned i=0;i<ntmp_vals;++i) { buffer[bstart] += myvals.get(myx+i); bstart += nspace; }

  if( !doNotCalculateDerivatives() ) {
      bstart = bufstart; 
      for(unsigned i=0;i<ntmp_vals;++i) {
          for(unsigned k=0; k<myvals.getNumberActive(myx+i); ++k) {
              unsigned kindex = myvals.getActiveIndex(myx+i,k);
              plumed_dbg_assert( bstart + 1 + kindex<buffer.size() );
              buffer[bstart + 1 + kindex] += myvals.getDerivative( myx+i, kindex );
          }
          bstart += nspace;
      }
  }
}

void WeightedAtomAverage::transformFinalValueAndDerivatives( const std::vector<double>& buffer ) {
  unsigned ntmp_vals = getNumberOfStoredQuantities(); 
  double ww = buffer[bufstart + ntmp_vals*nspace]; if( unorm ) ww = 1.0;
  // This finalizes the value 
  for(unsigned i=0;i<ntmp_vals;++i) final_vals[i] = buffer[bufstart+i*nspace]/ww;
  finalizeValue( final_vals );
  if( !doNotCalculateDerivatives() ) {
      for(unsigned i=0; i<getNumberOfDerivatives(); ++ i ) {
          for(unsigned j=0; j<ntmp_vals; ++j) {
              final_deriv[j][i] = buffer[bufstart + j*nspace + 1 + i] / ww;
          }
          weight_deriv[i] = buffer[bufstart + ntmp_vals*nspace + 1 + i ] / ww;
      }
      finalizeDerivatives( final_vals, final_deriv, weight_deriv, val_deriv );
  }
}

void WeightedAtomAverage::applyForcesToValue( const std::vector<double>& fff ) {
  val_forces.assign( val_forces.size(), 0.0 );
  for(unsigned j=0;j<fff.size();++j) {
      for(unsigned k=0;k<val_deriv[j].size();++k ) val_forces[k] += fff[j]*val_deriv[j][k]; 
  }
  if( actionInChain() ) {
      unsigned start=0; ActionWithArguments::setForcesOnActionChain( val_forces, start, val_weights->getPntrToAction() );
  } else {
      unsigned narg_v = val_weights->getNumberOfValues();
      for(unsigned j=0; j<narg_v; ++j) val_weights->addForce( j, val_forces[j] ); 
  }
}

}
}
