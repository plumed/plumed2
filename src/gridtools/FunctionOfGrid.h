/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2020 The plumed team
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
#ifndef __PLUMED_gridtools_FunctionOfGrid_h
#define __PLUMED_gridtools_FunctionOfGrid_h

#include "core/ActionWithValue.h"
#include "core/ActionWithArguments.h"
#include "function/Custom.h"
#include "tools/Matrix.h"

namespace PLMD {
namespace gridtools {

template <class T>
class FunctionOfGrid : 
public ActionWithValue,
public ActionWithArguments
{
private:
/// The function that is being computed
  T myfunc;
/// Is this the first step
  bool firststep;
/// The number of derivatives
  unsigned nderivatives;
/// Forces if we are doing the integral
  std::vector<double> forcesToApply;
/// Function for setting up the output
  void reshapeOutput();
/// Function for running a calculation on the whole grid
  void runSingleTaskCalculation( const Value* arg, ActionWithValue* action, T& f );
public:
  static void registerKeywords(Keywords&);
  explicit FunctionOfGrid(const ActionOptions&);
/// Get the number of derivatives for this action
  unsigned getNumberOfDerivatives() const override ;
/// Get the label to write in the graph
  std::string writeInGraph() const override { return myfunc.getGraphInfo( getName() ); }
  void calculate() override;
  void update() override;
  void runFinalJobs() override;
/// Calculate the function
  void performTask( const unsigned& current, MultiValue& myvals ) const override ;
///
  void gatherStoredValue( const unsigned& valindex, const unsigned& code, const MultiValue& myvals,
                          const unsigned& bufstart, std::vector<double>& buffer ) const override ;
/// Add the forces 
  void apply() override;
};

template <class T>
void FunctionOfGrid<T>::registerKeywords(Keywords& keys ) {
  Action::registerKeywords(keys); ActionWithValue::registerKeywords(keys); ActionWithArguments::registerKeywords(keys); keys.use("ARG");
  keys.reserve("compulsory","PERIODIC","if the output of your function is periodic then you should specify the periodicity of the function.  If the output is not periodic you must state this using PERIODIC=NO");
  T tfunc; tfunc.registerKeywords( keys ); if( typeid(tfunc)==typeid(function::Custom()) ) keys.add("hidden","NO_ACTION_LOG","suppresses printing from action on the log");
}

template <class T>
FunctionOfGrid<T>::FunctionOfGrid(const ActionOptions&ao):
Action(ao),
ActionWithValue(ao),
ActionWithArguments(ao),
firststep(true),
nderivatives(0)
{
  if( getNumberOfArguments()==0 ) error("found no arguments");
  bool foundgrid=false; std::vector<double> gspacing; std::vector<unsigned> nbin; std::vector<bool> pbc;
  unsigned argstart=myfunc.getArgStart(); unsigned npoints=0; 
  std::string gtype; std::vector<std::string> gargn, min, max; double volume;
  for(unsigned i=argstart; i<getNumberOfArguments(); ++i) {
    if( getPntrToArgument(i)->getRank()>0 && getPntrToArgument(i)->hasDerivatives() ) {
        foundgrid=true; npoints=getPntrToArgument(i)->getNumberOfValues();
        nderivatives = getPntrToArgument(i)->getRank() + getNumberOfArguments() - argstart;
        gspacing.resize( getPntrToArgument(i)->getRank() ); nbin.resize( getPntrToArgument(i)->getRank() );
        min.resize( getPntrToArgument(i)->getRank() ); max.resize( getPntrToArgument(i)->getRank() );
        gargn.resize( getPntrToArgument(i)->getRank() ); pbc.resize( getPntrToArgument(i)->getRank() );
        (getPntrToArgument(i)->getPntrToAction())->getInfoForGridHeader( gtype, gargn, min, max, nbin, gspacing, pbc, false );
        if( gtype=="flat" ) {
            volume=1; for(unsigned j=0;j<gspacing.size();++j) volume *= gspacing[j];
        } else volume=4*pi/ static_cast<double>( npoints );
        break;
    }
  }
  if( !foundgrid ) error("found no grid in input");
  // Now check that all grids have the same size
  std::vector<unsigned> shape( min.size() ); std::vector<unsigned> gnbin( min.size() ); std::vector<bool> gpbc( min.size() );
  std::vector<std::string> ggargn( min.size() ), gmin( min.size() ), gmax( min.size() ); std::string ggtype;
  if( arg_ends.size()==0 && getNumberOfArguments()==1 ) { arg_ends.push_back(0); arg_ends.push_back(1); }
  for(unsigned j=argstart; j<getNumberOfArguments(); ++j) {
    if( getPntrToArgument(j)->getRank()!=0 ) {
      if( getPntrToArgument(j)->getNumberOfValues()!=npoints || !getPntrToArgument(j)->hasDerivatives() ) error("mismatch in input arguments");
      (getPntrToArgument(j)->getPntrToAction())->getInfoForGridHeader( ggtype, ggargn, gmin, gmax, gnbin, gspacing, gpbc, false );
      if( gtype!=ggtype ) error("mismatch between grid types");
      for(unsigned k=0;k<min.size();++k) {
          if( min[k]!=gmin[k] ) error("mismatch between input grid domains");
          if( max[k]!=gmax[k] ) error("mismatch between input grid domains");
          if( pbc[k]!=gpbc[k] ) error("mismatch between input grid domains");
          if( nbin[k]!=gnbin[k] ) error("mismatch between input grid domains");
      }
      // Make sure elements of shape are set correctly
      std::vector<unsigned> ss( getPntrToArgument(j)->getShape() ); for(unsigned i=0;i<ss.size();++i) shape[i]=ss[i];   
    } 
  }
  // Read the input and do some checks
  myfunc.read( this );
  // Check we are not calculating an integral
  if( myfunc.zeroRank() ) { 
      shape.resize(0); forcesToApply.resize( npoints ); 
      for(unsigned j=0; j<npoints; ++j) addTaskToList(j);
  } 
  // This sets the prefactor to the volume which converts integrals to sums
  myfunc.setPrefactor( this, volume );
  // Check that derivatives are available 
  if( !myfunc.derivativesImplemented() ) error("derivatives have not been implemended for " + getName() );
  // Get the names of the components
  std::vector<std::string> components( keywords.getAllOutputComponents() );
  // Create the values to hold the output
  if( components.size()==0 && myfunc.zeroRank() ) addValueWithDerivatives();
  else if( components.size()==0 ) { 
     addValueWithDerivatives( shape ); getPntrToOutput(0)->alwaysStoreValues();
  } else error("functions of grid should only output one grid");
  // Set the periodicities of the output components
  myfunc.setPeriodicityForOutputs( this );
}

template <class T>
unsigned FunctionOfGrid<T>::getNumberOfDerivatives() const {
  if( myfunc.zeroRank() ) return getPntrToArgument(0)->getNumberOfValues();
  return nderivatives;
}

template <class T>
void FunctionOfGrid<T>::performTask( const unsigned& current, MultiValue& myvals ) const {
  unsigned argstart=myfunc.getArgStart(); std::vector<double> args( getNumberOfArguments() - argstart ); 
  for(unsigned i=argstart;i<getNumberOfArguments();++i) {
      if( getPntrToArgument(i)->getRank()==0 ) args[i-argstart]=getPntrToArgument(i)->get();
      else args[i-argstart] = getPntrToArgument(i)->get(current);
  }
  // Calculate the function and its derivatives
  std::vector<double> vals(1); Matrix<double> derivatives( 1, getNumberOfArguments()-argstart );
  myfunc.calc( this, args, vals, derivatives ); unsigned np = myvals.getTaskIndex();
  // And set the values and derivatives
  unsigned ostrn = getPntrToOutput(0)->getPositionInStream();
  myvals.addValue( ostrn, vals[0] ); 
  if( !myfunc.zeroRank() ) {  
      // Add the derivatives for a grid
      for(unsigned j=argstart;j<getNumberOfArguments();++j) {
          // We store all the derivatives of all the input values - i.e. the grid points these are used in apply
          myvals.addDerivative( ostrn, getPntrToOutput(0)->getRank()+j-argstart, derivatives(0,j-argstart) );
          // And now we calculate the derivatives of the value that is stored on the grid correctly so that we can interpolate functions 
          if( getPntrToArgument(j)->getRank()!=0 ) {
              for(unsigned k=0; k<getPntrToArgument(j)->getRank(); ++k) myvals.addDerivative( ostrn, k, derivatives(0,j-argstart)*getPntrToArgument(j)->getGridDerivative( np, k ) );
          }   
      }
      for(unsigned j=0; j<nderivatives; ++j) myvals.updateIndex( ostrn, j );
  } else if( !doNotCalculateDerivatives() ) { 
      // These are the derivatives of the integral
      myvals.addDerivative( ostrn, current, derivatives(0,0) ); myvals.updateIndex( ostrn, current ); 
  }
}

template <class T>
void FunctionOfGrid<T>::gatherStoredValue( const unsigned& valindex, const unsigned& code, const MultiValue& myvals,
                                           const unsigned& bufstart, std::vector<double>& buffer ) const {
  if( getPntrToOutput(0)->getRank()>0 && getPntrToOutput(0)->hasDerivatives() ) {
      plumed_dbg_assert( getNumberOfComponents()==1 && valindex==0 );
      unsigned nder = getPntrToOutput(0)->getRank(), ostr = getPntrToOutput(0)->getPositionInStream();
      unsigned kp = bufstart + code*(1+nderivatives); buffer[kp] += myvals.get( ostr );
      for(unsigned i=0; i<nderivatives; ++i) buffer[kp + 1 + i] += myvals.getDerivative( ostr, i );
  } else ActionWithValue::gatherStoredValue( valindex, code, myvals, bufstart, buffer );
}

template <class T>
void FunctionOfGrid<T>::reshapeOutput() {
  plumed_assert( firststep ); myfunc.setup( this ); 
  unsigned npoints=0, argstart=myfunc.getArgStart();
  for(unsigned i=argstart; i<getNumberOfArguments(); ++i) {
    if( getPntrToArgument(i)->getRank()>0 && getPntrToArgument(i)->hasDerivatives() ) { npoints=getPntrToArgument(i)->getNumberOfValues(); break; }
  }
  if( !myfunc.zeroRank() && getPntrToOutput(0)->getNumberOfValues()!=npoints ) {
      std::vector<unsigned> shape;
      for(unsigned i=argstart; i<getNumberOfArguments(); ++i) {
         if( getPntrToArgument(i)->getRank()>0 && getPntrToArgument(i)->hasDerivatives() ) {
             unsigned dim=getPntrToArgument(i)->getShape().size();
             shape.resize(dim); for(unsigned j=0;j<dim;++j) shape[j]=getPntrToArgument(i)->getShape()[j];
             break;
         }
      }
      plumed_assert( shape.size()>0 ); getPntrToOutput(0)->setShape( shape ); 
  }
  if( myfunc.doWithTasks() && getFullNumberOfTasks()<npoints ) {
      if( myfunc.zeroRank() ) { getPntrToOutput(0)->resizeDerivatives( npoints ); forcesToApply.resize( npoints ); } 
      for(unsigned j=getFullNumberOfTasks(); j<npoints; ++j) addTaskToList(j);
  }
  firststep=false;
}

template <class T>
void FunctionOfGrid<T>::runSingleTaskCalculation( const Value* arg, ActionWithValue* action, T& f ) {
  unsigned nv = arg->getNumberOfValues(); std::vector<double> args( nv );
  for(unsigned i=0;i<nv;++i) args[i] = arg->get(i);
  std::vector<double> vals( nv ); Matrix<double> derivatives( arg->getRank(), nv );
  ActionWithArguments* aa=dynamic_cast<ActionWithArguments*>(action); plumed_assert( aa ); f.calc( aa, args, vals, derivatives );
  for(unsigned i=0;i<vals.size();++i) action->copyOutput(0)->set( i, vals[i] );
  // Return if we are not computing derivatives
  if( action->doNotCalculateDerivatives() ) return;
  for(unsigned i=0; i<nv; ++i) {
      for(unsigned j=0; j<arg->getRank(); ++j) action->copyOutput(0)->setGridDerivative( i, j, derivatives(j,i) );
  }
}

template <class T>
void FunctionOfGrid<T>::calculate() {
  // This is done if we are calculating a function of multiple cvs
  if( !skipUpdate() ) return;
  if( firststep ) reshapeOutput();
  // This is done if we are calculating a function of multiple cvs
  if( getFullNumberOfTasks()>0 ) runAllTasks();
  // This is used if we are doing sorting actions on a single vector
  else if( !myfunc.doWithTasks() ) runSingleTaskCalculation( getPntrToArgument(0), this, myfunc );
}

template <class T>
void FunctionOfGrid<T>::update() {
  if( skipUpdate() ) return;
  if( firststep ) reshapeOutput();
  // This is done if we are calculating a function of multiple cvs
  if( getFullNumberOfTasks()>0 ) runAllTasks();
  // This is used if we are doing sorting actions on a single vector
  else if( !myfunc.doWithTasks() ) runSingleTaskCalculation( getPntrToArgument(0), this, myfunc );
}
  
template <class T>
void FunctionOfGrid<T>::runFinalJobs() {
  if( skipUpdate()  ) return;
  if( firststep ) reshapeOutput(); 
  // This is done if we are calculating a function of multiple cvs
  if( getFullNumberOfTasks()>0 ) runAllTasks();
  // This is used if we are doing sorting actions on a single vector
  else if( !myfunc.doWithTasks() ) runSingleTaskCalculation( getPntrToArgument(0), this, myfunc );
}

template <class T>
void FunctionOfGrid<T>::apply() {
  if( doNotCalculateDerivatives() ) return;
          
  if( !getPntrToOutput(0)->forcesWereAdded() ) return ;  

  // This applies forces for the integral
  if( myfunc.zeroRank() ) { 
      std::fill(forcesToApply.begin(),forcesToApply.end(),0); unsigned fstart=0;
      getForcesFromValues( forcesToApply ); setForcesOnArguments( 0, forcesToApply, fstart );
      return;
  }
 
  // Work out how to deal with arguments
  unsigned nscalars=0, argstart=myfunc.getArgStart();
  for(unsigned i=argstart; i<getNumberOfArguments(); ++i) {
    if( getPntrToArgument(i)->getRank()==0 ) { nscalars++; }
  }
  
  std::vector<double> totv(nscalars,0);
  for(unsigned i=0; i<getFullNumberOfTasks(); ++i) { 
    nscalars=0;
    for(unsigned j=argstart; j<getNumberOfArguments(); ++j) {
      double fforce = getPntrToOutput(0)->getForce(i);
      if( getPntrToArgument(j)->getRank()==0 ) {
        totv[nscalars] += fforce*getPntrToOutput(0)->getGridDerivative( i, getPntrToOutput(0)->getRank()+j ); nscalars++;
      } else { 
        double vval = getPntrToOutput(0)->getGridDerivative( i, getPntrToOutput(0)->getRank()+j  );
        getPntrToArgument(j)->addForce( i, fforce*vval );
      }
    }
  }
  nscalars=0; 
  for(unsigned i=argstart; i<getNumberOfArguments(); ++i) {
    if( getPntrToArgument(i)->getRank()==0 ) { getPntrToArgument(i)->addForce( 0, totv[nscalars] ); nscalars++; }
  }
}

}
}
#endif
