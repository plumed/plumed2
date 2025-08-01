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

#include "ActionWithGrid.h"
#include "function/FunctionSetup.h"
#include "function/Custom.h"
#include "EvaluateGridFunction.h"

namespace PLMD {
namespace gridtools {

template <class T>
class FunctionOfGrid : public ActionWithGrid {
private:
/// Set equal to one if we are doing EvaluateGridFunction
  unsigned argstart;
/// The function that is being computed
  T myfunc;
public:
  static void registerKeywords(Keywords&);
  explicit FunctionOfGrid(const ActionOptions&);
/// This does setup required on first step
  void setupOnFirstStep( const bool incalc ) override ;
/// Get the number of derivatives for this action
  unsigned getNumberOfDerivatives() override ;
/// Get the label to write in the graph
  std::string writeInGraph() const override ;
/// Get the underlying names
  std::vector<std::string> getGridCoordinateNames() const override ;
/// Get the underlying grid coordinates object
  const GridCoordinatesObject& getGridCoordinatesObject() const override ;
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
  ActionWithGrid::registerKeywords(keys);
  std::string name = keys.getDisplayName();
  std::size_t und=name.find("_GRID");
  keys.setDisplayName( name.substr(0,und) );
  keys.reserve("compulsory","PERIODIC","if the output of your function is periodic then you should specify the periodicity of the function.  If the output is not periodic you must state this using PERIODIC=NO");
  T tfunc;
  T::registerKeywords( keys );
  if( typeid(tfunc)==typeid(function::Custom()) ) {
    keys.add("hidden","NO_ACTION_LOG","suppresses printing from action on the log");
  }
  if( keys.outputComponentExists(".#!value") ) {
    keys.setValueDescription("grid","the grid obtained by doing an element-wise application of " + keys.getOutputComponentDescription(".#!value") + " to the input grid");
    keys.addInputKeyword("compulsory","ARG","scalar/grid","the labels of the scalars and functions on a grid that we are using to compute the required function");
  }
}

template <class T>
FunctionOfGrid<T>::FunctionOfGrid(const ActionOptions&ao):
  Action(ao),
  ActionWithGrid(ao),
  argstart(0) {
  if( getNumberOfArguments()==0 ) {
    error("found no arguments");
  }
  if( typeid(myfunc)==typeid(EvaluateGridFunction) ) {
    argstart=1;
  }
  // This will require a fix
  if( getPntrToArgument(argstart)->getRank()==0 || !getPntrToArgument(argstart)->hasDerivatives() ) {
    error("first input to this action must be a grid");
  }
  // Get the shape of the input grid
  std::vector<std::size_t> shape( getPntrToArgument(argstart)->getShape() );
  for(unsigned i=argstart+1; i<getNumberOfArguments(); ++i ) {
    if( getPntrToArgument(i)->getRank()==0 ) {
      continue;
    }
    std::vector<std::size_t> s( getPntrToArgument(i)->getShape() );
    if( s.size()!=shape.size() ) {
      error("mismatch between dimensionalities of input grids");
    }
  }
  // Create the values for this grid
  function::FunctionData<T>::setup( myfunc, keywords.getOutputComponents(), shape, true, this  );
  // And setup on first step
  setupOnFirstStep( false );
}

template <class T>
std::string FunctionOfGrid<T>::writeInGraph() const {
  std::size_t und = getName().find_last_of("_");
  return getName().substr(0,und);
}

template <class T>
void FunctionOfGrid<T>::setupOnFirstStep( const bool incalc ) {
  const GridCoordinatesObject& mygrid = getGridCoordinatesObject();
  unsigned npoints = getPntrToArgument(0)->getNumberOfValues();
  if( mygrid.getGridType()=="flat" ) {
    std::vector<std::size_t> shape( getGridCoordinatesObject().getNbin(true) );
    for(unsigned i=1; i<getNumberOfArguments(); ++i ) {
      if( getPntrToArgument(i)->getRank()==0 ) {
        continue;
      }
      std::vector<std::size_t> s( getPntrToArgument(i)->getShape() );
      for(unsigned j=0; j<shape.size(); ++j) {
        if( shape[j]!=s[j] ) {
          error("mismatch between sizes of input grids");
        }
      }
    }
    for(int i=0; i<getNumberOfComponents(); ++i) {
      if( getPntrToComponent(i)->getRank()>0 ) {
        getPntrToComponent(i)->setShape(shape);
      }
    }
  }
  // This resizes the scalars
  for(int i=0; i<getNumberOfComponents(); ++i) {
    if( getPntrToComponent(i)->getRank()==0 ) {
      getPntrToComponent(i)->resizeDerivatives( npoints );
    }
  }
}

template <class T>
const GridCoordinatesObject& FunctionOfGrid<T>::getGridCoordinatesObject() const {
  ActionWithGrid* ag=ActionWithGrid::getInputActionWithGrid( getPntrToArgument(0)->getPntrToAction() );
  plumed_assert( ag );
  return ag->getGridCoordinatesObject();
}

template <class T>
std::vector<std::string> FunctionOfGrid<T>::getGridCoordinateNames() const {
  ActionWithGrid* ag=ActionWithGrid::getInputActionWithGrid( getPntrToArgument(0)->getPntrToAction() );
  plumed_assert( ag );
  return ag->getGridCoordinateNames();
}

template <class T>
unsigned FunctionOfGrid<T>::getNumberOfDerivatives() {
  unsigned nder = getGridCoordinatesObject().getDimension();
  return getGridCoordinatesObject().getDimension() + getNumberOfArguments() - argstart;
}

template <class T>
void FunctionOfGrid<T>::performTask( const unsigned& current, MultiValue& myvals ) const {
  std::vector<double> args( getNumberOfArguments() - argstart );
  for(unsigned i=argstart; i<getNumberOfArguments(); ++i) {
    if( getPntrToArgument(i)->getRank()==0 ) {
      args[i-argstart]=getPntrToArgument(i)->get();
    } else {
      args[i-argstart] = getPntrToArgument(i)->get(current);
    }
  }
  // Calculate the function and its derivatives
  std::vector<double> vals(getNumberOfComponents()), deriv( getNumberOfComponents()*args.size() );
  auto funcout = function::FunctionOutput::create( getNumberOfComponents(),
                 vals.data(),
                 args.size(),
                 deriv.data() );
  T::calc( myfunc,
           false,
           View<const double>(args.data(), args.size()),
           funcout );
  unsigned np = myvals.getTaskIndex();
  // And set the values and derivatives
  myvals.addValue( 0, vals[0] );
  // Add the derivatives for a grid
  for(unsigned j=argstart; j<getNumberOfArguments(); ++j) {
    // We store all the derivatives of all the input values - i.e. the grid points these are used in apply
    myvals.addDerivative( 0, getConstPntrToComponent(0)->getRank()+j-argstart, funcout.derivs[0][j-argstart] );
    // And now we calculate the derivatives of the value that is stored on the grid correctly so that we can interpolate functions
    for(unsigned k=0; k<getPntrToArgument(j)->getRank(); ++k) {
      myvals.addDerivative( 0, k, funcout.derivs[0][j-argstart]*getPntrToArgument(j)->getGridDerivative( np, k ) );
    }
  }
  unsigned nderivatives = getConstPntrToComponent(0)->getNumberOfGridDerivatives();
  for(unsigned j=0; j<nderivatives; ++j) {
    myvals.updateIndex( 0, j );
  }
}

template <class T>
void FunctionOfGrid<T>::gatherStoredValue( const unsigned& valindex, const unsigned& code, const MultiValue& myvals,
    const unsigned& bufstart, std::vector<double>& buffer ) const {
  if( getConstPntrToComponent(0)->getRank()>0 && getConstPntrToComponent(0)->hasDerivatives() ) {
    plumed_dbg_assert( getNumberOfComponents()==1 && valindex==0 );
    unsigned nder = getConstPntrToComponent(0)->getNumberOfGridDerivatives();
    unsigned kp = bufstart + code*(1+nder);
    buffer[kp] += myvals.get( 0 );
    for(unsigned i=0; i<nder; ++i) {
      buffer[kp + 1 + i] += myvals.getDerivative( 0, i );
    }
  } else {
    ActionWithVector::gatherStoredValue( valindex, code, myvals, bufstart, buffer );
  }
}

template <class T>
void FunctionOfGrid<T>::apply() {
  if( doNotCalculateDerivatives() || !getPntrToComponent(0)->forcesWereAdded() ) {
    return;
  }

  // Work out how to deal with arguments
  unsigned nscalars=0;
  for(unsigned i=argstart; i<getNumberOfArguments(); ++i) {
    if( getPntrToArgument(i)->getRank()==0 ) {
      nscalars++;
    }
  }

  std::vector<double> totv(nscalars,0);
  Value* outval=getPntrToComponent(0);
  for(unsigned i=0; i<outval->getNumberOfValues(); ++i) {
    nscalars=0;
    for(unsigned j=argstart; j<getNumberOfArguments(); ++j) {
      double fforce = outval->getForce(i);
      if( getPntrToArgument(j)->getRank()==0 ) {
        totv[nscalars] += fforce*outval->getGridDerivative( i, outval->getRank()+j );
        nscalars++;
      } else {
        double vval = outval->getGridDerivative( i, outval->getRank()+j  );
        getPntrToArgument(j)->addForce( i, fforce*vval );
      }
    }
  }
  nscalars=0;
  for(unsigned i=argstart; i<getNumberOfArguments(); ++i) {
    if( getPntrToArgument(i)->getRank()==0 ) {
      getPntrToArgument(i)->addForce( 0, totv[nscalars] );
      nscalars++;
    }
  }

}

}
}
#endif
