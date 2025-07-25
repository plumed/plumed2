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
#include "core/ParallelTaskManager.h"
#include "function/FunctionSetup.h"
#include "function/Custom.h"
#include "EvaluateGridFunction.h"

namespace PLMD {
namespace gridtools {

template <class T>
class FunctionOfGrid : public ActionWithGrid {
public:
  using input_type = function::FunctionData<T>;
  using PTM = ParallelTaskManager<FunctionOfGrid<T>>;
private:
/// The parallel task manager
  PTM taskmanager;
//// Ensures we setup on first step
  bool firststep;
/// Set equal to one if we are doing EvaluateGridFunction
  unsigned argstart;
public:
  static void registerKeywords(Keywords&);
  explicit FunctionOfGrid(const ActionOptions&);
/// This does setup required on first step
  void setupOnFirstStep( const bool incalc );
/// Get the number of derivatives for this action
  unsigned getNumberOfDerivatives() override ;
/// Get the label to write in the graph
  std::string writeInGraph() const override ;
/// Get the underlying names
  std::vector<std::string> getGridCoordinateNames() const override ;
/// Get the underlying grid coordinates object
  const GridCoordinatesObject& getGridCoordinatesObject() const override ;
/// Calculate the function
  void performTask( const unsigned& current, MultiValue& myvals ) const override { plumed_error(); }
/// Get the input data for doing the parallel calculation
  void getInputData( std::vector<double>& inputdata ) const override ;
/// Do the calculation
  void calculate() override ;
// Calculate the value of the function at a grid point
  static void performTask( std::size_t task_index,
                           const function::FunctionData<T>& actiondata,
                           ParallelActionsInput& input,
                           ParallelActionsOutput& output );
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
  PTM::registerKeywords( keys );
}

template <class T>
FunctionOfGrid<T>::FunctionOfGrid(const ActionOptions&ao):
  Action(ao),
  ActionWithGrid(ao),
  taskmanager(this),
  firststep(true),
  argstart(0) {
  if( getNumberOfArguments()==0 ) {
    error("found no arguments");
  }
  if( typeid(taskmanager.getActionInput().f)==typeid(EvaluateGridFunction) ) {
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
  auto & myfunc = taskmanager.getActionInput();
  myfunc.argstart = argstart;
  myfunc.nscalars = 0;
  // Create the values for this grid
  function::FunctionData<T>::setup( myfunc.f, keywords.getOutputComponents(), shape, true, this  );
  // Setup the task manager
  taskmanager.setupParallelTaskManager( getNumberOfArguments()-argstart, 0 );
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
void FunctionOfGrid<T>::calculate() {
  if( firststep ) {
      setupOnFirstStep( true );
      firststep=false;
  } 
  taskmanager.runAllTasks();
}

template <class T>
void FunctionOfGrid<T>::getInputData( std::vector<double>& inputdata ) const {
  unsigned nargs = getNumberOfArguments();

  std::size_t ntasks = 0, ngder = 0;
  for(unsigned i=argstart; i<nargs; ++i) {
    if( getPntrToArgument(i)->getRank()>0 ) {
      ntasks = getPntrToArgument(i)->getNumberOfStoredValues();
      ngder = getPntrToArgument(i)->getNumberOfDerivatives();
      break;
    }
  }

  std::size_t ndata = static_cast<std::size_t>(nargs-argstart)*ntasks*(1+ngder); 
  if( inputdata.size()!=ndata ) {
    inputdata.resize( ndata );
  }

  for(unsigned j=argstart; j<nargs; ++j) {
    const Value* myarg =  getPntrToArgument(j);
    if( myarg->getRank()==0 ) {
      double val = myarg->get();
      for(unsigned i=0; i<ntasks; ++i) {
        inputdata[(nargs-argstart)*(1+ngder)*i + j-argstart] = val;
      }
    } else {
      for(unsigned i=0; i<ntasks; ++i) {
        inputdata[(nargs-argstart)*(1+ngder)*i + j-argstart] = myarg->get(i);
        for(unsigned k=0; k<ngder; ++k) {
            inputdata[(nargs-argstart)*(1+ngder)*i + (nargs-argstart) + (j-argstart)*ngder + k] = myarg->getGridDerivative( i, k );
        }
      }
    }
  }
}

template <class T>
void FunctionOfGrid<T>::performTask( std::size_t task_index,
                                     const function::FunctionData<T>& actiondata,
                                     ParallelActionsInput& input,
                                     ParallelActionsOutput& output ) {

  std::size_t rank=0; 
  for(unsigned j=actiondata.argstart; j<input.nargs; ++j) {
      if( input.ranks[j]>0 ) {
          rank=input.ranks[j];
          break;
      }
  } 
  std::size_t spacing = (input.nargs-actiondata.argstart)*(1+rank);
  auto funcout = function::FunctionOutput::create( input.ncomponents,
                                                   output.values.data(),
                                                   input.nargs-actiondata.argstart,
                                                   output.values.data()+1+rank );
  T::calc( actiondata.f,
           input.noderiv,
           View<const double>( input.inputdata + task_index*spacing,
                               input.nargs-actiondata.argstart ),
           funcout );

  for(unsigned j=actiondata.argstart; j<input.nargs; ++j) {
      double df = output.values[1+rank+j];
      View<const double> inders( input.inputdata + task_index*spacing + (input.nargs-actiondata.argstart) + (j-actiondata.argstart)*rank, rank );
      for(unsigned k=0; k<input.ranks[actiondata.argstart]; ++k) {
          output.values[1+k] += df*inders[k];
      }
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
