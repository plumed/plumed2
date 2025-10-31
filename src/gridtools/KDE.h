/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2017 The plumed team
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
#ifndef __PLUMED_gridtools_KDE_h
#define __PLUMED_gridtools_KDE_h

#include "ActionWithGrid.h"
#include "SumOfKernels.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "core/PbcAction.h"
#include "tools/HistogramBead.h"
#include "core/ParallelTaskManager.h"
#include "tools/SwitchingFunction.h"
#include "tools/Matrix.h"

namespace PLMD {
namespace gridtools {

template <class K, class P, class G>
class KDEHelper {
public:
  G g;
  bool fixed_width;
  std::size_t maxkernels;
  SumOfKernels<K, P> kernelsum;
  std::vector<unsigned> nneigh;
  std::vector<std::size_t> nkernels_per_point;
  std::vector<std::size_t> kernels_for_gridpoint;
  static void registerKeywords( Keywords& keys );
  static void read( KDEHelper<K,P,G>& func,
                    ActionWithArguments* action,
                    const std::vector<Value*>& args,
                    GridCoordinatesObject& gridobject,
                    std::vector<std::size_t>& shape,
                    function::FunctionOptions& options );
  static void readKernelParameters( std::string& value, ActionWithArguments* action, const std::string& outlab, bool rerequestargs );
  static void addArgument( const std::string& value, ActionWithArguments* action );
  static void setupGridBounds( KDEHelper<K,P,G>& func, const Tensor& box, GridCoordinatesObject& gridobject, const std::vector<Value*>& args, Value* myval );
  static void transferParamsToKernel( const std::vector<double>& argval, KDEHelper<K,P,G>& func, GridCoordinatesObject& gridobject, bool updateNeighborsOnEachKernel, std::size_t nkernels, unsigned kval, K& kp );
  static void transferKernels( KDEHelper<K,P,G>& func, const std::vector<Value*>& args, GridCoordinatesObject& gridobject );
};

template <class K, class P, class G>
void KDEHelper<K,P,G>::registerKeywords( Keywords& keys ) {
  SumOfKernels<K,P>::registerKeywords( keys );
  G::registerKeywords( keys );
}

template <class K, class P, class G>
void KDEHelper<K,P,G>::read( KDEHelper<K,P,G>& func,
                             ActionWithArguments* action,
                             const std::vector<Value*>& args,
                             GridCoordinatesObject& gridobject,
                             std::vector<std::size_t>& shape,
                             function::FunctionOptions& options ) {
  func.nneigh.resize( args.size() );
  SumOfKernels<K,P>::read( func.kernelsum, action, args, options );
  G::readBandwidthAndHeight( func.kernelsum.params, action );
  for(unsigned i=1; i<action->getNumberOfArguments(); ++i) {
    if( (action->getPntrToArgument(0))->getNumberOfValues()==0 && (action->getPntrToArgument(i))->isConstant() ) {
      continue;
    }

    if( (action->getPntrToArgument(0))->getNumberOfValues()!=1 || (action->getPntrToArgument(0))->getNumberOfValues()!=1 ) {
      if( (action->getPntrToArgument(0))->getRank()!=(action->getPntrToArgument(i))->getRank() ) {
        action->error("mismatch between ranks of input actions");
      }
      for(unsigned j=0; j<(action->getPntrToArgument(0))->getRank(); ++j) {
        if( (action->getPntrToArgument(0))->getShape()[j]!=(action->getPntrToArgument(i))->getShape()[j] ) {
          action->error("mismatch between shapes of input actions");
        }
      }
    }
  }
  G::readGridParameters( func.g, action, gridobject, shape );
}

template <class K, class P, class G>
void KDEHelper<K,P,G>::setupGridBounds( KDEHelper<K,P,G>& func, const Tensor& box, GridCoordinatesObject& gridobject, const std::vector<Value*>& args, Value* myval ) {
  // Setup the grid boundaries on first step
  G::setupGridBounds( func.g, box, gridobject, args, myval );
  // Check if the bandwidth changes during the simulation
  func.fixed_width = false;
  if( K::bandwidthIsConstant( gridobject.getDimension(), args ) && K::bandwidthsAllSame( gridobject.getDimension(), args ) ) {
    K myk;
    std::vector<double> myargs( args.size() );
    for(unsigned j=0; j<args.size(); ++j) {
      myargs[j] = args[j]->get(0);
    }
    K::setKernelAndCheckHeight( myk, gridobject.getDimension(), myargs );
    G::getDiscreteSupport( func.g, func.kernelsum.params, myk, func.nneigh, gridobject );
    func.fixed_width = true;
  }
  if( gridobject.getGridType()=="fibonacci" ) {
    return;
  }
  // Set the periodicity of the parameters
  for(unsigned i=0; i<gridobject.getDimension(); ++i) {
    P::setArgumentDomain( i, func.kernelsum.params, gridobject.getGridSpacing()[i], gridobject.isPeriodic(i), gridobject.getMin()[i], gridobject.getMax()[i] );
  }
}

template <class K, class P, class G>
void KDEHelper<K,P,G>::readKernelParameters( std::string& value, ActionWithArguments* action, const std::string& outlab, bool rerequestargs ) {
  std::size_t dot = value.find_first_of('.');
  ActionWithValue* av = action->plumed.getActionSet().selectWithLabel<ActionWithValue*>( value.substr(0,dot) );
  if( !av ) {
    std::string matstr, vals = "VALUES=" + value;
    if( (action->getPntrToArgument(0))->getRank()==2 ) {
      std::string nr, nc;
      Tools::convert( (action->getPntrToArgument(0))->getShape()[0], nr );
      Tools::convert( (action->getPntrToArgument(0))->getShape()[1], nc );
      matstr = " NROWS=" + nr + " NCOLS=" + nc;
    }
    for(unsigned i=1; i<(action->getPntrToArgument(0))->getNumberOfValues(); ++i) {
      vals += "," + value;
    }
    action->plumed.readInputWords( Tools::getWords(action->getLabel() + outlab + ": CONSTANT " + vals + matstr ), false );
    value = action->getLabel() + outlab;
  } else {
    Value* myval;
    if( dot!=std::string::npos ) {
      myval = av->copyOutput( value );
    } else {
      if( av->getNumberOfComponents()>1 ) {
        action->error("problem reading argument " + value );
      }
      myval = av->copyOutput(0);
    }
    if( myval->getRank()==0 ) {
      std::string nvals;
      if( (action->getPntrToArgument(0))->getRank()==2 ) {
        std::string nr, nc;
        Tools::convert( (action->getPntrToArgument(0))->getShape()[0], nr );
        Tools::convert( (action->getPntrToArgument(0))->getShape()[1], nc );
        nvals = nr + "," + nc;
      } else {
        Tools::convert( (action->getPntrToArgument(0))->getNumberOfValues(), nvals );
      }
      action->plumed.readInputWords( Tools::getWords(action->getLabel() + outlab + "_ones: ONES SIZE=" + nvals ), false );
      action->plumed.readInputWords( Tools::getWords(action->getLabel() + outlab + ": CUSTOM ARG=" + action->getLabel() + outlab + "_ones," + value ), false );
      value = action->getLabel() + outlab;
    }
  }
  if( !rerequestargs ) {
    return;
  }
  KDEHelper<K,P,G>::addArgument( value, action );
}

template <class K, class P, class G>
void KDEHelper<K,P,G>::addArgument( const std::string& value, ActionWithArguments* action ) {
  std::size_t dot = value.find_first_of(".");
  ActionWithValue* av = action->plumed.getActionSet().selectWithLabel<ActionWithValue*>( value.substr(0,dot) );
  plumed_assert( av );
  std::vector<Value*> args( action->getArguments() );
  args.push_back( av->copyOutput(0) );
  action->requestArguments( args );
}

template <class K, class P, class G>
void KDEHelper<K,P,G>::transferParamsToKernel( const std::vector<double>& argval, KDEHelper<K,P,G>& func, GridCoordinatesObject& gridobject, bool updateNeighborsOnEachKernel, std::size_t nkernels, unsigned kval, K& kp ) {
  // This sets the kernel parameters for the Kth kernel and checks that we want
  // to consider it
  if( !K::setKernelAndCheckHeight( kp, gridobject.getDimension(), argval ) ) {
    return;
  }
  // If the widths of each kernel are not all the same then get the discrete support
  if( updateNeighborsOnEachKernel ) {
    G::getDiscreteSupport( func.g, func.kernelsum.params, kp, func.nneigh, gridobject );
  }

  // Now get the grid points for this particular kernel
  unsigned num_neigh;
  std::vector<unsigned> neighbors;
  G::getNeighbors( func.kernelsum.params, kp, gridobject, func.nneigh, num_neigh, neighbors );

  // And transfer the neighbor information to the holders
  for(unsigned j=0; j<num_neigh; ++j) {
    func.kernels_for_gridpoint[ neighbors[j]*nkernels + func.nkernels_per_point[neighbors[j]] ] = kval;
    func.nkernels_per_point[ neighbors[j] ]++;
  }
}

template <class K, class P, class G>
void KDEHelper<K,P,G>::transferKernels( KDEHelper<K,P,G>& func, const std::vector<Value*>& args, GridCoordinatesObject& gridobject ) {
  // Resize the kernel sum if we need to
  // Number of kernels is determined based on sparsity pattern of matrix input as matrix of heights
  std::size_t nkernels = args[args.size()-1]->getNumberOfStoredValues();
  if( func.kernelsum.kernelParams.size()!=nkernels ) {
    func.kernelsum.kernelParams.resize( nkernels );
  }
  // And resize the grid counters if we need to
  std::size_t ngp = gridobject.getNumberOfPoints();
  if( func.nkernels_per_point.size()!=ngp ) {
    func.nkernels_per_point.resize( ngp );
    func.kernels_for_gridpoint.resize( ngp*nkernels );
  }
  std::fill( func.nkernels_per_point.begin(), func.nkernels_per_point.end(), 0 );

  bool updateNeighborsOnEachKernel = !func.fixed_width;
  if( !func.fixed_width && K::bandwidthsAllSame( gridobject.getDimension(), args ) ) {
    G::getDiscreteSupport( func.g, func.kernelsum.params, func.kernelsum.kernelParams[0], func.nneigh, gridobject );
    updateNeighborsOnEachKernel = false;
  }

  std::vector<double> argval( args.size() );
  if( args[args.size()-1]->getRank()==2 ) {
    const unsigned nc    = args[args.size()-1]->getShape()[1];
    const unsigned nrows = args[args.size()-1]->getShape()[0];
    const unsigned ncs   = args[args.size()-1]->getNumberOfColumns();
    for(unsigned i=0; i<nrows; ++i) {
      unsigned ncols = args[args.size()-1]->getRowLength(i);
      for(unsigned k=0; k<args.size(); ++k) {
        plumed_massert( args[k]->isConstant() || ncols==args[k]->getRowLength(i), "all input matrices must have same sparsity pattern" );
      }
      for(unsigned j=0; j<ncols; ++j) {
        unsigned jind = args[args.size()-1]->getRowIndex( i, j );
        for(unsigned k=0; k<args.size(); ++k) {
          if( jind==args[k]->getRowIndex( i, j ) ) {
            argval[k] = args[k]->get( i*ncs+ j, false );
          } else {
            argval[k] = args[k]->get( i*nc + jind );
          }
        }
        KDEHelper<K,P,G>::transferParamsToKernel( argval, func, gridobject, updateNeighborsOnEachKernel, nkernels, i*ncs+j, func.kernelsum.kernelParams[i*ncs+j] );
      }
    }
  } else {
    for(unsigned i=0; i<nkernels; ++i) {
      // Transfer the kernel parameters to local vector of doubles
      for(unsigned j=0; j<args.size(); ++j) {
        argval[j] = args[j]->get(i,false);
      }
      KDEHelper<K,P,G>::transferParamsToKernel( argval, func, gridobject, updateNeighborsOnEachKernel, nkernels, i, func.kernelsum.kernelParams[i] );
    }
  }
  // Get the maximum number of kernels for any given grid point (used for resizing derivatives)
  func.maxkernels = 0;
  for(unsigned i=0; i<ngp; ++i) {
    if( func.nkernels_per_point[i]>func.maxkernels ) {
      func.maxkernels = func.nkernels_per_point[i];
    }
  }
}

template <class K, class P, class G>
class KDE : public ActionWithGrid {
public:
  using input_type = KDEHelper<K, P, G>;
  using PTM = ParallelTaskManager<KDE<K,P,G>>;
private:
  bool firststep;
/// The parallel task manager
  PTM taskmanager;
  GridCoordinatesObject gridobject;
public:
  static void registerKeywords( Keywords& keys );
  explicit KDE(const ActionOptions&ao);
  std::vector<std::string> getGridCoordinateNames() const override ;
  const GridCoordinatesObject& getGridCoordinatesObject() const override ;
  unsigned getNumberOfDerivatives() override;
  int checkTaskIsActive( const unsigned& itask ) const override ;
  void prepare() override ;
  void calculate() override ;
  void getInputData( std::vector<double>& inputdata ) const override ;
  static void performTask( std::size_t task_index,
                           const KDEHelper<K, P, G>& actiondata,
                           ParallelActionsInput& input,
                           ParallelActionsOutput& output );
  void applyNonZeroRankForces( std::vector<double>& outforces ) override ;
  static int getNumberOfValuesPerTask( std::size_t task_index,
                                       const KDEHelper<K, P, G>& actiondata );
  static void getForceIndices( std::size_t task_index,
                               std::size_t colno,
                               std::size_t ntotal_force,
                               const KDEHelper<K, P, G>& actiondata,
                               const ParallelActionsInput& input,
                               ForceIndexHolder force_indices );
};

template <class K, class P, class G>
void KDE<K,P,G>::registerKeywords( Keywords& keys ) {
  ActionWithGrid::registerKeywords( keys );
  keys.addInputKeyword("compulsory","ARG","scalar/vector/matrix","the label for the value that should be used to construct the histogram");
  KDEHelper<K,P,G>::registerKeywords( keys );
  // Keywords for spherical KDE
  keys.add("hidden","MASKED_INPUT_ALLOWED","turns on that you are allowed to use masked inputs ");
  keys.setValueDescription("grid","a function on a grid that was obtained by doing a Kernel Density Estimation using the input arguments");
  if( keys.getDisplayName()!="SPHERICAL_KDE" ) {
    keys.setDisplayName("KDE");
  }
  PTM::registerKeywords( keys );
}

template <class K, class P, class G>
KDE<K,P,G>::KDE(const ActionOptions&ao):
  Action(ao),
  ActionWithGrid(ao),
  firststep(true),
  taskmanager(this) {

  std::vector<std::size_t> shape( getNumberOfArguments() );
  unsigned numberOfKernels=getPntrToArgument(0)->getNumberOfValues();
  for(unsigned i=1; i<shape.size(); ++i) {
    if( numberOfKernels!=getPntrToArgument(i)->getNumberOfValues() ) {
      error("mismatch between numbers of values in input arguments");
    }
  }

  function::FunctionOptions foptions;
  KDEHelper<K,P,G>::read( taskmanager.getActionInput(), this, getArguments(), gridobject, shape, foptions );
  addValueWithDerivatives( shape );
  setNotPeriodic();
  getPntrToComponent(0)->setDerivativeIsZeroWhenValueIsZero();
}

template <class K, class P, class G>
unsigned KDE<K,P,G>::getNumberOfDerivatives() {
  return gridobject.getDimension();
}

template <class K, class P, class G>
std::vector<std::string> KDE<K,P,G>::getGridCoordinateNames() const {
  std::vector<std::string> names( gridobject.getDimension() );
  for(unsigned i=0; i<names.size(); ++i) {
    names[i] = getPntrToArgument(i)->getName();
  }
  return names;
}

template <class K, class P, class G>
const GridCoordinatesObject& KDE<K,P,G>::getGridCoordinatesObject() const {
  return gridobject;
}

template <class K, class P, class G>
int KDE<K,P,G>::checkTaskIsActive( const unsigned& itask ) const {
  if( taskmanager.getActionInput().nkernels_per_point[itask]>0 ) {
    return 1;
  }
  return -1;
}

template <class K, class P, class G>
void KDE<K,P,G>::prepare() {
  ActionWithVector::prepare();
  std::size_t nkernels = getPntrToArgument(0)->getNumberOfValues();
  for(unsigned i=1; i<getNumberOfArguments(); ++i) {
    Value* myarg = getPntrToArgument(i);
    if( myarg->getNumberOfValues()!=nkernels ) {
      if( myarg->isConstant() && myarg->getNumberOfValues()==1 ) {
        myarg->reshapeConstantValue( getPntrToArgument(0)->getShape() );
      } else {
        plumed_merror("found mismatched numbers of arguments in input");
      }
    }
  }
}

template <class K, class P, class G>
void KDE<K,P,G>::calculate() {
  if( firststep ) {
    PbcAction* bv = plumed.getActionSet().template selectWithLabel<PbcAction*>("Box");
    KDEHelper<K,P,G>::setupGridBounds( taskmanager.getActionInput(), bv->getPbc().getBox(), gridobject, getArguments(), getPntrToComponent(0) );
    firststep=false;
  }
  KDEHelper<K,P,G>::transferKernels( taskmanager.getActionInput(), getArguments(), gridobject );
  taskmanager.setupParallelTaskManager( getNumberOfArguments()*taskmanager.getActionInput().maxkernels, getNumberOfForceDerivatives() );
  taskmanager.runAllTasks();
}

template <class K, class P, class G>
void KDE<K,P,G>::getInputData( std::vector<double>& inputdata ) const {
  std::size_t ndim = gridobject.getDimension();
  std::size_t nstored = getConstPntrToComponent(0)->getNumberOfStoredValues();
  std::vector<double> pos( ndim );
  if( inputdata.size()!=nstored*ndim ) {
    inputdata.resize( ndim*nstored );
  }

  for(unsigned i=0; i<nstored; ++i) {
    gridobject.getGridPointCoordinates( i, pos );
    for(unsigned j=0; j<ndim; ++j) {
      inputdata[ i*ndim + j ] = pos[j];
    }
  }
}

template <class K, class P, class G>
void KDE<K,P,G>::performTask( std::size_t task_index,
                              const KDEHelper<K, P, G>& actiondata,
                              ParallelActionsInput& input,
                              ParallelActionsOutput& output ) {
  std::size_t ndim = actiondata.nneigh.size();
  SumOfKernels<K,P>::calc( View<const std::size_t>( actiondata.kernels_for_gridpoint.data() + task_index*input.argstarts[1], actiondata.nkernels_per_point[task_index] ),
                           actiondata.kernelsum,
                           View<const double>(input.inputdata + task_index*actiondata.nneigh.size(),ndim),
                           View<double>(output.values.data(), 1),
                           View<double>(output.values.data()+1, ndim),
                           View<double>(output.derivatives.data(),actiondata.maxkernels*input.nargs) );

}

template <class K, class P, class G>
void KDE<K,P,G>::applyNonZeroRankForces( std::vector<double>& outforces ) {
  taskmanager.applyForces( outforces );
}

template <class K, class P, class G>
int KDE<K,P,G>::getNumberOfValuesPerTask( std::size_t task_index,
    const KDEHelper<K, P, G>& actiondata ) {
  return 1;
}

template <class K, class P, class G>
void KDE<K,P,G>::getForceIndices( std::size_t task_index,
                                  std::size_t colno,
                                  std::size_t ntotal_force,
                                  const KDEHelper<K, P, G>& actiondata,
                                  const ParallelActionsInput& input,
                                  ForceIndexHolder force_indices ) {
  force_indices.threadsafe_derivatives_end[0] = 0;
  std::size_t nparams = K::getNumberOfParameters( actiondata.kernelsum.kernelParams[0] );
  View<const std::size_t> kernellist( actiondata.kernels_for_gridpoint.data() + task_index*input.argstarts[1], actiondata.nkernels_per_point[task_index] );
  for(unsigned i=0; i<kernellist.size(); ++i) {
    for(unsigned j=0; j<nparams; ++j) {
      force_indices.indices[0][i*nparams+j] = input.argstarts[j] + kernellist[i];
    }
  }
  force_indices.tot_indices[0] = kernellist.size()*nparams;
}

}
}
#endif
