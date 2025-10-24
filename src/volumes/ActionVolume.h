/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2020 The plumed team
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
#ifndef __PLUMED_volumes_ActionVolume_h
#define __PLUMED_volumes_ActionVolume_h

#include "core/ActionWithVector.h"
#include "core/ParallelTaskManager.h"
#include "tools/ColvarOutput.h"

namespace PLMD {
namespace volumes {

template <class T>
struct VolumeData {
  bool not_in;
  std::size_t numberOfNonReferenceAtoms;
  T voldata;
#ifdef __PLUMED_USE_OPENACC
  void toACCDevice() const {
#pragma acc enter data copyin(this[0:1],not_in,numberOfNonReferenceAtoms)
    voldata.toACCDevice();
  }
  void removeFromACCDevice() const {
    voldata.removeFromACCDevice();
#pragma acc exit data delete(numberOfNonReferenceAtoms,not_in,this[0:1])
  }
#endif //__PLUMED_USE_OPENACC
};

struct VolumeInput {
  std::size_t task_index;
  const Pbc& pbc;
  View<const double,3> cpos;
  View2D<const double,helpers::dynamic_extent,3> refpos;
  VolumeInput( std::size_t t, unsigned nref, double* p, double* rp, const Pbc& box ) :
    task_index(t),pbc(box),cpos(p),refpos(rp,nref) {
  }
};

class VolumeOutput {
private:
  class RefderHelper {
  private:
    double* derivatives;
  public:
    RefderHelper( double* d ) : derivatives(d) {}
    View<double,3> operator[](std::size_t i) {
      return View<double,3>( derivatives + 3*(i+1) );
    }
  };
  ColvarOutput fulldata;
public:
  View<double>& values;
  ColvarOutput::VirialHelper& virial;
  View<double,3> derivatives;
  RefderHelper refders;
  VolumeOutput( View<double>& v, std::size_t nder, double* d ) :
    fulldata(v, nder, d),
    values(fulldata.values),
    virial(fulldata.virial),
    derivatives(d), refders(d) {
  }
};

/**
\ingroup INHERIT
This is the abstract base class to use for implementing a new way of defining a particular region of the simulation
box. You can use this to calculate the number of atoms inside that part or the average value of a quantity like the
coordination number inside that part of the cell.
*/

template <class T>
class ActionVolume : public ActionWithVector {
public:
  using input_type = VolumeData<T>;
  using PTM = ParallelTaskManager<ActionVolume<T>>;
private:
/// The parallel task manager
  PTM taskmanager;
public:
  static void registerKeywords( Keywords& keys );
  explicit ActionVolume(const ActionOptions&);
  unsigned getNumberOfDerivatives() override;
  void getInputData( std::vector<double>& inputdata ) const override ;
  void calculate() override ;
  void applyNonZeroRankForces( std::vector<double>& outforces ) override ;
  static void performTask( std::size_t task_index, const VolumeData<T>& actiondata, ParallelActionsInput& input, ParallelActionsOutput& output );
  static void gatherForces( std::size_t task_index, const VolumeData<T>& actiondata, const ParallelActionsInput& input, const ForceInput& fdata, ForceOutput& forces );
  static int getNumberOfValuesPerTask( std::size_t task_index, const VolumeData<T>& actiondata );
  static void getForceIndices( std::size_t task_index, std::size_t colno, std::size_t ntotal_force, const VolumeData<T>& actiondata, const ParallelActionsInput& input, ForceIndexHolder force_indices );
};

template <class T>
unsigned ActionVolume<T>::getNumberOfDerivatives() {
  return 3*getNumberOfAtoms()+9;
}

template <class T>
void ActionVolume<T>::registerKeywords( Keywords& keys ) {
  ActionWithVector::registerKeywords( keys );
  PTM::registerKeywords( keys );
  keys.add("atoms","ATOMS","the group of atoms that you would like to investigate");
  keys.addFlag("OUTSIDE",false,"calculate quantities for colvars that are on atoms outside the region of interest");
  keys.setValueDescription("scalar/vector","vector of numbers between 0 and 1 that measure the degree to which each atom is within the volume of interest");
  T::registerKeywords( keys );
}

template <class T>
ActionVolume<T>::ActionVolume(const ActionOptions&ao):
  Action(ao),
  ActionWithVector(ao),
  taskmanager(this) {
  std::vector<AtomNumber> atoms;
  parseAtomList("ATOMS",atoms);
  if( atoms.size()==0 ) {
    error("no atoms were specified");
  }
  log.printf("  examining positions of atoms ");
  for(unsigned i=0; i<atoms.size(); ++i) {
    log.printf(" %d", atoms[i].serial() );
  }
  log.printf("\n");
  std::vector<std::size_t> shape(1);
  shape[0]=atoms.size();

  std::vector<AtomNumber> refatoms;
  T::parseAtoms( this, refatoms );
  for(unsigned i=0; i<refatoms.size(); ++i) {
    atoms.push_back( refatoms[i] );
  }
  requestAtoms( atoms );
  VolumeData<T> actioninput;
  actioninput.voldata.parseInput(this);

  actioninput.numberOfNonReferenceAtoms=shape[0];
  parseFlag("OUTSIDE",actioninput.not_in);

  if( shape[0]==1 ) {
    ActionWithValue::addValueWithDerivatives();
  } else {
    ActionWithValue::addValue( shape );
    taskmanager.setupParallelTaskManager( 3*(1+refatoms.size())+9, 3*refatoms.size()+9 );
  }
  setNotPeriodic();
  getPntrToComponent(0)->setDerivativeIsZeroWhenValueIsZero();

  taskmanager.setActionInput( actioninput );
}

template <class T>
void ActionVolume<T>::getInputData( std::vector<double>& inputdata ) const {
  if( inputdata.size()!=3*getNumberOfAtoms() ) {
    inputdata.resize( 3*getNumberOfAtoms() );
  }

  for(unsigned i=0; i<getNumberOfAtoms(); ++i) {
    Vector ipos = getPosition(i);
    for(unsigned j=0; j<3; ++j) {
      inputdata[3*i+j] = ipos[j];
    }
  }
}

template <class T>
void ActionVolume<T>::calculate() {
  unsigned k=0;
  std::vector<Vector> positions( getNumberOfAtoms()-getPntrToComponent(0)->getNumberOfValues() );
  for(unsigned i=getPntrToComponent(0)->getNumberOfValues(); i<getNumberOfAtoms(); ++i) {
    positions[k] = getPosition( i );
    k++;
  }
  taskmanager.getActionInput().voldata.setupRegions( this, getPbc(), positions );

  if( getPntrToComponent(0)->getRank()==0 ) {
    std::size_t nref = getNumberOfAtoms() - 1;
    std::vector<double> posvec;
    getInputData( posvec );
    std::vector<double> deriv( getNumberOfDerivatives() );
    std::vector<double> val(1);
    View<double> valview(val.data(),1);
    VolumeOutput output( valview, getNumberOfDerivatives(), deriv.data() );
    T::calculateNumberInside( VolumeInput( 0, nref, posvec.data(), posvec.data()+3, getPbc() ), taskmanager.getActionInput().voldata, output );
    if( taskmanager.getActionInput().not_in ) {
      val[0] = 1.0 - val[0];
      for(unsigned i=0; i<deriv.size(); ++i) {
        deriv[i] *= -1;
      }
    }
    Value* v = getPntrToComponent(0);
    v->set( val[0] );
    for(unsigned i=0; i<deriv.size(); ++i) {
      v->addDerivative( i, deriv[i] );
    }
  } else {
    taskmanager.runAllTasks();
  }
}

template <class T>
void ActionVolume<T>::applyNonZeroRankForces( std::vector<double>& outforces ) {
  taskmanager.applyForces( outforces );
}

template <class T>
void ActionVolume<T>::performTask( std::size_t task_index, const VolumeData<T>& actiondata, ParallelActionsInput& input, ParallelActionsOutput& output ) {
  std::size_t nref = output.derivatives.size()/3 - 4; // This is the number of reference atoms
  VolumeOutput volout( output.values, output.derivatives.size(), output.derivatives.data() );
  T::calculateNumberInside( VolumeInput( task_index, nref, input.inputdata+3*task_index, input.inputdata+3*actiondata.numberOfNonReferenceAtoms, *input.pbc ), actiondata.voldata, volout );

  if( actiondata.not_in ) {
    output.values[0] = 1.0 - output.values[0];
    if( input.noderiv ) {
      return;
    }
    for(unsigned i=0; i<output.derivatives.size(); ++i) {
      output.derivatives[i] *= -1;
    }
  }
}

template <class T>
int ActionVolume<T>::getNumberOfValuesPerTask( std::size_t task_index,
    const VolumeData<T>& actiondata ) {
  return 1;
}

template<class T>
void ActionVolume<T>::getForceIndices( std::size_t task_index,
                                       std::size_t colno,
                                       std::size_t ntotal_force,
                                       const VolumeData<T>& actiondata,
                                       const ParallelActionsInput& input,
                                       ForceIndexHolder force_indices ) {
  std::size_t base = 3*task_index;
  force_indices.indices[0][0] = base;
  force_indices.indices[0][1] = base + 1;
  force_indices.indices[0][2] = base + 2;
  force_indices.threadsafe_derivatives_end[0]=3;
  std::size_t m=3;
  for(unsigned n=3*actiondata.numberOfNonReferenceAtoms; n<ntotal_force; ++n) {
    force_indices.indices[0][m] = n;
    ++m;
  }
  force_indices.tot_indices[0] = m;
}

}
}
#endif
