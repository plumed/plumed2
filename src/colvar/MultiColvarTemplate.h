/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2023 The plumed team
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
#ifndef __PLUMED_colvar_MultiColvarTemplate_h
#define __PLUMED_colvar_MultiColvarTemplate_h

#include "core/ActionWithVector.h"
#include "core/ParallelTaskManager.h"

namespace PLMD {

class Colvar;

namespace colvar {

class MultiColvarInput {
public:
  bool usepbc;
  unsigned mode;
  MultiColvarInput() : usepbc(false), mode(0)  {}
  MultiColvarInput( bool u, const unsigned m ) : usepbc(u), mode(m) {}
  MultiColvarInput& operator=( const MultiColvarInput& m ) { usepbc = m.usepbc; mode = m.mode; return *this; }
};

struct ColvarInput {
  unsigned mode;
  const Pbc& pbc;
  View2D<const double,helpers::dynamic_extent,3> pos;
  View<const double,helpers::dynamic_extent> mass;
  View<const double,helpers::dynamic_extent> charges;
  ColvarInput( unsigned m, unsigned natoms, const double* p, const double* w, const double* q, const Pbc& box );
  static ColvarInput createColvarInput( unsigned m, const std::vector<Vector>& p, const Colvar* colv );
};

class ColvarOutput {
private:
  std::size_t ncomponents;
  class DerivHelper {
  private:
    std::size_t nderiv;
    std::vector<double>& derivatives;
  public:
    DerivHelper( std::size_t n, std::vector<double>& d ) : nderiv(n), derivatives(d) {}
    View2D<double, helpers::dynamic_extent, 3> operator[](std::size_t i) {
      return View2D<double, helpers::dynamic_extent, 3>( derivatives.data() + i*nderiv, nderiv );
    }
    Vector getAtomDerivatives( std::size_t i, std::size_t a ) {
      std::size_t base = i*nderiv + 3*a;
      return Vector( derivatives[base], derivatives[base+1], derivatives[base+2] );
    }
  };
  class VirialHelper {
  private:
    std::size_t nderiv;
    std::vector<double>& derivatives;
  public:
    VirialHelper( std::size_t n, std::vector<double>& d ) : nderiv(n), derivatives(d) {}
    Tensor operator[](std::size_t i) const {
      std::size_t n=(i+1)*nderiv;
      return Tensor( derivatives[n-9], derivatives[n-8], derivatives[n-7], derivatives[n-6], derivatives[n-5], derivatives[n-4], derivatives[n-3], derivatives[n-2], derivatives[n-1] );
    }
    void set( std::size_t i, const Tensor& v ) {
      std::size_t n=(i+1)*nderiv;
      derivatives[n-9]=v[0][0]; derivatives[n-8]=v[0][1]; derivatives[n-7]=v[0][2];
      derivatives[n-6]=v[1][0]; derivatives[n-5]=v[1][1]; derivatives[n-4]=v[1][2];
      derivatives[n-3]=v[2][0]; derivatives[n-2]=v[2][1]; derivatives[n-1]=v[2][2];
    }
  };
public:
  View<double,helpers::dynamic_extent> values;
  DerivHelper derivs;
  VirialHelper virial;
  ColvarOutput( std::size_t n, double* v, std::size_t m, std::vector<double>& d );
  Vector getAtomDerivatives( std::size_t i, std::size_t a ) { return derivs.getAtomDerivatives(i,a); }
  void setBoxDerivativesNoPbc( const ColvarInput& inpt );
  static ColvarOutput createColvarOutput( std::vector<double>& v, std::vector<double>& d, Colvar* action );
};

template <class T>
class MultiColvarTemplate : public ActionWithVector {
private:
/// The parallel task manager
  ParallelTaskManager<MultiColvarTemplate<T>,MultiColvarInput> taskmanager;
/// An index that decides what we are calculating
  unsigned mode;
/// Are we using pbc to calculate the CVs
  bool usepbc;
/// Do we reassemble the molecule
  bool wholemolecules;
/// The number of atoms per task
  unsigned natoms_per_task;
public:
  static void registerKeywords(Keywords&);
  explicit MultiColvarTemplate(const ActionOptions&);
  unsigned getNumberOfDerivatives() override ;
  void addValueWithDerivatives( const std::vector<unsigned>& shape=std::vector<unsigned>() ) override ;
  void addComponentWithDerivatives( const std::string& name, const std::vector<unsigned>& shape=std::vector<unsigned>() ) override ;
  void getInputData( std::vector<double>& inputdata ) const override ;
  void performTask( const unsigned&, MultiValue& ) const override { plumed_error(); }
  void calculate() override;
  void applyNonZeroRankForces( std::vector<double>& outforces ) override ;
  static void performTask( unsigned task_index, ParallelActionsInput<MultiColvarInput>& input, ParallelActionsOutput& output );
  static void gatherForces( unsigned task_index, const ParallelActionsInput<MultiColvarInput>& input, const ForceInput& fdata, ForceOutput& forces );
  static void gatherThreads( ForceOutput& forces );
  static void transferToValue( unsigned task_index, const std::vector<double>& values, std::vector<double>& value_mat );
};

template <class T>
void MultiColvarTemplate<T>::registerKeywords(Keywords& keys ) {
  T::registerKeywords( keys );
  ParallelTaskManager<MultiColvarTemplate<T>,MultiColvarInput>::registerKeywords( keys );
  keys.add("optional","MASK","the label for a sparse matrix that should be used to determine which elements of the matrix should be computed");
  unsigned nkeys = keys.size();
  for(unsigned i=0; i<nkeys; ++i) {
    if( keys.style( keys.get(i), "atoms" ) ) keys.reset_style( keys.get(i), "numbered" );
  }
  if( keys.outputComponentExists(".#!value") ) keys.setValueDescription("the " + keys.getDisplayName() + " for each set of specified atoms");
}

template <class T>
MultiColvarTemplate<T>::MultiColvarTemplate(const ActionOptions&ao):
  Action(ao),
  ActionWithVector(ao),
  taskmanager(this),
  mode(0),
  usepbc(true),
  wholemolecules(false)
{
  std::vector<AtomNumber> all_atoms;
  if( getName()=="POSITION_VECTOR" || getName()=="MASS_VECTOR" || getName()=="CHARGE_VECTOR" ) parseAtomList( "ATOMS", all_atoms );
  if( all_atoms.size()>0 ) {
    natoms_per_task=1;
  } else {
    std::vector<AtomNumber> t;
    for(int i=1;; ++i ) {
      T::parseAtomList( i, t, this );
      if( t.empty() ) break;

      if( i==1 ) { natoms_per_task=t.size(); }
      if( t.size()!=natoms_per_task ) {
        std::string ss; Tools::convert(i,ss);
        error("ATOMS" + ss + " keyword has the wrong number of atoms");
      }
      for(unsigned j=0; j<natoms_per_task; ++j) all_atoms.push_back( t[j] );
      t.resize(0);
    }
  }
  if( all_atoms.size()==0 ) error("No atoms have been specified");
  requestAtoms(all_atoms);
  if( keywords.exists("NOPBC") ) {
    bool nopbc=!usepbc; parseFlag("NOPBC",nopbc);
    usepbc=!nopbc;
  }
  if( keywords.exists("WHOLEMOLECULES") ) {
    parseFlag("WHOLEMOLECULES",wholemolecules);
    if( wholemolecules ) usepbc=false;
  }
  if( usepbc ) log.printf("  using periodic boundary conditions\n");
  else    log.printf("  without periodic boundary conditions\n");

  // Setup the values
  mode = T::getModeAndSetupValues( this );
  // This sets up an array in the parallel task manager to hold all the indices
  // Sets up the index list in the task manager
  taskmanager.setNumberOfIndicesAndDerivativesPerTask( natoms_per_task, 3*natoms_per_task + 9 );
  taskmanager.setNumberOfThreadedForces( 9 ); taskmanager.setActionInput( MultiColvarInput( usepbc, mode ) );
}

template <class T>
unsigned MultiColvarTemplate<T>::getNumberOfDerivatives() {
  return 3*getNumberOfAtoms()+9;
}

template <class T>
void MultiColvarTemplate<T>::calculate() {
  if( wholemolecules ) makeWhole();
  taskmanager.runAllTasks();
}

template <class T>
void MultiColvarTemplate<T>::applyNonZeroRankForces( std::vector<double>& outforces ) {
  taskmanager.applyForces( outforces );
}

template <class T>
void MultiColvarTemplate<T>::addValueWithDerivatives( const std::vector<unsigned>& shape ) {
  std::vector<unsigned> s(1); s[0]=getNumberOfAtoms() / natoms_per_task; addValue( s );
}

template <class T>
void MultiColvarTemplate<T>::addComponentWithDerivatives( const std::string& name, const std::vector<unsigned>& shape ) {
  std::vector<unsigned> s(1); s[0]=getNumberOfAtoms() / natoms_per_task; addComponent( name, s );
}

template <class T>
void MultiColvarTemplate<T>::getInputData( std::vector<double>& inputdata ) const {
  unsigned ntasks = getConstPntrToComponent(0)->getNumberOfStoredValues();
  if( inputdata.size()!=5*natoms_per_task*ntasks ) inputdata.resize( 5*natoms_per_task*ntasks );

  std::size_t k=0;
  for(unsigned i=0; i<ntasks; ++i) {
    for(unsigned j=0; j<natoms_per_task; ++j) {
      Vector mypos( getPosition( natoms_per_task*i + j ) );
      inputdata[k] = mypos[0]; k++;
      inputdata[k] = mypos[1]; k++;
      inputdata[k] = mypos[2]; k++;
    }
    for(unsigned j=0; j<natoms_per_task; ++j) { inputdata[k] = getMass( natoms_per_task*i + j ); k++; }
    for(unsigned j=0; j<natoms_per_task; ++j) { inputdata[k] = getCharge( natoms_per_task*i + j ); k++; }
  }
}

template <class T>
void MultiColvarTemplate<T>::performTask( unsigned task_index, ParallelActionsInput<MultiColvarInput>& input, ParallelActionsOutput& output ) {
  std::size_t pos_start = 5*input.nindices_per_task*task_index;
  if( input.actiondata.usepbc ) {
    if( input.nindices_per_task==1 ) {
      Vector fpos=input.pbc.distance(Vector(0.0,0.0,0.0),Vector(input.inputdata[pos_start], input.inputdata[pos_start+1], input.inputdata[pos_start+2]) );
      input.inputdata[pos_start]=fpos[0]; input.inputdata[pos_start+1]=fpos[1]; input.inputdata[pos_start+2]=fpos[2];
    } else {
      std::size_t apos_start = pos_start;
      for(unsigned j=0; j<input.nindices_per_task-1; ++j) {
        Vector first(input.inputdata[apos_start], input.inputdata[apos_start+1], input.inputdata[apos_start+2]);
        Vector second(input.inputdata[apos_start+3], input.inputdata[apos_start+4], input.inputdata[apos_start+5]);
        second=first+input.pbc.distance(first,second);
        input.inputdata[apos_start+3]=second[0]; input.inputdata[apos_start+4]=second[1]; input.inputdata[apos_start+5]=second[2];
        apos_start += 3;
      }
    }
  } else if( input.nindices_per_task==1 ) {
    Vector fpos=delta(Vector(0.0,0.0,0.0),Vector(input.inputdata[pos_start], input.inputdata[pos_start+1], input.inputdata[pos_start+2]));
    input.inputdata[pos_start]=fpos[0]; input.inputdata[pos_start+1]=fpos[1]; input.inputdata[pos_start+2]=fpos[2];
  }

  std::size_t mass_start = pos_start + 3*input.nindices_per_task;
  std::size_t charge_start = mass_start + input.nindices_per_task;
  ColvarOutput cvout = ColvarOutput( input.ncomponents, output.values.data(), 3*input.nindices_per_task+9, output.derivatives );
  T::calculateCV( ColvarInput(input.actiondata.mode, input.nindices_per_task, input.inputdata.data()+pos_start, input.inputdata.data()+mass_start, input.inputdata.data()+charge_start, input.pbc), cvout );
}

template <class T>
void MultiColvarTemplate<T>::transferToValue( unsigned task_index, const std::vector<double>& values, std::vector<double>& value_mat ) {
  unsigned ncomponents = values.size();
  for(unsigned i=0; i<ncomponents; ++i) value_mat[ncomponents*task_index+i] = values[i];
}

template <class T>
void MultiColvarTemplate<T>::gatherForces( unsigned task_index, const ParallelActionsInput<MultiColvarInput>& input, const ForceInput& fdata, ForceOutput& forces ) {
  std::size_t base = 3*task_index*input.nindices_per_task;
  for(unsigned i=0; i<input.ncomponents; ++i) {
    unsigned m = 0; double ff = fdata.force[i];
    for(unsigned j=0; j<input.nindices_per_task; ++j) {
      forces.thread_unsafe[base + m] += ff*fdata.deriv[i][m]; m++;
      forces.thread_unsafe[base + m] += ff*fdata.deriv[i][m]; m++;
      forces.thread_unsafe[base + m] += ff*fdata.deriv[i][m]; m++;
    }
    for(unsigned n=0; n<9; ++n) { forces.thread_safe[n] += ff*fdata.deriv[i][m]; m++; }
  }
}

template <class T>
void MultiColvarTemplate<T>::gatherThreads( ForceOutput& forces ) {
  unsigned k=0; for(unsigned n=forces.thread_unsafe.size()-9; n<forces.thread_unsafe.size(); ++n) { forces.thread_unsafe[n] += forces.thread_safe[k]; k++; }
}

}
}
#endif
