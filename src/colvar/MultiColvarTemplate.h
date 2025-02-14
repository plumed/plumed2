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
  MultiColvarInput( const bool& u, const unsigned& m ) : usepbc(u), mode(m) {}
  MultiColvarInput& operator=( const MultiColvarInput& m ) { usepbc = m.usepbc; mode = m.mode; return *this; }
};

class ColvarInput {
public:
  unsigned mode;
  const Pbc& pbc;
  const std::vector<Vector>& pos;
  const std::vector<double>& mass;
  const std::vector<double>& charges;
  ColvarInput( const unsigned& m, const std::vector<Vector>& p, const std::vector<double>& w, const std::vector<double>& q, const Pbc& box );
  static ColvarInput createColvarInput( const unsigned& m, const std::vector<Vector>& p, const Colvar* colv );
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
  static void performTask( unsigned task_index, const ParallelActionsInput<MultiColvarInput>& input, ParallelActionsOutput& output );
  static void gatherForces( unsigned task_index, const ParallelActionsInput<MultiColvarInput>& input, const Matrix<double>& force_in, const Matrix<double>& derivs, std::vector<double>& thred_unsafe_force_out, std::vector<double>& thred_safe_force_out );
  static void gatherThreads( std::vector<double>& thred_unsafe_force_out, std::vector<double>& thred_safe_force_out );
  static void transferToValue( unsigned task_index, const std::vector<double>& values, Matrix<double>& value_mat );
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
  taskmanager.setActionInput( MultiColvarInput( usepbc, mode ) );
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
      inputdata[k] = mypos[0];
      inputdata[k+1] = mypos[1];
      inputdata[k+2] = mypos[2];
      inputdata[k+3] = getMass( natoms_per_task*i + j );
      inputdata[k+4] = getCharge( natoms_per_task*i + j );
      k+=5;
    }
  }
}

template <class T>
void MultiColvarTemplate<T>::performTask( unsigned task_index, const ParallelActionsInput<MultiColvarInput>& input, ParallelActionsOutput& output ) {
  std::vector<double> mass( input.nindices_per_task );
  std::vector<double> charge( input.nindices_per_task );
  std::vector<Vector> fpositions( input.nindices_per_task );
  for(unsigned i=0; i<fpositions.size(); ++i) {
    std::size_t base = 5*fpositions.size()*task_index + 5*i;
    fpositions[i][0] = input.inputdata[base + 0];
    fpositions[i][1] = input.inputdata[base + 1];
    fpositions[i][2] = input.inputdata[base + 2];
    mass[i] = input.inputdata[base + 3];
    charge[i] = input.inputdata[base + 4];
  }
  if( input.actiondata.usepbc ) {
    if( fpositions.size()==1 ) {
      fpositions[0]=input.pbc.distance(Vector(0.0,0.0,0.0),fpositions[0]);
    } else {
      for(unsigned j=0; j<fpositions.size()-1; ++j) {
        const Vector & first (fpositions[j]); Vector & second (fpositions[j+1]);
        second=first+input.pbc.distance(first,second);
      }
    }
  } else if( fpositions.size()==1 ) fpositions[0]=delta(Vector(0.0,0.0,0.0),fpositions[0]);

  std::vector<Tensor> virial( input.ncomponents );
  Matrix<Vector> derivs( input.ncomponents, fpositions.size() );
  T::calculateCV( ColvarInput( input.actiondata.mode, fpositions, mass, charge, input.pbc ), output.values, derivs, virial );
  if( input.noderiv ) return;

  for(unsigned i=0; i<input.ncomponents; ++i) {
    unsigned k=0;
    for(unsigned j=0; j<fpositions.size(); ++j) {
      output.derivatives[i][k] = derivs[i][j][0]; k++;
      output.derivatives[i][k] = derivs[i][j][1]; k++;
      output.derivatives[i][k] = derivs[i][j][2]; k++;
    }
    for(unsigned j=0; j<3; ++j) {
      for(unsigned n=0; n<3; ++n) { output.derivatives[i][k] = virial[i][j][n]; k++; }
    }
  }
}

template <class T>
void MultiColvarTemplate<T>::transferToValue( unsigned task_index, const std::vector<double>& values, Matrix<double>& value_mat ) {
  for(unsigned i=0; i<values.size(); ++i) value_mat[task_index][i] = values[i];
}

template <class T>
void MultiColvarTemplate<T>::gatherForces( unsigned task_index, const ParallelActionsInput<MultiColvarInput>& input, const Matrix<double>& force_in, const Matrix<double>& derivs, std::vector<double>& thred_unsafe_force_out, std::vector<double>& thred_safe_force_out ) {
  std::size_t base = 3*task_index*input.nindices_per_task;
  for(unsigned i=0; i<force_in.ncols(); ++i) {
    unsigned m = 0; double ff = force_in[task_index][i];
    for(unsigned j=0; j<input.nindices_per_task; ++j) {
      thred_unsafe_force_out[base + m] += ff*derivs[i][m]; m++;
      thred_unsafe_force_out[base + m] += ff*derivs[i][m]; m++;
      thred_unsafe_force_out[base + m] += ff*derivs[i][m]; m++;
    }
    for(unsigned n=thred_safe_force_out.size()-9; n<thred_safe_force_out.size(); ++n) { thred_safe_force_out[n] += ff*derivs[i][m]; m++; }
  }
}

template <class T>
void MultiColvarTemplate<T>::gatherThreads( std::vector<double>& thred_safe_force_out, std::vector<double>& thred_unsafe_force_out ) {
  for(unsigned n=thred_safe_force_out.size()-9; n<thred_safe_force_out.size(); ++n) thred_unsafe_force_out[n] += thred_safe_force_out[n];
}

}
}
#endif
