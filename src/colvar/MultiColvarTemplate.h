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
  ParallelTaskManager<MultiColvarTemplate<T> > taskmanager;
/// An index that decides what we are calculating
  unsigned mode;
/// Are we using pbc to calculate the CVs
  bool usepbc;
/// Do we reassemble the molecule
  bool wholemolecules;
/// Blocks of atom numbers
  std::vector< std::vector<unsigned> > ablocks;
public:
  static void registerKeywords(Keywords&);
  explicit MultiColvarTemplate(const ActionOptions&);
  unsigned getNumberOfDerivatives() override ;
  unsigned getNumberOfAtomsPerTask() const override ;
  void addValueWithDerivatives( const std::vector<unsigned>& shape=std::vector<unsigned>() ) override ;
  void addComponentWithDerivatives( const std::string& name, const std::vector<unsigned>& shape=std::vector<unsigned>() ) override ;
  void getInputData( std::vector<double>& inputdata ) const override ;
  void performTask( const unsigned&, MultiValue& ) const override ;
  void calculate() override;
  static std::pair<std::vector<double>,Matrix<double> > performTask( const unsigned& task_index, const ParallelActionsInput& input );
  static void performTask( const unsigned& m, const std::vector<std::size_t>& der_indices, const bool noderiv, const bool haspbc, const Pbc& pbc, MultiValue& myvals );
  static void transferToValue( const unsigned& task_index, const std::vector<double>& values, Matrix<double>& value_mat );
};

template <class T>
void MultiColvarTemplate<T>::registerKeywords(Keywords& keys ) {
  T::registerKeywords( keys );
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
    ablocks.resize(1); ablocks[0].resize( all_atoms.size() );
    for(unsigned i=0; i<all_atoms.size(); ++i) ablocks[0][i] = i;
  } else {
    std::vector<AtomNumber> t;
    for(int i=1;; ++i ) {
      T::parseAtomList( i, t, this );
      if( t.empty() ) break;

      if( i==1 ) { ablocks.resize(t.size()); }
      if( t.size()!=ablocks.size() ) {
        std::string ss; Tools::convert(i,ss);
        error("ATOMS" + ss + " keyword has the wrong number of atoms");
      }
      for(unsigned j=0; j<ablocks.size(); ++j) {
        ablocks[j].push_back( ablocks.size()*(i-1)+j ); all_atoms.push_back( t[j] );
      }
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
  std::vector<std::size_t> ind( ablocks.size()*ablocks[0].size() );
  for(unsigned i=0; i<ablocks[0].size(); ++i) {
    for(unsigned j=0; j<ablocks.size(); ++j) ind[i*ablocks.size() + j] = ablocks[j][i];
  }
  // Sets up the index list in the task manager
  taskmanager.setupIndexList( ablocks.size(), ind );
  taskmanager.setPbcFlag( usepbc );
  taskmanager.setMode( mode );
}

template <class T>
unsigned MultiColvarTemplate<T>::getNumberOfDerivatives() {
  return 3*getNumberOfAtoms()+9;
}

template <class T>
void MultiColvarTemplate<T>::calculate() {
  if( wholemolecules ) makeWhole();
  setForwardPass(true);
  taskmanager.runAllTasks( ablocks.size() );
  setForwardPass(false);
}

template <class T>
void MultiColvarTemplate<T>::addValueWithDerivatives( const std::vector<unsigned>& shape ) {
  std::vector<unsigned> s(1); s[0]=ablocks[0].size(); addValue( s );
}

template <class T>
void MultiColvarTemplate<T>::addComponentWithDerivatives( const std::string& name, const std::vector<unsigned>& shape ) {
  std::vector<unsigned> s(1); s[0]=ablocks[0].size(); addComponent( name, s );
}

template <class T>
unsigned MultiColvarTemplate<T>::getNumberOfAtomsPerTask() const {
  return ablocks.size();
}

template <class T>
void MultiColvarTemplate<T>::getInputData( std::vector<double>& inputdata ) const {
  unsigned ntasks = ablocks[0].size(); std::size_t k=0;
  if( inputdata.size()!=5*ablocks.size()*ntasks ) inputdata.resize( 5*ablocks.size()*ntasks );
  for(unsigned i=0; i<ntasks; ++i) {
    for(unsigned j=0; j<ablocks.size(); ++j) {
      Vector mypos( getPosition( ablocks[j][i] ) );
      inputdata[k] = mypos[0];
      inputdata[k+1] = mypos[1];
      inputdata[k+2] = mypos[2];
      inputdata[k+3] = getMass( ablocks[j][i] );
      inputdata[k+4] = getCharge( ablocks[j][i] );
      k+=5;
    }
  }
}

template <class T>
void MultiColvarTemplate<T>::performTask( const unsigned& task_index, MultiValue& myvals ) const {
  // Retrieve the positions
  std::vector<double> & mass( myvals.getTemporyVector(0) );
  std::vector<double> & charge( myvals.getTemporyVector(1) );
  std::vector<Vector> & fpositions( myvals.getFirstAtomVector() );
  for(unsigned i=0; i<ablocks.size(); ++i) {
    fpositions[i] = getPosition( ablocks[i][task_index] );
    mass[i]=getMass( ablocks[i][task_index] );
    charge[i]=getCharge( ablocks[i][task_index] );
  }
  std::vector<std::size_t> der_indices( ablocks.size() );
  for(unsigned i=0; i<der_indices.size(); ++i) der_indices[i] = ablocks[i][task_index];
  performTask( mode, der_indices, doNotCalculateDerivatives(), usepbc, getPbc(), myvals );
}

template <class T>
std::pair<std::vector<double>,Matrix<double> > MultiColvarTemplate<T>::performTask( const unsigned& task_index, const ParallelActionsInput& input ) {
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
  if( input.usepbc ) {
    if( fpositions.size()==1 ) {
      fpositions[0]=input.pbc.distance(Vector(0.0,0.0,0.0),fpositions[0]);
    } else {
      for(unsigned j=0; j<fpositions.size()-1; ++j) {
        const Vector & first (fpositions[j]); Vector & second (fpositions[j+1]);
        second=first+input.pbc.distance(first,second);
      }
    }
  } else if( fpositions.size()==1 ) fpositions[0]=delta(Vector(0.0,0.0,0.0),fpositions[0]);
  std::vector<double> values( input.ncomponents );
  std::vector<Tensor> virial( input.ncomponents );
  Matrix<Vector> derivs( values.size(), fpositions.size() );
  Matrix<double> derivatives( values.size(), 3*fpositions.size() + 9 );
  T::calculateCV( ColvarInput( input.mode, fpositions, mass, charge, input.pbc ), values, derivs, virial );
  if( input.noderiv ) return {values, derivatives};
}

template <class T>
void MultiColvarTemplate<T>::performTask( const unsigned& m, const std::vector<std::size_t>& der_indices, const bool noderiv, const bool haspbc, const Pbc& pbc, MultiValue& myvals ) {
  // Retrieve the inputs
  std::vector<double> & mass( myvals.getTemporyVector(0) );
  std::vector<double> & charge( myvals.getTemporyVector(1) );
  std::vector<Vector> & fpositions( myvals.getFirstAtomVector() );
  // If we are using pbc make whole
  if( haspbc ) {
    if( fpositions.size()==1 ) {
      fpositions[0]=pbc.distance(Vector(0.0,0.0,0.0),fpositions[0]);
    } else {
      for(unsigned j=0; j<fpositions.size()-1; ++j) {
        const Vector & first (fpositions[j]); Vector & second (fpositions[j+1]);
        second=first+pbc.distance(first,second);
      }
    }
  } else if( fpositions.size()==1 ) fpositions[0]=delta(Vector(0.0,0.0,0.0),fpositions[0]);
  // Make some space to store various things
  std::vector<double> values( myvals.getNumberOfValues() );
  std::vector<Tensor> & virial( myvals.getFirstAtomVirialVector() );
  Matrix<Vector> derivs( values.size(), fpositions.size() );
  // Calculate the CVs using the method in the Colvar
  T::calculateCV( ColvarInput(m, fpositions, mass, charge, pbc ), values, derivs, virial );
  for(unsigned i=0; i<values.size(); ++i) myvals.setValue( i, values[i] );
  // Finish if there are no derivatives
  if( noderiv ) return;

  // Now transfer the derivatives to the underlying MultiValue
  for(unsigned i=0; i<der_indices.size(); ++i) {
    unsigned base=3*der_indices[i];
    for(int j=0; j<values.size(); ++j) {
      myvals.addDerivative( j, base + 0, derivs[j][i][0] );
      myvals.addDerivative( j, base + 1, derivs[j][i][1] );
      myvals.addDerivative( j, base + 2, derivs[j][i][2] );
    }
    // Check for duplicated indices during update to avoid double counting
    bool newi=true;
    for(unsigned j=0; j<i; ++j) {
      if( der_indices[j]==der_indices[i] ) { newi=false; break; }
    }
    if( !newi ) continue;
    for(int j=0; j<values.size(); ++j) {
      myvals.updateIndex( j, base );
      myvals.updateIndex( j, base + 1 );
      myvals.updateIndex( j, base + 2 );
    }
  }
  unsigned nvir=myvals.getNumberOfDerivatives() - 9;
  for(int j=0; j<values.size(); ++j) {
    for(unsigned i=0; i<3; ++i) {
      for(unsigned k=0; k<3; ++k) {
        myvals.addDerivative( j, nvir + 3*i + k, virial[j][i][k] );
        myvals.updateIndex( j, nvir + 3*i + k );
      }
    }
  }
}

template <class T>
void MultiColvarTemplate<T>::transferToValue( const unsigned& task_index, const std::vector<double>& values, Matrix<double>& value_mat ) {
  for(unsigned i=0; i<values.size(); ++i) value_mat[task_index][i] = values[i];
}

}
}
#endif
