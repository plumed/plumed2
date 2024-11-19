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

namespace PLMD {
namespace colvar {

namespace multiColvars {
struct emptyMode {};
enum class components {
  withCompontents,
  noComponents
};
enum class plainOrScaled {
  scaled,
  plain
};
enum class scaledComponents {
  withCompontents,
  scaledComponents,
  noComponents
};
}
#define MULTICOLVAR_SETTINGS_MODE(type) using  Modetype=type
#define MULTICOLVAR_SETTINGS_SETUPF() static Modetype getModeAndSetupValues( ActionWithValue* av )
#define MULTICOLVAR_SETTINGS_CALCULATE_CONST() static void calculateCV( Modetype mode, \
                          const std::vector<double>& masses, \
                          const std::vector<double>& charges, \
                          const std::vector<Vector>& pos, \
                          std::vector<double>& vals, \
                          std::vector<std::vector<Vector> >& derivs,\
                          std::vector<Tensor>& virial, \
                          const ActionAtomistic* aa )

#define MULTICOLVAR_SETTINGS(type) MULTICOLVAR_SETTINGS_MODE(type); \
          MULTICOLVAR_SETTINGS_SETUPF(); \
          MULTICOLVAR_SETTINGS_CALCULATE_CONST()



//^no ';' here so that the macro will "consume" it and not generate a double semicolon error

template <class Colvar>
class MultiColvarTemplate : public ActionWithVector {
private:
  using Modetype =typename Colvar::Modetype;
/// An index that decides what we are calculating
  Modetype mode;
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
  void addValueWithDerivatives( const std::vector<unsigned>& shape=std::vector<unsigned>() ) override ;
  void addComponentWithDerivatives( const std::string& name, const std::vector<unsigned>& shape=std::vector<unsigned>() ) override ;
  void performTask( const unsigned&, MultiValue& ) const override ;
  void calculate() override;
};

template <class Colvar>
void MultiColvarTemplate<Colvar>::registerKeywords(Keywords& keys ) {
  Colvar::registerKeywords( keys );
  keys.add("optional","MASK","the label for a sparse matrix that should be used to determine which elements of the matrix should be computed");
  unsigned nkeys = keys.size();
  for(unsigned i=0; i<nkeys; ++i) {
    if( keys.style( keys.get(i), "atoms" ) ) keys.reset_style( keys.get(i), "numbered" );
  }
  if( keys.outputComponentExists(".#!value") ) keys.setValueDescription("the " + keys.getDisplayName() + " for each set of specified atoms");
}

template <class Colvar>
MultiColvarTemplate<Colvar>::MultiColvarTemplate(const ActionOptions&ao):
  Action(ao),
  ActionWithVector(ao),
  mode(),
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
      Colvar::parseAtomList( i, t, this );
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
  mode = Colvar::getModeAndSetupValues( this );
}

template <class Colvar>
unsigned MultiColvarTemplate<Colvar>::getNumberOfDerivatives() {
  return 3*getNumberOfAtoms()+9;
}

template <class Colvar>
void MultiColvarTemplate<Colvar>::calculate() {
  if( wholemolecules ) makeWhole();
  runAllTasks();
}

template <class Colvar>
void MultiColvarTemplate<Colvar>::addValueWithDerivatives( const std::vector<unsigned>& shape ) {
  std::vector<unsigned> s(1); s[0]=ablocks[0].size(); addValue( s );
}

template <class Colvar>
void MultiColvarTemplate<Colvar>::addComponentWithDerivatives( const std::string& name, const std::vector<unsigned>& shape ) {
  std::vector<unsigned> s(1); s[0]=ablocks[0].size(); addComponent( name, s );
}

template <class Colvar>
void MultiColvarTemplate<Colvar>::performTask( const unsigned& task_index, MultiValue& myvals ) const {
  // Retrieve the positions
  std::vector<Vector> & fpositions( myvals.getFirstAtomVector() );
  if( fpositions.size()!=ablocks.size() ) {
    fpositions.resize( ablocks.size() );
  }
  for(unsigned i=0; i<ablocks.size(); ++i) {
    fpositions[i] = getPosition( ablocks[i][task_index] );
  }
  // If we are using pbc make whole
  if( usepbc ) {
    if( fpositions.size()==1 ) {
      fpositions[0]=pbcDistance(Vector(0.0,0.0,0.0),getPosition( ablocks[0][task_index] ) );
    } else {
      for(unsigned j=0; j<fpositions.size()-1; ++j) {
        const Vector & first (fpositions[j]); Vector & second (fpositions[j+1]);
        second=first+pbcDistance(first,second);
      }
    }
  } else if( fpositions.size()==1 ) {
    fpositions[0]=delta(Vector(0.0,0.0,0.0),getPosition( ablocks[0][task_index] ) );
  }
  // Retrieve the masses and charges
  myvals.resizeTemporyVector(2);
  std::vector<double> & mass( myvals.getTemporyVector(0) );
  std::vector<double> & charge( myvals.getTemporyVector(1) );
  if( mass.size()!=ablocks.size() ) {
    mass.resize(ablocks.size());
    charge.resize(ablocks.size());
  }
  for(unsigned i=0; i<ablocks.size(); ++i) {
    mass[i]=getMass( ablocks[i][task_index] );
    charge[i]=getCharge( ablocks[i][task_index] );
  }
  // Make some space to store various things
  std::vector<double> values( getNumberOfComponents() );
  std::vector<Tensor> & virial( myvals.getFirstAtomVirialVector() );
  std::vector<std::vector<Vector> > & derivs( myvals.getFirstAtomDerivativeVector() );
  if( derivs.size()!=values.size() ) {
    derivs.resize( values.size() );
    virial.resize( values.size() );
  }
  for(unsigned i=0; i<derivs.size(); ++i) {
    if( derivs[i].size()<ablocks.size() ) {
      derivs[i].resize( ablocks.size() );
    }
  }
  // Calculate the CVs using the method in the Colvar
  Colvar::calculateCV( mode, mass, charge, fpositions, values, derivs, virial, this );
  for(unsigned i=0; i<values.size(); ++i) myvals.setValue( i, values[i] );
  // Finish if there are no derivatives
  if( doNotCalculateDerivatives() ) return;

  // Now transfer the derivatives to the underlying MultiValue
  for(unsigned i=0; i<ablocks.size(); ++i) {
    unsigned base=3*ablocks[i][task_index];
    for(int j=0; j<getNumberOfComponents(); ++j) {
      myvals.addDerivative( j, base + 0, derivs[j][i][0] );
      myvals.addDerivative( j, base + 1, derivs[j][i][1] );
      myvals.addDerivative( j, base + 2, derivs[j][i][2] );
    }
    // Check for duplicated indices during update to avoid double counting
    bool newi=true;
    for(unsigned j=0; j<i; ++j) {
      if( ablocks[j][task_index]==ablocks[i][task_index] ) { newi=false; break; }
    }
    if( !newi ) continue;
    for(int j=0; j<getNumberOfComponents(); ++j) {
      myvals.updateIndex( j, base );
      myvals.updateIndex( j, base + 1 );
      myvals.updateIndex( j, base + 2 );
    }
  }
  unsigned nvir=3*getNumberOfAtoms();
  for(int j=0; j<getNumberOfComponents(); ++j) {
    for(unsigned i=0; i<3; ++i) {
      for(unsigned k=0; k<3; ++k) {
        myvals.addDerivative( j, nvir + 3*i + k, virial[j][i][k] );
        myvals.updateIndex( j, nvir + 3*i + k );
      }
    }
  }
}

}
}
#endif
