/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

   This file is part of plumed, version 2.0.

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
#ifndef __PLUMED_multicolvar_MultiColvarFunction_h
#define __PLUMED_multicolvar_MultiColvarFunction_h

#include "MultiColvarBase.h"
#include "StoreCentralAtomsVessel.h"

namespace PLMD {
namespace multicolvar {

class MultiColvarFunction : public MultiColvarBase {
private:
/// The multicolvar from which we construct these quantities
  multicolvar::MultiColvarBase* mycolv;
/// The central atom positions
  multicolvar::StoreCentralAtomsVessel* catoms;
protected:
/// Get the number of functions in the multicolvar we are operating on
  unsigned getNumberOfBaseFunctions() const;
/// Return a pointer to the multicolvar we are using as a base function
  MultiColvarBase* getPntrToMultiColvar();
/// Finish off the setup of the VectorFunction
  void completeSetup();
/// Get the position of one of the central atoms
  Vector getPositionOfCentralAtom(const unsigned&) const;
/// Add derivatives of value wrt to an atomic position 
  void addCentralAtomsDerivatives( const unsigned& , const Vector& );
/// Add derivatives of weight wrt to an atomic position
  void addCentralAtomsDerivativesOfWeight( const unsigned& , const Vector& );
/// We use the distance from mycolv as otherwise numerical derivatives dont work
  Vector getSeparation( const Vector& vec1, const Vector& vec2 ) const ;
/// Add derivative wrt to the position of the central atom
  void addDerivativeOfCentralAtomPos( const unsigned& iatom, const Tensor& der );
public:
  MultiColvarFunction(const ActionOptions&);
  static void registerKeywords( Keywords& keys );
/// Used to make sure we are calculating everything during neighbor list update step
  void unlockContributors();
  void lockContributors();
/// Active element in atoms_with_derivatives
  void atomHasDerivative( const unsigned& iatom );
/// Resize the dynamic arrays 
  void resizeDynamicArrays();
/// Regular calculate
  void calculate();
/// Calculate the numerical derivatives for this action
  void calculateNumericalDerivatives( ActionWithValue* a=NULL );
/// Do the calculation
  double doCalculation( const unsigned& j );
/// Actually compute the colvar
  virtual double compute( const unsigned& j )=0;
/// Calculate the position of the central atom
  Vector calculateCentralAtomPosition();
/// Get the position of the central atom
  virtual Vector getCentralAtom()=0;
};

inline
unsigned MultiColvarFunction::getNumberOfBaseFunctions() const {
  return mycolv->colvar_atoms.size();
} 

inline
MultiColvarBase* MultiColvarFunction::getPntrToMultiColvar(){
  return mycolv;
}

inline
Vector MultiColvarFunction::getPositionOfCentralAtom( const unsigned& iatom ) const {
  plumed_dbg_assert( iatom<colvar_atoms[current].getNumberActive() );
  return catoms->getPosition( colvar_atoms[current][iatom] );
}

inline
void MultiColvarFunction::addCentralAtomsDerivatives( const unsigned& iatom, const Vector& der ){
  plumed_dbg_assert( iatom<colvar_atoms[current].getNumberActive() );
  catoms->addAtomsDerivatives( colvar_atoms[current][iatom], der, this ); 
}

inline
void MultiColvarFunction::addCentralAtomsDerivativesOfWeight( const unsigned& iatom, const Vector& der ){
  plumed_dbg_assert( iatom<colvar_atoms[current].getNumberActive() );
  catoms->addAtomsDerivativeOfWeight( colvar_atoms[current][iatom], der, this );
}

inline
void MultiColvarFunction::atomHasDerivative( const unsigned& iatom ){
  atoms_with_derivatives.activate( iatom );
}

inline
void MultiColvarFunction::addDerivativeOfCentralAtomPos( const unsigned& iatom, const Tensor& der ){
  catoms->addDerivativeOfCentralAtomPos( colvar_atoms[current][iatom], der, this );
}


}
}
#endif
