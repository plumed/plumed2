/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013 The plumed team
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
#ifndef __PLUMED_multicolvar_StoreCentralAtomsVessel_h
#define __PLUMED_multicolvar_StoreCentralAtomsVessel_h

#include "tools/DynamicList.h"
#include "vesselbase/Vessel.h" 

namespace PLMD {
namespace multicolvar {

class MultiColvarBase;
class MultiColvarFunction;

class StoreCentralAtomsVessel : public vesselbase::Vessel {
private:
  MultiColvarBase* mycolv;
  std::vector<unsigned> start;
  unsigned nspace;
  std::vector< DynamicList<unsigned> > active_atoms;
public:
/// Constructor
  StoreCentralAtomsVessel( const vesselbase::VesselOptions& );
/// Return the number of terms
  unsigned getNumberOfTerms(){ return 2; }
/// This does the resizing of the buffer
  void resize();
/// This does nothing
  std::string description(){ return ""; }
/// This should mpi gather the active atoms
  void finish();
/// This does nothing
  bool applyForce(std::vector<double>&){ return false; }
/// This makes sure all vectors are stored
  bool calculate();
/// Get the orientation of the ith vector
  Vector getPosition( const unsigned&  ) const ;
/// Add derivatives to central atom position
  void addAtomsDerivatives( const unsigned& iatom, const Vector& df, MultiColvarFunction* funcout ) const ;
/// Add derivatives of the weight wrt to the central atom position
  void addAtomsDerivativeOfWeight( const unsigned& iatom, const Vector& df, MultiColvarFunction* funcout  ) const ; 
/// Add derivative to the central atom position
  void addDerivativeOfCentralAtomPos( const unsigned& iatom, const Tensor& df, MultiColvarFunction* funcout ) const ;
};

}
}
#endif

