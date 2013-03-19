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
#ifndef __PLUMED_multicolvar_StoreCentralAtomsVessel_h
#define __PLUMED_multicolvar_StoreCentralAtomsVessel_h

#include "tools/DynamicList.h"
#include "vesselbase/Vessel.h" 
#include "MultiColvar.h"

namespace PLMD {
namespace multicolvar {

class MultiColvar;

class StoreCentralAtomsVessel : public vesselbase::Vessel {
private:
  MultiColvar* mycolv;
  std::vector<unsigned> start;
  bool wasforced;
  std::vector< DynamicList<unsigned> > active_der;
  std::vector<double> forces;
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
/// Add some force to the atoms
  void addForces( const std::vector<double>& );
/// Copy forces to local arrays
  bool applyForce(std::vector<double>&);
/// Set all forces equal to zero before main loop
  void prepare();
/// This makes sure all vectors are stored
  bool calculate();
/// Get the orientation of the ith vector
  Vector getPosition( const unsigned&  ) const ;
/// Chain rule for central atom
  void chainRuleForCentralAtom( const unsigned& iatom, const unsigned& iderno, const Vector& df, vesselbase::ActionWithVessel* act ) const ;
};

}
}
#endif

