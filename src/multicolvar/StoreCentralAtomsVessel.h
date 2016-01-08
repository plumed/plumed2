/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2014 The plumed team
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
#ifndef __PLUMED_multicolvar_StoreCentralAtomsVessel_h
#define __PLUMED_multicolvar_StoreCentralAtomsVessel_h

#include "vesselbase/StoreDataVessel.h" 

namespace PLMD {
namespace multicolvar {

class MultiColvarBase;
class MultiColvarFunction;

class StoreCentralAtomsVessel : public vesselbase::StoreDataVessel {
private:
/// The base multicolvar
  MultiColvarBase* mycolv;
/// A vector that is used to store derivatives
  std::vector<double> tmpdf;
public:
/// Constructor
  StoreCentralAtomsVessel( const vesselbase::VesselOptions& );
/// This does nothing
  std::string description(){ return ""; }
/// Get the orientation of the ith vector
  Vector getPosition( const unsigned& );
/// Recalculate the central atom position
  void performTask( const unsigned& );
  void finishTask( const unsigned& );
/// Get the indices
  void getIndexList( const unsigned& , const unsigned& , const unsigned& , std::vector<unsigned>& );
/// Add derivatives to central atom position
  void addAtomsDerivatives( const unsigned& iatom, const unsigned& jout, const unsigned& base_cv_no, const Vector& df, MultiColvarFunction* funcout );
};

}
}
#endif

