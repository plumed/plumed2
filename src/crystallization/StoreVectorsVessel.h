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
#ifndef __PLUMED_crystallization_StoreVectorsVessel_h
#define __PLUMED_crystallization_StoreVectorsVessel_h

#include "tools/DynamicList.h"
#include "vesselbase/StoreDataVessel.h" 

namespace PLMD {

class multicolvar::MultiColvarFunction;

namespace crystallization {

class VectorMultiColvar;

class StoreVectorsVessel : public vesselbase::StoreDataVessel {
friend class VectorMultiColvar;
private:
/// We want to store the director rather than the value
  bool store_director; 
  unsigned ncomponents;
  std::vector<double> myfvec;
  VectorMultiColvar* vecs;
  void normalizeVector( const int& );
public:
  static void registerKeywords( Keywords& keys );
/// Constructor
  StoreVectorsVessel( const vesselbase::VesselOptions& );
/// This turns on the full use of this action for storage
  void usedInFunction( const bool& );
/// This makes sure vectors are normalized (they are already stored)
  bool calculate();
/// This reperforms a calculation
  void recompute( const unsigned& , const unsigned& );
/// This does nothing
  std::string description(){ return ""; }
/// Get the orientation of the ith vector
  void getVector( const unsigned& , std::vector<double>& );
/// Chain rule for component
  void chainRuleForComponent( const unsigned& , const unsigned& , const unsigned& jout, const unsigned& , const double& , multicolvar::MultiColvarFunction* );
/// Chain rule for whole vector
  void chainRuleForVector( const unsigned& , const unsigned& , const unsigned& , const std::vector<double>& , multicolvar::MultiColvarFunction* );
};

inline
void StoreVectorsVessel::getVector( const unsigned& imol, std::vector<double>& vec ){
  plumed_dbg_assert( vec.size()==getNumberOfComponents() );
  for(unsigned i=0;i<getNumberOfComponents();++i) vec[i]=getComponent( imol, i );
}



}
}
#endif

