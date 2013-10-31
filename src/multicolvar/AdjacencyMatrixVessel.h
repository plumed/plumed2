/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

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
#ifndef __PLUMED_multicolvar_AdjacencyMatrixVessel_h
#define __PLUMED_multicolvar_AdjacencyMatrixVessel_h

#include "vesselbase/StoreDataVessel.h" 

namespace PLMD {
namespace multicolvar {

class AdjacencyMatrixAction;

class AdjacencyMatrixVessel : public vesselbase::StoreDataVessel {
friend class VectorMultiColvar;
private:
  unsigned nrows;
/// Pointer to underlying action
  AdjacencyMatrixAction* function;
/// Tempory vector for chain rule
  std::vector<double> tmpdf;
public:
  static void registerKeywords( Keywords& keys );
/// Constructor
  AdjacencyMatrixVessel( const vesselbase::VesselOptions& );
/// This does nothing
  std::string description(){ return ""; }
/// This recomputes the colvar
  void recompute( const unsigned& ivec, const unsigned& jstore );
/// Get the i,j th element of the matrix
  double getElement( const unsigned& ivec ); 
/// Finish the calculation
  void finish();
};

inline
double AdjacencyMatrixVessel::getElement( const unsigned& ivec ){
  return getComponent( ivec, 0 );
} 

}
}
#endif

