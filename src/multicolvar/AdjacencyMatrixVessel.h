/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2015 The plumed team
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
/// Pointer to underlying action
  AdjacencyMatrixAction* function;
/// Has the vessel been finished
  bool finished;
public:
  static void registerKeywords( Keywords& keys );
/// Constructor
  explicit AdjacencyMatrixVessel( const vesselbase::VesselOptions& );
/// Ensures we use less memory for buffer in final loop
  void setBufferStart( unsigned& start );
/// Ensures that finish is set properly
  void prepare();
/// Set the finished flag true
  void setFinishedTrue();
/// An overwrite of calculate to stop this being done more than once
  bool calculate( const unsigned& current, MultiValue& myvals, std::vector<double>& buffer, std::vector<unsigned>& der_index ) const ;
/// Finish the calculation
  void finish( const std::vector<double>& buffer );
};

inline
void AdjacencyMatrixVessel::setBufferStart( unsigned& start ){
  if( finished ){ bufstart=start; }
  else { Vessel::setBufferStart( start ); } 
}


}
}
#endif

