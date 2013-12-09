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
#ifndef __PLUMED_multicolvar_AdjacencyMatrixAction_h
#define __PLUMED_multicolvar_AdjacencyMatrixAction_h

#include "MultiColvarFunction.h"
#include "tools/SwitchingFunction.h"
#include "tools/Matrix.h"
#include "AdjacencyMatrixVessel.h"

namespace PLMD {
namespace multicolvar {

class AdjacencyMatrixAction : public MultiColvarFunction {
friend class AdjacencyMatrixVessel;
private:
/// Flag to make sure derivatives are calculated infrequently
  bool dertime;
/// This is the vessel that stores the adjacency matrix
  AdjacencyMatrixVessel* mat;
///  Tempory vectors for storing vectors
  std::vector<double> tmpdf, orient0, orient1;
/// switching function
  Matrix<SwitchingFunction> switchingFunction;
/// Which matrix elements have value
  bool gathered;
  DynamicList<unsigned> active_elements;
protected:
/// Get the adjacency matrix
  void retrieveMatrix( Matrix<double>& mymatrix );
/// Get number of active matrix elements
  unsigned getNumberOfActiveMatrixElements();
/// Put the indices of the matrix elements in current atoms
  void setMatrixIndexesForTask( const unsigned& ii );
/// Get an element of the adjacency matrix
  double getMatrixElement( const unsigned& ielem ) const ;
/// Add derivatives to a matrix element
  void addDerivativesOnMatrixElement( const unsigned& ielem, const unsigned& jrow, const double& df, Matrix<double>& der );
public:
  static void registerKeywords( Keywords& keys );
  AdjacencyMatrixAction(const ActionOptions&);
  double compute();
  void calculateWeight();
  void doJobsRequiredBeforeTaskList();
/// Finish the calculation
  virtual void completeCalculation()=0;
/// None of these things are allowed
  bool isPeriodic(){ return false; }
  Vector getCentralAtom(){ plumed_merror("cannot find central atoms for adjacency matrix actions"); Vector dum; return dum; }
};

inline
unsigned AdjacencyMatrixAction::getNumberOfActiveMatrixElements(){
  if(!gathered) active_elements.mpi_gatherActiveMembers( comm );
  gathered=true; return active_elements.getNumberActive();
}

inline
double AdjacencyMatrixAction::getMatrixElement( const unsigned& ielem ) const {
  return mat->getComponent( active_elements[ielem], 0 );
}

}
}

#endif
