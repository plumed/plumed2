/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2020 The plumed team
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
#ifndef __PLUMED_volumes_ActionVolume_h
#define __PLUMED_volumes_ActionVolume_h

#include "tools/HistogramBead.h"
#include "core/ActionWithVector.h"

namespace PLMD {
namespace volumes {

/**
\ingroup INHERIT
This is the abstract base class to use for implementing a new way of defining a particular region of the simulation
box. You can use this to calculate the number of atoms inside that part or the average value of a quantity like the
coordination number inside that part of the cell.
*/

class ActionVolume : public ActionWithVector {
private:
/// The value of sigma
  double sigma;
/// Are we interested in the area outside the colvar
  bool not_in;
/// The kernel type for this histogram
  std::string kerneltype;
protected:
  double getSigma() const ;
  std::string getKernelType() const ;
  Vector getPosition( const unsigned& index ) const ;
  void requestAtoms( const std::vector<AtomNumber> & a );
public:
  static void registerKeywords( Keywords& keys );
  explicit ActionVolume(const ActionOptions&);
  unsigned getNumberOfDerivatives();
  void areAllTasksRequired( std::vector<ActionWithVector*>& task_reducing_actions ) override;
  void getNumberOfTasks( unsigned& ntasks ) override ;
  int checkTaskStatus( const unsigned& taskno, int& flag ) const override;
  void calculate();
  virtual void setupRegions() = 0;
  bool isInSubChain( unsigned& nder ) override ;
  void performTask( const unsigned&, MultiValue& ) const ;
  virtual double calculateNumberInside( const Vector& cpos, Vector& derivatives, Tensor& vir, std::vector<Vector>& refders ) const=0;
};

inline
unsigned ActionVolume::getNumberOfDerivatives() {
  return 3*getNumberOfAtoms()+9;
}

inline
double ActionVolume::getSigma() const {
  return sigma;
}

inline
std::string ActionVolume::getKernelType() const {
  return kerneltype;
}

inline
Vector ActionVolume::getPosition( const unsigned& index ) const {
  if( getConstPntrToComponent(0)->getRank()==0 ) return ActionAtomistic::getPosition( 1 + index );
  return ActionAtomistic::getPosition( getConstPntrToComponent(0)->getShape()[0] + index );
}

}
}
#endif
