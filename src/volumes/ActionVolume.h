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
#include "core/ActionAtomistic.h"
#include "core/ActionWithValue.h"

namespace PLMD {
namespace volumes {

/**
\ingroup INHERIT
This is the abstract base class to use for implementing a new way of defining a particular region of the simulation
box. You can use this to calculate the number of atoms inside that part or the average value of a quantity like the
coordination number inside that part of the cell.
*/

class ActionVolume :
  public ActionAtomistic,
  public ActionWithValue
{
private:
/// The value of sigma
  double sigma;
/// Are we interested in the area outside the colvar
  bool not_in;
/// The kernel type for this histogram
  std::string kerneltype;
/// The forces that we collect and apply
  std::vector<double> forcesToApply;
protected:
  double getSigma() const ;
  std::string getKernelType() const ;
  Vector getPosition( const unsigned& index ) const ;
  void requestAtoms( const std::vector<AtomNumber> & a );
public:
  static void registerKeywords( Keywords& keys );
  explicit ActionVolume(const ActionOptions&);
  unsigned getNumberOfDerivatives() const ;
  void setupCurrentTaskList();
  void calculate();
  void apply();
  virtual void setupRegions() = 0;
  void performTask( const unsigned&, MultiValue& ) const ;
  virtual double calculateNumberInside( const Vector& cpos, Vector& derivatives, Tensor& vir, std::vector<Vector>& refders ) const=0;
};

inline
unsigned ActionVolume::getNumberOfDerivatives() const {
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
  if( getPntrToOutput(0)->getRank()==0 ) return ActionAtomistic::getPosition( 1 + index ); 
  return ActionAtomistic::getPosition( getPntrToOutput(0)->getShape()[0] + index );
}

}
}
#endif
