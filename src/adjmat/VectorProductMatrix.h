/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015-2017 The plumed team
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
#ifndef __PLUMED_adjmat_VectorProductMatrixBase_h

#include "core/ActionAtomistic.h"
#include "core/ActionWithArguments.h"
#include "core/ActionWithValue.h"

namespace PLMD {
namespace adjmat {

class VectorProductMatrix :
  public ActionAtomistic,
  public ActionWithArguments,
  public ActionWithValue
{
private:
  void updateCentralMatrixIndex( const unsigned& ind, MultiValue& myvals ) const ;
protected:
  unsigned ncol_args;
  unsigned nspAa,nspAb;
  unsigned narg_derivatives;
  std::vector<double> forcesToApply;
public:
  static void registerKeywords( Keywords& keys );
  explicit VectorProductMatrix(const ActionOptions&);
  unsigned getNumberOfDerivatives() const ;
  bool canBeAfterInChain( ActionWithValue* av ) const ;
  bool mustBeTreatedAsDistinctArguments() const ;
  void lockRequests();
  void unlockRequests();
  void calculateNumericalDerivatives( ActionWithValue* a=NULL );
  void calculate();
  void performTask( const unsigned& task_index, MultiValue& myvals ) const ;
  bool performTask( const std::string& controller, const unsigned& index1, const unsigned& index2, MultiValue& myvals ) const ;
  virtual double computeVectorProduct( const unsigned& index1, const unsigned& index2,
                                       const std::vector<double>& vec1, const std::vector<double>& vec2,
                                       std::vector<double>& dvec1, std::vector<double>& dvec2, MultiValue& myvals ) const = 0;
  void apply();
};

inline
bool VectorProductMatrix::mustBeTreatedAsDistinctArguments() const {
  return true;
}


}
}
#endif

