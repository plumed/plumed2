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
  std::vector<double> forcesToApply;
  Value* convertStringToValue( const std::string& name );
  void updateCentralMatrixIndex( const unsigned& ind, MultiValue& myvals ) const ;
protected:
  unsigned ncol_args;
public:
  static void registerKeywords( Keywords& keys );
  explicit VectorProductMatrix(const ActionOptions&);
  bool mustBeTreatedAsDistinctArguments() const ;
  void lockRequests();
  void unlockRequests();
  void calculateNumericalDerivatives( ActionWithValue* a=NULL );
  unsigned getNumberOfDerivatives() const ;
  void buildCurrentTaskList( std::vector<unsigned>& tflags );
  void calculate();
  void performTask( const unsigned& task_index, MultiValue& myvals ) const ;
  void performTask( const std::string& controller, const unsigned& index1, const unsigned& index2, MultiValue& myvals ) const ;
  virtual double computeVectorProduct( const unsigned& index1, const unsigned& index2, 
                                       const std::vector<double>& vec1, const std::vector<double>& vec2, 
                                       std::vector<double>& dvec1, std::vector<double>& dvec2, MultiValue& myvals ) const = 0;
  void apply();
};

inline
unsigned VectorProductMatrix::getNumberOfDerivatives() const  {
  unsigned nat_der = 0;
  if( getNumberOfAtoms()>0 ) nat_der = 3*getNumberOfAtoms()+9;
  if( ncol_args>0 ) return nat_der + (getPntrToArgument(0)->getShape()[0]+getPntrToArgument(ncol_args)->getShape()[0])*getNumberOfArguments()/2;
  return nat_der + getPntrToArgument(0)->getShape()[0]*getNumberOfArguments();
}

inline
bool VectorProductMatrix::mustBeTreatedAsDistinctArguments() const {
  return true;
}


}
}
#endif

