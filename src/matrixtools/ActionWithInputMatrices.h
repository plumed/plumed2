/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015-2020 The plumed team
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
#ifndef __PLUMED_matrixtools_ActionWithInputMatrices_h
#define __PLUMED_matrixtools_ActionWithInputMatrices_h

#include "core/ActionWithArguments.h"
#include "core/ActionWithValue.h"
#include "tools/Matrix.h"

namespace PLMD {
namespace matrixtools {

class ActionWithInputMatrices :
  public ActionWithArguments,
  public ActionWithValue
{
protected:
/// These are used to hold the matrix
  std::vector<double> vals;
  std::vector<std::pair<unsigned,unsigned> > pairs;  
/// Add an output value to this action
  void addValue( const std::vector<unsigned>& shape );
/// This sets the values in vals equal to the non-zero elements fo the matrix
/// pairs is then used to keep track of the indices of these non-zero values
  void retrieveEdgeList( const unsigned& imat, unsigned& nedge );
/// Retrieve a dense version of the ith matrix that is used by this action
  void retrieveFullMatrix( const unsigned& imat, Matrix<double>& mymatrix );
public:
  static void registerKeywords( Keywords& keys );
///
  explicit ActionWithInputMatrices(const ActionOptions&);
/// Get the number of derivatives
  unsigned getNumberOfDerivatives() const ;
/// Get the number of columns in the output matrix
  unsigned getNumberOfColumns() const ;
///
  void calculate() override;
///
  void update() override;
///
  void runFinalJobs() override;
///
  virtual void completeMatrixOperations()=0;
};

}
}
#endif
