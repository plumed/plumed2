/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2023 The plumed team
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
#ifndef __PLUMED_matrixtools_MatrixOperationBase_h
#define __PLUMED_matrixtools_MatrixOperationBase_h

#include "core/ActionWithValue.h"
#include "core/ActionWithArguments.h"
#include "tools/Matrix.h"

namespace PLMD {
namespace matrixtools {

class MatrixOperationBase :
  public ActionWithArguments,
  public ActionWithValue {
private:
/// These are used to hold the matrix
  std::vector<double> MOBvals;
  std::vector<std::pair<unsigned,unsigned> > MOBpairs;
protected:
/// Retrieve a dense version of the ith matrix that is used by this action
  void retrieveFullMatrix( Matrix<double>& mymatrix );
public:
  static void registerKeywords( Keywords& keys );
///
  explicit MatrixOperationBase(const ActionOptions&);
/// Apply the forces
  virtual void apply() override;
/// Get the force on a matrix element
  virtual double getForceOnMatrixElement( const unsigned& jrow, const unsigned& krow ) const = 0 ;
};

}
}
#endif
