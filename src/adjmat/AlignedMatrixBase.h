/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2019 The plumed team
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
#ifndef __PLUMED_adjmat_AlignedMatrixBase_h
#define __PLUMED_adjmat_AlignedMatrixBase_h

#include "tools/SwitchingFunction.h"
#include "AdjacencyMatrixBase.h"


namespace PLMD {
namespace adjmat {

class AlignedMatrixBase : public AdjacencyMatrixBase {
private:
  unsigned ncol_t;
/// switching function
  Matrix<SwitchingFunction> switchingFunction;
public:
/// Create manual
  static void registerKeywords( Keywords& keys );
/// Constructor
  explicit AlignedMatrixBase(const ActionOptions&);
///
  void setupConnector( const unsigned& id, const unsigned& i, const unsigned& j, const std::vector<std::string>& desc );
  virtual void readOrientationConnector( const unsigned& i, const unsigned& j, const std::vector<std::string>& desc ) = 0;
/// This actually calculates the value of the contact function
  double calculateWeight( const unsigned& taskCode, const double& weight, multicolvar::AtomValuePack& myatoms ) const ;
  double compute( const unsigned& tindex, multicolvar::AtomValuePack& myatoms ) const ;
/// This transforms the dot product
  virtual double computeVectorFunction( const unsigned& iv, const unsigned& jv,
                                        const Vector& conn, const std::vector<double>& vec1, const std::vector<double>& vec2,
                                        Vector& dconn, std::vector<double>& dvec1, std::vector<double>& dvec2 ) const = 0;
};


}
}
#endif
