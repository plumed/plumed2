/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2020 The plumed team
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
#ifndef __PLUMED_dimred_PCA_h
#define __PLUMED_dimred_PCA_h

#include "DimensionalityReductionBase.h"
#include "reference/ReferenceConfiguration.h"
#include "reference/Direction.h"
#include "tools/PDB.h"

namespace PLMD {
namespace dimred {

class PCA : public DimensionalityReductionBase {
  friend class OutputPCAProjection;
private:
/// The way we are measuring distances
  std::string mtype;
/// The position of the reference configuration (the one we align to)
  PDB mypdb;
/// The eigenvectors for the displacements in argument space
  std::string ofilename, fmt;
/// The eigenvectors that we are using
  std::unique_ptr<ReferenceConfiguration> myref;
  std::vector<Direction> directions;
public:
  static void registerKeywords( Keywords& keys );
  explicit PCA(const ActionOptions&ao);
  void performAnalysis() override;
  void getProjection( const unsigned& idata, std::vector<double>& point, double& weight ) override;
  void getProjection( analysis::DataCollectionObject& myidata, std::vector<double>& point );
  void calculateProjections( const Matrix<double>&, Matrix<double>& ) override { plumed_error(); }
  void setTargetDistance( const unsigned&, const double& ) override { plumed_error(); }
  double calculateStress( const std::vector<double>& pp, std::vector<double>& der ) override { plumed_error(); }
};

}
}
#endif
