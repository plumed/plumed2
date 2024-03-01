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
#ifndef __PLUMED_colvar_RMSDVector_h
#define __PLUMED_colvar_RMSDVector_h

#include "core/ActionWithVector.h"
#include "tools/RMSD.h"

namespace PLMD {
namespace colvar {

class RMSDVector : public ActionWithVector {

  bool firststep;
  bool squared;
  bool displacement;
  bool norm_weights;
  std::string type;
  std::vector<PLMD::RMSD> myrmsd;
  std::vector<double> align, displace, sqrtdisplace;
  double calculateRMSD( const unsigned& current, std::vector<Vector>& pos, std::vector<Vector>& der, std::vector<Vector>& direction ) const ;
public:
  static void registerKeywords(Keywords& keys);
  explicit RMSDVector(const ActionOptions&);
  unsigned getNumberOfDerivatives() override ;
  void performTask( const unsigned& current, MultiValue& myvals ) const override ;
  void gatherStoredValue( const unsigned& valindex, const unsigned& code, const MultiValue& myvals,
                          const unsigned& bufstart, std::vector<double>& buffer ) const override ;
  void gatherForcesOnStoredValue( const Value* myval, const unsigned& itask, const MultiValue& myvals, std::vector<double>& forces ) const override ;
  void setReferenceConfigurations();
  void calculate() override ;
  bool checkForTaskForce( const unsigned& itask, const Value* myval ) const override ;
  void apply() override ;
};

}
}
#endif
