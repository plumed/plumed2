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
#include "core/ParallelTaskManager.h"
#include "tools/RMSD.h"

namespace PLMD {
namespace colvar {

class RMSDVectorData {
public:
  bool squared;
  bool displacement;
  std::string type;
  std::vector<PLMD::RMSD> myrmsd;
  std::vector<double> align, displace, sqrtdisplace;
};

class RMSDVector : public ActionWithVector {
public:
  using input_type = RMSDVectorData;
  using PTM = ParallelTaskManager<RMSDVector>;
private:
  bool firststep;
  bool norm_weights;
  ParallelActionsInput input;
  std::vector<double> force_stash;
  std::vector<double> input_buffer;
  PTM taskmanager;
  static void getPositionsFromInputData( const ParallelActionsInput& input, std::vector<Vector>& pos );
public:
  static void registerKeywords(Keywords& keys);
  explicit RMSDVector(const ActionOptions&);
  unsigned getNumberOfDerivatives() override ;
  int checkTaskIsActive( const unsigned& itask ) const override ;
  void transferStashToValues( const std::vector<unsigned>& partialTaskList, const std::vector<double>& stash ) override ;
  void transferForcesToStash( const std::vector<unsigned>& partialTaskList, std::vector<double>& stash ) const override ;
  static void performTask( std::size_t task_index,
                           const RMSDVectorData& actiondata,
                           ParallelActionsInput& input,
                           ParallelActionsOutput& output );
  static void gatherForces( std::size_t task_index,
                            const RMSDVectorData& actiondata,
                            const ParallelActionsInput& input,
                            View<const double,helpers::dynamic_extent> f,
                            const std::vector<double>& deriv,
                            std::vector<double>& outforces );
  void setReferenceConfigurations();
  void calculate() override ;
  void applyNonZeroRankForces( std::vector<double>& outforces ) override ;
  void apply() override ;
};

}
}
#endif
