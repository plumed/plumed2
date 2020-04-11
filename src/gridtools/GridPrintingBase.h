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
#ifndef __PLUMED_gridtools_GridPrintingBase_h
#define __PLUMED_gridtools_GridPrintingBase_h

#include "core/ActionPilot.h"
#include "GridVessel.h"
#include "tools/OFile.h"

namespace PLMD {
namespace gridtools {

class GridPrintingBase : public ActionPilot {
protected:
  GridVessel* ingrid;
  std::string fmt, filename;
  bool output_for_all_replicas;
  std::vector<unsigned> preps;
public:
  static void registerKeywords( Keywords& keys );
  explicit GridPrintingBase(const ActionOptions&ao);
  void calculate() override {}
  void apply() override {}
  void update() override;
  void runFinalJobs() override;
  virtual void printGrid( OFile& ofile ) const=0;
};

}
}
#endif
