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
#ifndef __PLUMED_mapping_TrigonometricPathVessel_h
#define __PLUMED_mapping_TrigonometricPathVessel_h

#include "vesselbase/StoreDataVessel.h"
#include "reference/ReferenceValuePack.h"
#include "reference/Direction.h"
#include "core/Value.h"
#include "Mapping.h"

namespace PLMD {
namespace mapping {

class TrigonometricPathVessel : public vesselbase::StoreDataVessel {
  friend class AdaptivePath;
private:
  Value* sp;
  Value* zp;
  Mapping* mymap;
  double dx;
  std::vector<double> cargs;
  unsigned iclose1, iclose2;
  Direction projdir;
  std::vector<double> mypack1_stashd_args;
  std::vector<Vector> mypack1_stashd_atoms;
  MultiValue mydpack1, mydpack2, mydpack3;
  ReferenceValuePack mypack1, mypack2, mypack3;
public:
  static void registerKeywords( Keywords& keys );
  static void reserveKeyword( Keywords& keys );
  explicit TrigonometricPathVessel( const vesselbase::VesselOptions& da );
  std::string description() override;
  void resize() override;
  void finish( const std::vector<double>& buffer ) override;
  bool applyForce(std::vector<double>&) override;
};

}
}
#endif
