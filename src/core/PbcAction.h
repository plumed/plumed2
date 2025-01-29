/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2023 The plumed team
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
#ifndef __PLUMED_core_PbcAction_h
#define __PLUMED_core_PbcAction_h

#include "ActionToPutData.h"
#include "tools/ForwardDecl.h"

#include <vector>
#include <string>

namespace PLMD {

class Pbc;
class DomainDecomposition;

class PbcAction : public ActionToPutData {
  friend class ActionAtomistic;
private:
  DomainDecomposition* interface;
  ForwardDecl<Pbc> pbc_fwd;
  Pbc&   pbc=*pbc_fwd;
  void setPbc();
public:
  explicit PbcAction(const ActionOptions&);
// active methods:
  static void registerKeywords( Keywords& keys );
  Pbc& getPbc();
  void wait() override;
  void readBinary(std::istream&i) override;
  PbcAction* castToPbcAction() noexcept final {
    return this;
  }
};

inline
Pbc& PbcAction::getPbc() {
  return pbc;
}

}
#endif
