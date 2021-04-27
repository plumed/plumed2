/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2020 The plumed team
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
#ifndef __PLUMED_colvar_RMSD_h
#define __PLUMED_colvar_RMSD_h

#include "Colvar.h"
#include "tools/RMSD.h"

namespace PLMD {
namespace colvar {

class RMSD : public Colvar {
private:
  bool fixed_reference;
  Tensor rot;
  Matrix<std::vector<Vector> > DRotDPos;
  std::vector<Vector> pos, der, direction, centeredpos, centeredreference;
  bool squared;
  bool nopbc;
  bool displacement;
  bool norm_weights;
  std::string type;
  std::vector<double> align,displace,sqrtdisplace;
  PLMD::RMSD myrmsd;
  std::vector<Vector> forcesToApply;
  void setReferenceConfiguration();
public:
  explicit RMSD(const ActionOptions&);
  virtual void calculate() override;
  void apply() override;
  static void registerRMSD(Keywords& keys );
  static void registerKeywords(Keywords& keys);
};

}
}

#endif

