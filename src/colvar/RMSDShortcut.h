/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2024 The plumed team
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
#ifndef __PLUMED_colvar_RMSDShortcut_h
#define __PLUMED_colvar_RMSDShortcut_h

#include "core/ActionShortcut.h"
#include "tools/RMSD.h"

namespace PLMD {
namespace colvar {

class RMSDShortcut : public ActionShortcut {
public:
  static void registerRMSD(Keywords& keys);
  static void registerKeywords(Keywords& keys);
  static void readAlignAndDisplace( ActionWithArguments* action, const bool& norm_weights, std::vector<double>& align, std::vector<double>& displace, std::vector<double>& sqrtdisplace );
  static void setReferenceConfiguration( const unsigned& num, const Value* refarg, const std::vector<double>& align, const std::vector<double>& displace, const std::string& type, const bool& norm_weights, PLMD::RMSD& myrmsd ); 
  explicit RMSDShortcut(const ActionOptions&);
  void createPosVector( const std::string& lab, const PDB& pdb );
  static double calculateDisplacement( const std::string& type, const std::vector<double>& align, const std::vector<double>& displace, const std::vector<double>& sqrtdisplace,
                                       const std::vector<Vector>& pos, PLMD::RMSD& myrmsd, std::vector<Vector>& direction, std::vector<Vector>& der, const bool& squared );
  static void addDisplacementForces( const std::string& type, const std::vector<double>& align, const std::vector<double>& displace, const std::vector<double>& sqrtdisplace,
                                     const std::vector<Vector>& pos, PLMD::RMSD& myrmsd, std::vector<Vector>& direction, std::vector<Vector>& der, Value* myval, const bool& squared );
};

}
}
#endif
