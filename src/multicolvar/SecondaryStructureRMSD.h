/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

   This file is part of plumed, version 2.0.

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
#ifndef __PLUMED_multicolvar_SecondaryStructureRMSD_h
#define __PLUMED_multicolvar_SecondaryStructureRMSD_h

#include "MultiColvar.h"
#include <vector>

namespace PLMD {

class RMSD;
class DRMSD;

namespace multicolvar {

/// Base action for calculating things like AlphRMSD, AntibetaRMSD, etc

class SecondaryStructureRMSD : public MultiColvar {
private:
  std::string alignType;
  std::vector<Vector> deriv, new_deriv;
  std::vector<RMSD*> secondary_rmsd;
  std::vector<DRMSD*> secondary_drmsd;
protected:
  void setSecondaryStructure( std::vector<Vector>& structure, double bondlength, double units );
  bool usingRMSD() const ;
public:
  static void registerKeywords( Keywords& keys );
  SecondaryStructureRMSD(const ActionOptions&);
  virtual ~SecondaryStructureRMSD();
  virtual double compute( const unsigned& j, const std::vector<Vector>& pos );
  unsigned getNumberOfFieldDerivatives();
  bool isPeriodic(){ return false; }
};

}
}

#endif
