/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2017 The plumed team
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
#ifndef __PLUMED_core_Average_h
#define __PLUMED_core_Average_h
#include "ActionPilot.h"
#include "ActionWithValue.h"
#include "ActionWithArguments.h"

namespace PLMD {

class Average :
public ActionPilot,
public ActionWithValue,
public ActionWithArguments {
private:
  enum {t,f,ndata} normalization;
  bool firststep, clearnextstep;
  unsigned clearstride;
  double lbound, pfactor;
public:
  static void registerKeywords( Keywords& keys );
  explicit Average( const ActionOptions& );
  void clearDerivatives( const bool& force=false ){}
  unsigned getNumberOfDerivatives() const ;
  bool allowComponentsAndValue() const { return true; }
  void getInfoForGridHeader( std::vector<std::string>& argn, std::vector<std::string>& min,
                             std::vector<std::string>& max, std::vector<unsigned>& nbin, 
                             std::vector<double>& spacing, std::vector<bool>& pbc ) const ;
  void getGridPointIndicesAndCoordinates( const unsigned& ind, std::vector<unsigned>& indices, std::vector<double>& coords ) const ;
  void calculate() {}
  void apply() {}
  void update();
};

}
#endif
