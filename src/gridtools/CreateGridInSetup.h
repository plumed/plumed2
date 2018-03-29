/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2018 The plumed team
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
#ifndef __PLUMED_gridtools_ActionWithIntegral_h
#define __PLUMED_gridtools_ActionWithIntegral_h

#include "core/ActionSetup.h"
#include "core/ActionWithValue.h"
#include "GridCoordinatesObject.h"

namespace PLMD {
namespace gridtools {

class CreateGridInSetup :
  public ActionSetup,
  public ActionWithValue
{
protected:
  std::vector<std::string> labels;
  GridCoordinatesObject gridobject;
public:
  static void registerKeywords( Keywords& keys );
  explicit CreateGridInSetup(const ActionOptions&ao);
  void createGridAndValue( const std::string& gtype, const std::vector<bool>& ipbc, const unsigned& nfermi,
                           const std::vector<std::string>& gmin, const std::vector<std::string>& gmax,
                           const std::vector<unsigned>& gbin );
  unsigned getNumberOfDerivatives() const ;
  void clearDerivatives( const bool& force=false ) {}
  void getInfoForGridHeader( std::string& gtype, std::vector<std::string>& argn, std::vector<std::string>& min,
                             std::vector<std::string>& max, std::vector<unsigned>& out_nbin,
                             std::vector<double>& spacing, std::vector<bool>& pbc, const bool& dumpcube ) const ;
  void getGridPointIndicesAndCoordinates( const unsigned& ind, std::vector<unsigned>& indices, std::vector<double>& coords ) const ;
  void getGridPointAsCoordinate( const unsigned& ind, const bool& setlength, std::vector<double>& coords ) const ;
};

}
}
#endif

