/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2017-2023 The plumed team
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
#ifndef __PLUMED_core_DataPassingTools_h
#define __PLUMED_core_DataPassingTools_h

#include <string>
#include <memory>
#include <map>
#include "tools/TypesafePtr.h"
#include "tools/Units.h"

namespace PLMD {

class PlumedMain;

class DataPassingTools {
  friend class PlumedMain;
private:
/// The units used in the MD code and PLUMED
  Units units;
  Units MDUnits;
/// Is the code using natural units
  bool usingNaturalUnits;
public:
/// Virtual destructor, just to allow inheritance.
  virtual ~DataPassingTools() {}
  static std::unique_ptr<DataPassingTools> create(unsigned n);
  virtual int getRealPrecision() const = 0;
  double getUnitConversion( const std::string& unit ) const ;
  virtual double MD2double(const TypesafePtr & m) const=0;
  virtual void double2MD(const double&,const TypesafePtr & m) const=0;
};

}
#endif
