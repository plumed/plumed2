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
#include "DataPassingTools.h"
#include "PlumedMain.h"
#include "tools/Tools.h"

namespace PLMD {

double DataPassingTools::getUnitConversion( const std::string& unit ) const {
  if( unit=="energy" ) {
    return MDUnits.getEnergy()/units.getEnergy();
  }
  if( unit=="length" ) {
    return MDUnits.getLength()/units.getLength();
  }
  if( unit=="mass" ) {
    return  MDUnits.getMass()/units.getMass();
  }
  if( unit=="charge" ) {
    return MDUnits.getCharge()/units.getCharge();
  }
  if( unit=="time" ) {
    return MDUnits.getTime()/units.getTime();
  }
  if( unit=="number" ) {
    return 1;
  }
  plumed_error();
}

template <class T>
class DataPassingToolsTyped : public DataPassingTools {
public:
  int getRealPrecision() const override;
  double MD2double(const TypesafePtr & m)const override;
  void double2MD(const double&d,const TypesafePtr & m) const override;
};

std::unique_ptr<DataPassingTools> DataPassingTools::create(unsigned n) {
  if(n==sizeof(double)) {
    return std::make_unique<DataPassingToolsTyped<double>>();
  } else  if(n==sizeof(float)) {
    return std::make_unique<DataPassingToolsTyped<float>>();
  }
  std::string pp;
  Tools::convert(n,pp);
  plumed_merror("cannot create an MD interface with sizeof(real)=="+ pp);
  return NULL;
}

template <class T>
int DataPassingToolsTyped<T>::getRealPrecision() const {
  return sizeof(T);
}

template <class T>
double DataPassingToolsTyped<T>::MD2double(const TypesafePtr & m) const {
  return double(m.template get<T>());
}

template <class T>
void DataPassingToolsTyped<T>::double2MD(const double&d,const TypesafePtr & m) const {
  m.set(T(d));
}

}
