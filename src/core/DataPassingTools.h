/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2017-2020 The plumed team
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

namespace PLMD {

class PlumedMain;

class DataPassingTools {
public:
  static std::unique_ptr<DataPassingTools> create(unsigned n);
  virtual int getRealPrecision() const = 0;
  virtual double MD2double(const void*) const=0;
  virtual void double2MD(const double&,void*) const=0;
};

}
#endif
