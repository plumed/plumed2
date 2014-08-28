/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2014 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

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
#ifndef __PLUMED_tools_Units_h
#define __PLUMED_tools_Units_h

#include <string>

namespace PLMD{

/// Small utility class.
/// It has no implemented methods, and all its data are public.
/// It just simplify the syntax of functions which should pass the
/// value of all the units.
class Units{
/// Units for energy, expressed in kj/mol (e.g. 4.184 means kcal/mol)
  double energy;
  std::string energyString;
/// Units for length, expressed in nm (e.g. 0.1 means A)
  double length;
  std::string lengthString;
/// Units for time, expressed in ps (e.g. 0.001 means fs)
  double time;
  std::string timeString;
public:
// Constructor, setting default values (1.0)
  Units();
  void setEnergy(const std::string &);
  void setTime(const std::string &);
  void setLength(const std::string &);
  void setEnergy(const double);
  void setTime(const double);
  void setLength(const double);
  const double & getEnergy()const;
  const double & getLength()const;
  const double & getTime()const;
  const std::string & getEnergyString()const;
  const std::string & getLengthString()const;
  const std::string & getTimeString()const;
};

inline
const double & Units::getEnergy()const{
  return energy;
}

inline
const double & Units::getLength()const{
  return length;
}

inline
const double & Units::getTime()const{
  return time;
}

inline
const std::string & Units::getEnergyString()const{
  return energyString;
}

inline
const std::string & Units::getLengthString()const{
  return lengthString;
}

inline
const std::string & Units::getTimeString()const{
  return timeString;
}


}

#endif
