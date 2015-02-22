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

/**
\ingroup TOOLBOX
Small utility class that contains information about units.

This class can be used to contain in a single place all the
information about units. Units are expressed in terms of
standard PLUMED units, i.e. kj/mol, nm, and ps.
Units can be set as double or as string. In the latter case,
one can also use strings such as kcal/mol.


*/
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
/// Constructor, setting default values (1.0)
  Units();
/// Set energy units from string.
/// Also understands the following strings:
/// kj/mol, kcal/mol, j/mol, and eV.
  void setEnergy(const std::string &);
/// Set time units from string.
/// Also understands the following strings:
/// ps, ns, fs.
  void setTime(const std::string &);
/// Set lengh units from string.
/// Also understands the following strings:
/// nm, A, um.
  void setLength(const std::string &);
/// Set energy units from double.
/// Should be specified in units of kj/mol (e.g. 4.184 means kcal/mol)
  void setEnergy(double);
/// Set time units from double.
/// Should be specified in units of ps (e.g. 0.001 means fs)
  void setTime(double);
/// Set lenght units from double.
/// Should be specified in units of nm (e.g. 0.1 means A)
  void setLength(double);
/// Get energy units as double.
  const double & getEnergy()const;
/// Get length units as double.
  const double & getLength()const;
/// Get time units as double.
  const double & getTime()const;
/// Get energy units as string.
  const std::string & getEnergyString()const;
/// Get length units as string.
  const std::string & getLengthString()const;
/// Get time units as string.
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
