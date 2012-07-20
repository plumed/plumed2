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
#ifndef __PLUMED_Units_h
#define __PLUMED_Units_h

namespace PLMD{

/// Small utility class.
/// It has no implemented methods, and all its data are public.
/// It just simplify the syntax of functions which should pass the
/// value of all the units.
class Units{
public:
/// Units for energy, expressed in kj/mol (e.g. 4.184 means kcal/mol)
  double energy;
/// Units for length, expressed in nm (e.g. 0.1 means A)
  double length;
/// Units for time, expressed in ps (e.g. 0.001 means fs)
  double time;
// Constructor, setting default values (1.0)
  Units();
};

inline
Units::Units():
  energy(1.0),
  length(1.0),
  time(1.0)
{ 
}

}

#endif
