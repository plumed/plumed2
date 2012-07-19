/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012 The PLUMED team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

   This file is part of PLUMED, version 2.0.

   PLUMED is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   PLUMED is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with PLUMED.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#ifndef __PLUMED_SwitchingFunction_h
#define __PLUMED_SwitchingFunction_h

#include <string>

namespace PLMD {

class Log;

/// \ingroup TOOLBOX
/// Small class to compure switching functions in the form 
/// In the future we might extend it so as to be set using
/// a string:
/// void set(std::string);
/// which can then be parsed for more complex stuff, e.g. exponentials
/// tabulated functions from file, matheval, etc...
class SwitchingFunction{
  bool init;
  enum {spline,exponential,gaussian} type;
  int nn,mm;
  double invr0,d0,dmax;
public:
  static std::string documentation();
  SwitchingFunction();
  void set(int nn,int mm,double r_0,double d_0);
  void set(const std::string& definition, std::string& errormsg);
  std::string description() const ;
  double calculate(double x,double&df)const;
  double get_r0() const;
  void printKeywords( Log& log ) const ;
};

}

#endif

