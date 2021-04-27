/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2018 The VES code team
   (see the PEOPLE-VES file at the root of this folder for a list of names)

   See http://www.ves-code.org for more information.

   This file is part of VES code module.

   The VES code module is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   The VES code module is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with the VES code module.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#ifndef __PLUMED_ves_FermiSwitchingFunction_h
#define __PLUMED_ves_FermiSwitchingFunction_h

#include <string>


namespace PLMD {

class Keywords;

namespace ves {


/// \ingroup TOOLBOX
/// Small class to compute fermi switching function.
/// kept similar to the original SwitchingFunction class.
class FermiSwitchingFunction {
/// This is to check that switching function has been initialized
  bool init;
/// Type of function
  enum {fermi} type;
/// Parameters for fermi function
  double r0_;
  double invr0_;
  double fermi_lambda_;
  double fermi_exp_max_;
  FermiSwitchingFunction& operator=(const FermiSwitchingFunction&);
public:
  static void registerKeywords( Keywords& keys );
/// Constructor
  FermiSwitchingFunction();
/// Destructor
  ~FermiSwitchingFunction();
/// Copy constructor
  FermiSwitchingFunction(const FermiSwitchingFunction&);
  //
  void set(const double, const double, const double fermi_rdist_max=-1.0);
  //
  void set(const std::string& definition, std::string& errormsg);
  //
  std::string description() const ;
  // Compute the switching function.
  // Returns s(x). df will be set to the value of the derivative
  // of the switching function _divided_by_x
  double calculate(double x,double&df)const;
  //
  double get_r0() const;
};

}
}

#endif
