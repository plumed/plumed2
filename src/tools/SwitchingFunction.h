/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2023 The plumed team
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
#ifndef __PLUMED_tools_SwitchingFunction_h
#define __PLUMED_tools_SwitchingFunction_h

#include <string>
#include <vector>
#include <memory>
#include "lepton/Lepton.h"

namespace PLMD {

class Log;
class Keywords;
namespace switchContainers {
/// container for the actual switching guntion used by PLMD::SwitchingFunction
class baseSwitch {
protected:
  /// Minimum distance (before this, function is one)
  double d0=0.0;
  /// Maximum distance (after this, function is zero)
  double dmax=0.0;
  /// Square of dmax, useful in calculateSqr()
  double dmax_2=0.0;
  /// We store the inverse to avoid a division
  double invr0=0.0;
  /// Square of invr0, useful in calculateSqr()
  double invr0_2=0.0;
  /// Parameters for stretching the function to zero at d_max
  double stretch=1.0;
  double shift=0.0;
  std::string mytype;
  virtual std::string specificDescription() const;
  //
  virtual double function(double rdist, double& dfunc) const=0;
public:
  baseSwitch(double D0,double DMAX, double R0, std::string name);
  virtual ~baseSwitch();
  ///the driver for the function (prepares rdist or returns 1 or 0 automatically)
  virtual double calculate(double distance, double& dfunc) const;
  virtual double calculateSqr(double distance2, double& dfunc) const;
  void setupStretch();
  void removeStretch();
  std::string description() const;
  double get_d0() const;
  double get_r0() const;
  double get_dmax() const;
  double get_dmax2() const;
};
} // namespace switchContainers

/// \ingroup TOOLBOX
/// Small class to compute switching functions.
/// Switching functions are created using set() and
/// then can be used with function calculate() or calculateSqr().
/// Since this is typically computed on a distance vector,
/// the second all (calculateSqr()) allows to skip the calculation
/// of a square root in some case, thus potentially increasing
/// performances.
class SwitchingFunction {
/// This is to check that switching function has been initialized
  bool init=false;
  std::unique_ptr<switchContainers::baseSwitch> function;
public:
  static void registerKeywords( Keywords& keys );
/// Set a "rational" switching function.
/// Notice that a d_max is set automatically to a value such that
/// f(d_max)=0.00001.
  void set(int nn,int mm,double r_0,double d_0);
/// Set an arbitrary switching function.
/// Parse the string in definition and possibly returns errors
/// in the errormsg string
  void set(const std::string& definition, std::string& errormsg);
/// Returns a string with a description of the switching function
  std::string description() const ;
/// Compute the switching function.
/// Returns s(x). df will be set to the value of the derivative
/// of the switching function _divided_by_x
  double calculate(double x,double&df)const;
/// Compute the switching function.
/// Returns \f$ s(\sqrt{x})\f$ .
/// df will be set to the \f$ \frac{1}{\sqrt{x}}\frac{ds}{d\sqrt{x}}= 2 \frac{ds}{dx}\f$
/// (same as calculate()).
/// The advantage is that in some case the expensive square root can be avoided
/// (namely for rational functions, if nn and mm are even and d0 is zero)
  double calculateSqr(double distance2,double&dfunc)const;
/// Returns d0
  double get_d0() const;
/// Returns r0
  double get_r0() const;
/// Return dmax
  double get_dmax() const;
/// Return dmax squared
  double get_dmax2() const;
};

} //namespace PLMD

#endif
