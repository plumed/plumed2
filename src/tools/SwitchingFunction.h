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

enum class switchType {
  rationalfix12,
  rationalfix10,
  rationalfix8,
  rationalfix6,
  rationalfix4,
  rationalfix2,
  rational,
  rationalFast,
  rationalSimple,
  rationalSimpleFast,
  exponential,
  gaussian,
  fastgaussian,
  smap,
  cubic,
  tanh,
  cosinus,
  nativeq,
  lepton,
  not_initialized
};
struct Data {
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
  //Rational stuff
  int nn=6;
  int mm=12;
  double preRes=0.0;
  double preDfunc=0.0;
  double preSecDev=0.0;
  int nnf=3;
  int mmf=6;
  double preDfuncF=0.0;
  double preSecDevF=0.0;
  //smap stuff
  int a=0;
  int b=0;
  double c=0.0;
  double d=0.0;
  //nativeq
  double beta = 50.0;  // nm-1
  double lambda = 1.8; // unitless
  double ref=0.0;
  static Data init(double D0,double DMAX, double R0);
  void toACCDevice() const;
  void removeFromACCDevice() const;
};

struct Switch {
  virtual double calculate(double distance, double& dfunc) const = 0;
  virtual double calculateSqr(double distance2, double& dfunc) const = 0;
  virtual const Data& getData() const = 0;
  virtual switchType getType() const = 0;
  virtual std::string description() const = 0;
  virtual void setupStretch() = 0;
  virtual ~Switch();
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
  std::unique_ptr<switchContainers::Switch> function{nullptr};
  void copyFunction(const SwitchingFunction&);
public:
  SwitchingFunction();
  SwitchingFunction(const SwitchingFunction&);
  SwitchingFunction(SwitchingFunction&&);
  SwitchingFunction& operator=(const SwitchingFunction&);
  SwitchingFunction& operator=(SwitchingFunction&&);
  ~SwitchingFunction();
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

/// \ingroup TOOLBOX
/// Small class to compute switching functions.
/// Switching functions are created using set() and
/// then can be used with function calculate() or calculateSqr().
/// Since this is typically computed on a distance vector,
/// the second all (calculateSqr()) allows to skip the calculation
/// of a square root in some case, thus potentially increasing
/// performances.
class SwitchingFunctionAccelerable {
  switchContainers::Data switchData;
  switchContainers::switchType type=switchContainers::switchType::not_initialized;
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
#pragma acc routine seq
  double calculate(double x,double&df)const;
/// Compute the switching function.
/// Returns \f$ s(\sqrt{x})\f$ .
/// df will be set to the \f$ \frac{1}{\sqrt{x}}\frac{ds}{d\sqrt{x}}= 2 \frac{ds}{dx}\f$
/// (same as calculate()).
/// The advantage is that in some case the expensive square root can be avoided
/// (namely for rational functions, if nn and mm are even and d0 is zero)
#pragma acc routine seq
  double calculateSqr(double distance2,double&dfunc)const;
/// Returns d0
  double get_d0() const;
/// Returns r0
  double get_r0() const;
/// Return dmax
  double get_dmax() const;
/// Return dmax squared
  double get_dmax2() const;
  void toACCDevice() const;
  void removeFromACCDevice() const;
};


} //namespace PLMD

#endif
