/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2016 The plumed team
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

namespace PLMD {

class Log;
class Keywords;

/// \ingroup TOOLBOX
/// Small class to compure switching functions.
/// Switching functions are created using set() and
/// then can be used with function calculate() or calculateSqr().
/// Since this is typically computed on a distance vector,
/// the second all (calculateSqr()) allows to skip the calculation
/// of a square root in some case, thus potentially increasing
/// performances.
class SwitchingFunction{
/// This is to check that switching function has been initialized
  bool init;
/// Type of function
  enum {rational,exponential,gaussian,smap,cubic,tanh,matheval,nativeq} type;
/// Inverse of scaling length.
/// We store the inverse to avoid a division
  double invr0;
/// Minimum distance (before this, function is one)
  double d0;
/// Maximum distance (after this, function is zero)
  double dmax;
/// Exponents for rational function
  int nn,mm;
/// Parameters for smap function
  int a,b;
  double c,d;
// nativeq
  double lambda, beta, ref;
/// Square of invr0, useful in calculateSqr()
  double invr0_2;
/// Square of dmax, useful in calculateSqr()
  double dmax_2;
/// Parameters for stretching the function to zero at d_max
  double stretch,shift;
/// Low-level tool to compute rational functions.
/// It is separated since it is called both by calculate() and calculateSqr()
  double do_rational(double rdist,double&dfunc,int nn,int mm)const;
/// Evaluator for matheval:
  void* evaluator;
/// Evaluator for matheval:
  void* evaluator_deriv;
public:
  static void registerKeywords( Keywords& keys );
/// Constructor
  SwitchingFunction();
/// Destructor
  ~SwitchingFunction();
/// Copy constructor
  SwitchingFunction(const SwitchingFunction&);
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

}

#endif

