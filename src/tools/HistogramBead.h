/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2023 The plumed team
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
#ifndef __PLUMED_tools_HistogramBead_h
#define __PLUMED_tools_HistogramBead_h

#include <string>
#include <vector>
#include "Exception.h"
#include "Tools.h"

namespace PLMD {

class Keywords;
class Log;

/**
\ingroup TOOLBOX
A class for calculating whether or not values are within a given range using : \f$ \sum_i \int_a^b G( s_i, \sigma*(b-a) ) \f$
*/

class HistogramBead {
public:
  enum class KernelType {gaussian,triangular};
private:
  bool init{false};
  double lowb{0.0};
  double highb{0.0};
  double width{0.0};
  double cutoff{std::numeric_limits<double>::max()};
  //given how this is initialized, it is not possible to get to an unknown state of the kernel type
  //so, it is not necessary to throw an exception if the type is not set
  KernelType type{KernelType::gaussian};
  enum class Periodicity {unset,periodic,notperiodic};
  Periodicity periodicity{Periodicity::unset};
  double min{0.0};
  double max{0.0};
  double max_minus_min{0.0};
  double inv_max_minus_min{0.0};

  double difference(double d1, double d2 ) const ;
public:
  static void registerKeywords( Keywords& keys );
  static void generateBins( const std::string& params, std::vector<std::string>& bins );
  std::string description() const ;
  // HistogramBead()
  //Non periodic constructor
#pragma acc routine seq
  explicit HistogramBead(KernelType);
  //with period constructor
  HistogramBead(KernelType, double mlow, double mhigh);
  HistogramBead(const HistogramBead&);
  HistogramBead(HistogramBead&&);
  HistogramBead& operator=(const HistogramBead&);
  HistogramBead& operator=(HistogramBead&&);
  // ~HistogramBead();
  bool hasBeenSet() const;
#pragma acc routine seq
  void isNotPeriodic();
  void isPeriodic( double mlow, double mhigh );
  static KernelType getKernelType( const std::string& ktype );
  void setKernelType( const std::string& ktype );
#pragma acc routine seq
  void setKernelType( KernelType ktype );
  void set(const std::string& params, std::string& errormsg);
#pragma acc routine seq
  void set(double l, double h, double w);
#pragma acc routine seq
  double calculate(double x, double&df) const;
  double calculateWithCutoff( double x, double& df ) const;
  double lboundDerivative( const double& x ) const;
  double uboundDerivative( const double& x ) const;
  double getlowb() const ;
  double getbigb() const ;
  double getCutoff() const ;
};

inline
bool HistogramBead::hasBeenSet() const {
  return init;
}

inline
void HistogramBead::isNotPeriodic() {
  periodicity=Periodicity::notperiodic;
}

inline
void HistogramBead::isPeriodic( const double mlow, const double mhigh ) {
  periodicity=Periodicity::periodic;
  min=mlow;
  max=mhigh;
  max_minus_min=max-min;
  plumed_massert(max_minus_min>0, "your function has a very strange domain?");
  inv_max_minus_min=1.0/max_minus_min;
}

inline
double HistogramBead::getlowb() const {
  return lowb;
}

inline
double HistogramBead::getbigb() const {
  return highb;
}

inline
double HistogramBead::getCutoff() const {
  return cutoff*width;
}

}

#endif
