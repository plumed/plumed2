/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2019 The plumed team
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
private:
  bool init;
  double lowb;
  double highb;
  double width;
  double cutoff;
  enum {gaussian,triangular} type;
  enum {unset,periodic,notperiodic} periodicity;
  double min, max, max_minus_min, inv_max_minus_min;
  double difference( const double& d1, const double& d2 ) const ;
public:
  static void registerKeywords( Keywords& keys );
  static void generateBins( const std::string& params, std::vector<std::string>& bins );
  HistogramBead();
  std::string description() const ;
  bool hasBeenSet() const;
  void isNotPeriodic();
  void isPeriodic( const double& mlow, const double& mhigh );
  void setKernelType( const std::string& ktype );
  void set(const std::string& params, std::string& errormsg);
  void set(double l, double h, double w);
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
  periodicity=notperiodic;
}

inline
void HistogramBead::isPeriodic( const double& mlow, const double& mhigh ) {
  periodicity=periodic; min=mlow; max=mhigh;
  max_minus_min=max-min;
  plumed_massert(max_minus_min>0, "your function has a very strange domain?");
  inv_max_minus_min=1.0/max_minus_min;
}

inline
double HistogramBead::getlowb() const { return lowb; }

inline
double HistogramBead::getbigb() const { return highb; }

inline
double HistogramBead::getCutoff() const { return cutoff*width; }

inline
double HistogramBead::difference( const double& d1, const double& d2 ) const {
  if(periodicity==notperiodic) {
    return d2-d1;
  } else if(periodicity==periodic) {
    // Make sure the point is in the target range
    double newx=d1*inv_max_minus_min;
    newx=Tools::pbc(newx);
    newx*=max_minus_min;
    return d2-newx;
  } else plumed_merror("periodicty was not set");
  return 0;
}

}

#endif
