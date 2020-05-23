/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2020 The plumed team
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

#include "Between.h"
#include "VesselRegister.h"

namespace PLMD {
namespace vesselbase {

PLUMED_REGISTER_VESSEL(Between,"BETWEEN")

void Between::registerKeywords( Keywords& keys ) {
  FunctionVessel::registerKeywords( keys );
  HistogramBead::registerKeywords( keys );
  keys.addFlag("NORM",false,"calculate the fraction of values rather than the number");
}

void Between::reserveKeyword( Keywords& keys ) {
  keys.reserve("vessel","BETWEEN","calculate the number of values that are within a certain range. "
               "These quantities are calculated using kernel density estimation as described on "
               "\\ref histogrambead.");
  keys.addOutputComponent("between","BETWEEN","the number/fraction of values within a certain range. This is calculated using one of the "
                          "formula described in the description of the keyword so as to make it continuous. "
                          "You can calculate this quantity multiple times using different parameters.");
}

Between::Between( const VesselOptions& da ) :
  FunctionVessel(da)
{
  usetol=true;
  bool isPeriodic=getAction()->isPeriodic();
  double min, max; std::string str_min, str_max;
  if( isPeriodic ) {
    getAction()->retrieveDomain( str_min, str_max );
    Tools::convert(str_min,min); Tools::convert(str_max,max);
  }

  parseFlag("NORM",norm); std::string errormsg;

  hist.set( getAllInput(),errormsg );
  if( !isPeriodic ) hist.isNotPeriodic();
  else hist.isPeriodic( min, max );
  if( errormsg.size()!=0 ) error( errormsg );
}

std::string Between::value_descriptor() {
  if(norm) return "the fraction of values " + hist.description();
  return "the number of values " + hist.description();
}

double Between::calcTransform( const double& val, double& dv ) const {
  double f = hist.calculate(val, dv); return f;
}

double Between::getCutoff() {
  return std::numeric_limits<double>::max();
}

}
}
