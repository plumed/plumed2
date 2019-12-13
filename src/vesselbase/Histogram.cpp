/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2019 The plumed team
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
#include "tools/HistogramBead.h"
#include "VesselRegister.h"
#include "ShortcutVessel.h"

namespace PLMD {
namespace vesselbase {

class Histogram : public ShortcutVessel {
public:
  static void registerKeywords( Keywords& keys );
  static void reserveKeyword( Keywords& keys );
  explicit Histogram( const VesselOptions& da );
};

PLUMED_REGISTER_VESSEL(Histogram,"HISTOGRAM")

void Histogram::registerKeywords( Keywords& keys ) {
  ShortcutVessel::registerKeywords( keys );
  HistogramBead::registerKeywords( keys );
  keys.add("compulsory","NBINS","The number of equal width bins you want to divide the range into");
  keys.addFlag("NORM",false,"calculate the fraction of values rather than the number");
  keys.add("compulsory","COMPONENT","1","the component of the vector for which to calculate this quantity");
}

void Histogram::reserveKeyword( Keywords& keys ) {
  keys.reserve("vessel","HISTOGRAM","calculate how many of the values fall in each of the bins of a histogram. "
               "This shortcut allows you to calculates NBIN quantities like BETWEEN.");
}

Histogram::Histogram( const VesselOptions& da ):
  ShortcutVessel(da)
{
  bool norm; parseFlag("NORM",norm); std::string normstr="";
  if(norm) normstr=" NORM";
  std::string compstr; parse("COMPONENT",compstr);
  normstr+=" COMPONENT=" + compstr;
  std::vector<std::string> bins; HistogramBead::generateBins( getAllInput(), bins );
  for(unsigned i=0; i<bins.size(); ++i) addVessel("BETWEEN",bins[i] + normstr);
}

}
}
