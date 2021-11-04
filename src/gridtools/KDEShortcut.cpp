/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2017 The plumed team
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
#include "KDEShortcut.h"
#include "core/ActionRegister.h"

namespace PLMD {
namespace gridtools {
    
PLUMED_REGISTER_ACTION(KDEShortcut,"KDE")
    
void KDEShortcut::registerKeywords( Keywords& keys ) {
  HistogramBase::registerKeywords( keys );
  keys.add("compulsory","KERNEL","GAUSSIAN","the kernel function you are using.  More details on  the kernels available "
           "in plumed plumed can be found in \\ref kernelfunctions.");
  keys.add("compulsory","BANDWIDTH","the bandwidths for kernel density esimtation");
}

KDEShortcut::KDEShortcut(const ActionOptions& ao):
Action(ao),
ActionShortcut(ao)
{
  std::string height, height_str, kernel=""; 
  parse("HEIGHTS",height); std::string ktype; parse("KERNEL",ktype);
  if( ktype!="DISCRETE" ) { 
      std::vector<std::string> bwidths; parseVector("BANDWIDTH",bwidths); double bnum;
      if( bwidths.size()>1 || Tools::convert( bwidths[0], bnum ) ) {
          KDEShortcut::convertBandwiths( getShortcutLabel(), bwidths, this ); kernel = " KERNEL=" + ktype + " METRIC=" + getShortcutLabel() + "_icov";
          // Compute the normalizing constant
          std::string pstr; Tools::convert( sqrt(pow(2*pi,bwidths.size())), pstr );
          if( ktype=="gaussian" || ktype=="GAUSSIAN" ) {
              readInputLine( getShortcutLabel() + "_bwprod: PRODUCT ARG=" + getShortcutLabel() + "_cov");
              readInputLine( getShortcutLabel() + "_vol: CUSTOM ARG=" + getShortcutLabel() + "_bwprod FUNC=(sqrt(x)*" + pstr + ") PERIODIC=NO"); 
              if( height.length()>0 ) readInputLine( getShortcutLabel() + "_height: CUSTOM ARG1=" + height + " ARG2=" + getShortcutLabel() + "_vol FUNC=x/y PERIODIC=NO");
              else readInputLine( getShortcutLabel() + "_height: CUSTOM ARG1=" + getShortcutLabel() + "_vol FUNC=1/x PERIODIC=NO");
              height_str = " HEIGHTS=" + getShortcutLabel() + "_height";
          } else if( height.length()>0 ) {
              height_str = " HEIGHTS=" + height;
          } else if( ktype.find("bin")!=std::string::npos ) {
              readInputLine( getShortcutLabel() + "_height: CONSTANT VALUE=1.0"); height_str = " HEIGHTS=" + getShortcutLabel() + "_height";
          } else error("you need to set the heights of the kernel functions you are using from the covariance so they are normalised");
      } else if( height.length()>0 ) height_str = " HEIGHTS=" + height;
  } else { kernel = " KERNEL=DISCRETE"; if( height.length()>0 ) height_str = " HEIGHTS=" + height; }
  HistogramBase::createKDEObject( getShortcutLabel(), "KDE", height, height_str + kernel, this );
}

void KDEShortcut::convertBandwiths( const std::string& lab, const std::vector<std::string>& bwidths, ActionShortcut* action ) {
  double bw; std::string band="VALUES=" + bwidths[0];
  for(unsigned i=0;i<bwidths.size();++i) {
      if( !Tools::convert( bwidths[i], bw ) ) action->error("could not convert input bandwidth to real number");
      if( i>0 ) band += "," + bwidths[i]; 
  }
  action->readInputLine( lab + "_sigma: CONSTANT_VALUE " + band );
  action->readInputLine( lab + "_cov: CUSTOM ARG=" + lab + "_sigma FUNC=x*x PERIODIC=NO" );
  action->readInputLine( lab + "_icov: CUSTOM ARG=" + lab + "_cov FUNC=1/x PERIODIC=NO" );
}

}
}
