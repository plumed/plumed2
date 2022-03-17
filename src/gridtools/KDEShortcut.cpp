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
#include "Histogram.h"
#include "KDE.h"
#include "SphericalKDE.h"
#include "core/ActionRegister.h"

namespace PLMD {
namespace gridtools {
 
template<class T>
class KDEShortcut : public ActionShortcut {
public:
  static void registerKeywords(Keywords& keys);
  explicit KDEShortcut(const ActionOptions&);
};
   
typedef KDEShortcut<KDE> kde;
PLUMED_REGISTER_ACTION(kde,"KDE")
typedef KDEShortcut<SphericalKDE> skde;
PLUMED_REGISTER_ACTION(skde,"SPHERICAL_KDE")
 
template<class T>   
void KDEShortcut<T>::registerKeywords( Keywords& keys ) {
  T::registerKeywords( keys );
  if( keys.exists("KERNEL") ) keys.add("compulsory","BANDWIDTH","the bandwidths for kernel density esimtation");
}

template<class T>
KDEShortcut<T>::KDEShortcut(const ActionOptions& ao):
Action(ao),
ActionShortcut(ao)
{
  std::string height, height_str, kernel=""; parse("HEIGHTS",height);
  if( getName()=="KDE") {
      std::string ktype; parse("KERNEL",ktype);
      if( ktype!="DISCRETE" ) { 
          std::vector<std::string> bwidths; parseVector("BANDWIDTH",bwidths); double bnum;
          if( bwidths.size()>1 || Tools::convert( bwidths[0], bnum ) ) {
              HistogramTools::convertBandwiths( getShortcutLabel(), bwidths, this ); kernel = " KERNEL=" + ktype + " METRIC=" + getShortcutLabel() + "_icov";
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
  } else if( getName()=="SPHERICAL_KDE" ) {
      if( height.length()>0 ) height_str = " HEIGHTS=" + height;
  }
  std::string inp; bool uflag; parseFlag("UNORMALIZED",uflag);
  // Deal with the weights if we are doing averages on a grid
  if( height.length()>0 && !uflag ) {
    inp = getShortcutLabel() + "_unorm: " + getName() + "_CALC " + convertInputLineToString();
    readInputLine( getShortcutLabel() + "_hsum: SUM ARG=" + height + " PERIODIC=NO");
    inp = inp + " UNORMALIZED";
  } else if( !uflag ) {
     inp = getShortcutLabel() + ": " + getName() + "_CALC " + convertInputLineToString();
  } else {
     inp = getShortcutLabel() + ": " + getName() + "_CALC UNORMALIZED " + convertInputLineToString();
  }
  readInputLine( inp + height_str + kernel );
  if( height.length()>0 && !uflag ) {
    readInputLine(  getShortcutLabel() + ": MATHEVAL ARG1=" + getShortcutLabel() + "_unorm ARG2=" + getShortcutLabel() + "_hsum FUNC=x/y PERIODIC=NO");
  }
}

}
}
