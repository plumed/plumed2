/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016,2017 The plumed team
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
#include "ActionWithIntegral.h"
#include "core/ActionShortcut.h"
#include "core/ActionRegister.h"

namespace PLMD {
namespace gridtools {

class Gradient : public ActionWithIntegral {
public:
  static void registerKeywords( Keywords& keys );
  explicit Gradient(const ActionOptions&ao);
  void performTask( const unsigned& current, MultiValue& myvals ) const ;
};

PLUMED_REGISTER_ACTION(Gradient,"INTEGRATE_GRADIENT")

void Gradient::registerKeywords( Keywords& keys ) {
  ActionWithIntegral::registerKeywords( keys );
}

Gradient::Gradient(const ActionOptions&ao):
  Action(ao),
  ActionWithIntegral(ao)
{
  if( getGridObject().getDimension()!=1 ) error("gradient should only be used on one dimensional grids");
}

void Gradient::performTask( const unsigned& current, MultiValue& myvals ) const {
  unsigned jval=0; double diff=0;
  if ( getGridObject().isPeriodic(0) && current==getGridObject().getNbin(false)[0]-1 ) {
    diff = getFunctionValue(current) - getFunctionValue(0); jval = 0;
  } else if( current<getGridObject().getNbin(false)[0] ) {
    diff = getFunctionValue(current) - getFunctionValue(current+1); jval = current + 1;
  }
  myvals.setValue( getPntrToOutput(0)->getPositionInStream(), diff*diff );
  if( !doNotCalculateDerivatives() ) {
    myvals.addDerivative( getPntrToOutput(0)->getPositionInStream(), current, 2.0*diff );
    myvals.updateIndex( getPntrToOutput(0)->getPositionInStream(), current );
    myvals.addDerivative( getPntrToOutput(0)->getPositionInStream(), jval, -2.0*diff );
    myvals.updateIndex( getPntrToOutput(0)->getPositionInStream(), jval );
  }
}

class GradientShortcut : public ActionShortcut {
public:
  static void registerKeywords(Keywords& keys);
  explicit GradientShortcut(const ActionOptions&);
};
    
PLUMED_REGISTER_ACTION(GradientShortcut,"GRADIENT")

void GradientShortcut::registerKeywords( Keywords& keys ) {
  ActionShortcut::registerKeywords( keys );
  keys.add("compulsory","ORIGIN","we will use the position of this atom as the origin in our calculation");
  keys.add("compulsory","NBINS","number of bins to use in each direction for the calculation of the gradient");
  keys.add("compulsory","DIR","xyz","the directions in which we are calculating the graident.  Should be x, y, z, xy, xz, yz or xyz");
  keys.add("compulsory","SIGMA","the width of the function to be used for kernel density estimation");
  keys.add("compulsory","KERNEL","gaussian-bin","the type of kernel function to be used in the grids");
  keys.add("compulsory","ATOMS","calculate the gradient of these atoms");
}

GradientShortcut::GradientShortcut(const ActionOptions&ao):
Action(ao),
ActionShortcut(ao)
{
  std::string atom_str; parse("ATOMS",atom_str);
  std::string dir; parse("DIR",dir);
  std::string origin_str; parse("ORIGIN",origin_str);
  std::string nbin_str; parse("NBINS",nbin_str); 
  std::string band_str; parse("SIGMA",band_str);
  std::string kernel_str; parse("KERNEL",kernel_str);
  // First get positions of all atoms relative to origin
  readInputLine( getShortcutLabel() + "_dist: DISTANCES ORIGIN=" + origin_str + " ATOMS=" + atom_str + " COMPONENTS"); 
  // Now constrcut the histograms
  if( dir=="x" || dir=="xy" || dir=="xz" || dir=="xyz" ) {
    readInputLine( getShortcutLabel() + "_xhisto: KDE ARG1=" + getShortcutLabel() + "_dist.x GRID_BIN=" + nbin_str + " KERNEL=" + kernel_str + " BANDWIDTH=" + band_str + " UNORMALIZED");
    std::string thislab = getShortcutLabel() + "_xgrad:"; if( dir=="x" ) thislab = getShortcutLabel() + ":";
    readInputLine( thislab + " INTEGRATE_GRADIENT ARG=" + getShortcutLabel() + "_xhisto"); 
  }
  if( dir=="y" || dir=="xy" || dir=="yz" || dir=="xyz" ) {
    readInputLine( getShortcutLabel() + "_yhisto: KDE ARG1=" + getShortcutLabel() + "_dist.y GRID_BIN=" + nbin_str + " KERNEL=" + kernel_str + " BANDWIDTH=" + band_str + " UNORMALIZED");
    std::string thislab = getShortcutLabel() + "_ygrad:"; if( dir=="y" ) thislab = getShortcutLabel() + ":";
    readInputLine( thislab + " INTEGRATE_GRADIENT ARG=" + getShortcutLabel() + "_yhisto"); 
  }
  if( dir=="z" || dir=="yz" || dir=="xz" || dir=="xyz" ) {
        readInputLine( getShortcutLabel() + "_zhisto: KDE ARG1=" + getShortcutLabel() + "_dist.z GRID_BIN=" + nbin_str + " KERNEL=" + kernel_str + " BANDWIDTH=" + band_str + " UNORMALIZED");
    std::string thislab = getShortcutLabel() + "_zgrad:"; if( dir=="z" ) thislab = getShortcutLabel() + ":";
    readInputLine( thislab + " INTEGRATE_GRADIENT ARG=" + getShortcutLabel() + "_zhisto"); 
  }
  if( dir=="xy" ) {
    readInputLine( getShortcutLabel() + ": COMBINE ARG=" + getShortcutLabel() + "_xgrad," + getShortcutLabel() + "_ygrad PERIODIC=NO");
  } else if( dir=="xz" ) {
    readInputLine( getShortcutLabel() + ": COMBINE ARG=" + getShortcutLabel() + "_xgrad," + getShortcutLabel() + "_zgrad PERIODIC=NO");
  } else if( dir=="yz" ) {
    readInputLine( getShortcutLabel() + ": COMBINE ARG=" + getShortcutLabel() + "_ygrad," + getShortcutLabel() + "_zgrad PERIODIC=NO");
  } else if( dir=="xyz" ) {
    readInputLine( getShortcutLabel() + ": COMBINE ARG=" + getShortcutLabel() + "_xgrad," + getShortcutLabel() + "_ygrad," + getShortcutLabel() + "_zgrad PERIODIC=NO");
  }
}



}
}
