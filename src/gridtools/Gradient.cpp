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
#include "core/ActionRegister.h"

namespace PLMD {
namespace gridtools {

class Gradient : public ActionWithIntegral {
public:
  static void shortcutKeywords( Keywords& keys );
  static void expandShortcut( const std::string& lab, const std::vector<std::string>& words,
                              const std::map<std::string,std::string>& keys,
                              std::vector<std::vector<std::string> >& actions );
  static void registerKeywords( Keywords& keys );
  explicit Gradient(const ActionOptions&ao);
  void performTask( const unsigned& current, MultiValue& myvals ) const ;
};

PLUMED_REGISTER_ACTION(Gradient,"GRADIENT")
PLUMED_REGISTER_SHORTCUT(Gradient,"GRADIENT")

void Gradient::shortcutKeywords( Keywords& keys ) {
  keys.add("compulsory","ORIGIN","we will use the position of this atom as the origin in our calculation");
  keys.add("compulsory","NBINS","number of bins to use in each direction for the calculation of the gradient");
  keys.add("compulsory","DIR","xyz","the directions in which we are calculating the graident.  Should be x, y, z, xy, xz, yz or xyz");
  keys.add("compulsory","SIGMA","the width of the function to be used for kernel density estimation");
  keys.add("compulsory","KERNEL","gaussian-bin","the type of kernel function to be used in the grids");
  keys.add("compulsory","ATOMS","calculate the gradient of these atoms");
}

void Gradient::expandShortcut( const std::string& lab, const std::vector<std::string>& words,
                               const std::map<std::string,std::string>& keys,
                               std::vector<std::vector<std::string> >& actions ) {
  // First get positions of all atoms relative to origin
  std::vector<std::string> dist_vec; dist_vec.push_back( lab + "_dist:"); dist_vec.push_back("DISTANCE");
  dist_vec.push_back("ORIGIN=" + keys.find("ORIGIN")->second ); dist_vec.push_back("ATOMS=" + keys.find("ATOMS")->second );
  dist_vec.push_back("COMPONENTS"); actions.push_back( dist_vec ); std::string dir = keys.find("DIR")->second;
  // Now constrcut the histograms
  if( dir=="x" || dir=="xy" || dir=="xz" || dir=="xyz" ) {
    std::vector<std::string> xhisto; xhisto.push_back( lab + "_xhisto:"); xhisto.push_back("KDE");
    xhisto.push_back("ARG1=" + lab + "_dist.x"); xhisto.push_back("GRID_BIN=" + keys.find("NBINS")->second );
    xhisto.push_back("KERNEL=" + keys.find("KERNEL")->second ); xhisto.push_back("BANDWIDTH=" + keys.find("SIGMA")->second );
    xhisto.push_back("UNORMALIZED"); actions.push_back( xhisto );
    std::vector<std::string> xgrad; if( dir=="x" ) xgrad.push_back( lab + ":" ); else xgrad.push_back( lab + "_xgrad:" );
    xgrad.push_back("GRADIENT"); xgrad.push_back("ARG=" + lab + "_xhisto"); actions.push_back( xgrad );
  }
  if( dir=="y" || dir=="xy" || dir=="yz" || dir=="xyz" ) {
    std::vector<std::string> xhisto; xhisto.push_back( lab + "_yhisto:"); xhisto.push_back("KDE");
    xhisto.push_back("ARG1=" + lab + "_dist.y"); xhisto.push_back("GRID_BIN=" + keys.find("NBINS")->second );
    xhisto.push_back("KERNEL=" + keys.find("KERNEL")->second ); xhisto.push_back("BANDWIDTH=" + keys.find("SIGMA")->second );
    xhisto.push_back("UNORMALIZED"); actions.push_back( xhisto );
    std::vector<std::string> xgrad; if( dir=="y" ) xgrad.push_back( lab + ":" ); else xgrad.push_back( lab + "_ygrad:" );
    xgrad.push_back("GRADIENT"); xgrad.push_back("ARG=" + lab + "_yhisto"); actions.push_back( xgrad );
  }
  if( dir=="z" || dir=="yz" || dir=="xz" || dir=="xyz" ) {
    std::vector<std::string> xhisto; xhisto.push_back( lab + "_zhisto:"); xhisto.push_back("KDE");
    xhisto.push_back("ARG1=" + lab + "_dist.z"); xhisto.push_back("GRID_BIN=" + keys.find("NBINS")->second );
    xhisto.push_back("KERNEL=" + keys.find("KERNEL")->second ); xhisto.push_back("BANDWIDTH=" + keys.find("SIGMA")->second );
    xhisto.push_back("UNORMALIZED"); actions.push_back( xhisto );
    std::vector<std::string> xgrad; if( dir=="z" ) xgrad.push_back( lab + ":" ); else xgrad.push_back( lab + "_zgrad:" );
    xgrad.push_back("GRADIENT"); xgrad.push_back("ARG=" + lab + "_zhisto"); actions.push_back( xgrad );
  }
  if( dir=="xy" ) {
    std::vector<std::string> combi; combi.push_back( lab + ":" ); combi.push_back("COMBINE");
    combi.push_back("ARG=" + lab + "_xgrad," + lab + "_ygrad"); combi.push_back("PERIODIC=NO");
    actions.push_back( combi );
  } else if( dir=="xz" ) {
    std::vector<std::string> combi; combi.push_back( lab + ":" ); combi.push_back("COMBINE");
    combi.push_back("ARG=" + lab + "_xgrad," + lab + "_zgrad"); combi.push_back("PERIODIC=NO");
    actions.push_back( combi );
  } else if( dir=="yz" ) {
    std::vector<std::string> combi; combi.push_back( lab + ":" ); combi.push_back("COMBINE");
    combi.push_back("ARG=" + lab + "_ygrad," + lab + "_zgrad"); combi.push_back("PERIODIC=NO");
    actions.push_back( combi );
  } else if( dir=="xyz" ) {
    std::vector<std::string> combi; combi.push_back( lab + ":" ); combi.push_back("COMBINE");
    combi.push_back("ARG=" + lab + "_xgrad," + lab + "_ygrad," + lab + "_zgrad"); combi.push_back("PERIODIC=NO");
    actions.push_back( combi );
  }
}

void Gradient::registerKeywords( Keywords& keys ) {
  ActionWithIntegral::registerKeywords( keys );
}

Gradient::Gradient(const ActionOptions&ao):
  Action(ao),
  ActionWithIntegral(ao)
{
  if( gridobject.getDimension()!=1 ) error("gradient should only be used on one dimensional grids");
}

void Gradient::performTask( const unsigned& current, MultiValue& myvals ) const {
  unsigned jval=0; double diff=0;
  if ( gridobject.isPeriodic(0) && current==gridobject.getNbin(false)[0]-1 ) {
    diff = getFunctionValue(current) - getFunctionValue(0); jval = 0;
  } else if( current<gridobject.getNbin(false)[0] ) {
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

}
}
