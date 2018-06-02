/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2018 The plumed team
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
#include "core/ActionRegister.h"
#include "ActionWithIntegral.h"

//+PLUMEDOC GRIDANALYSIS INTEGRATE_GRID
/*
Calculate the total integral of the function on the input grid

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace gridtools {

class IntegrateGrid : public ActionWithIntegral {
public:
  static void shortcutKeywords( Keywords& keys );
  static void expandShortcut( const std::string& lab, const std::vector<std::string>& words,
                              const std::map<std::string,std::string>& keys,
                              std::vector<std::vector<std::string> >& actions );
  static void registerKeywords( Keywords& keys );
  explicit IntegrateGrid(const ActionOptions&ao);
  void performTask( const unsigned& current, MultiValue& myvals ) const ;
};

PLUMED_REGISTER_ACTION(IntegrateGrid,"INTEGRATE_GRID")
PLUMED_REGISTER_SHORTCUT(IntegrateGrid,"KL_ENTROPY")

void IntegrateGrid::shortcutKeywords( Keywords& keys ) {
  keys.add("compulsory","REFERENCE","a file containing the reference density in grid format");
  keys.add("compulsory","VALUE","the name of the value that should be read from the grid");
}

void IntegrateGrid::expandShortcut( const std::string& lab, const std::vector<std::string>& words,
                                    const std::map<std::string,std::string>& keys,
                                    std::vector<std::vector<std::string> >& actions ) {
  if( words[0]=="KL_ENTROPY" ) {
      // Reference grid object
      std::vector<std::string> ref_keys; ref_keys.push_back( lab + "_ref:"); ref_keys.push_back("REFERENCE_GRID");
      if( !keys.count("VALUE") ) plumed_merror("missing name of VALUE that should be read from input grid");
      ref_keys.push_back("FILE=" + keys.find("REFERENCE")->second ); ref_keys.push_back("VALUE=" + keys.find("VALUE")->second );
      actions.push_back( ref_keys );
      // Compute KL divergence
      std::string input_g=""; 
      for(unsigned i=0;i<words.size();++i) {
          if( words[i].find("ARG=")!=std::string::npos ){ 
              std::size_t eq=words[i].find_first_of("="); input_g=words[i].substr(eq+1); break;
          }
      }
      if( input_g=="") plumed_merror("could not find ARG keyword in input to KL_ENTROPY");
      std::vector<std::string> mt_keys; mt_keys.push_back( lab + "_kl:"); mt_keys.push_back("MATHEVAL");
      mt_keys.push_back("ARG1=" + input_g ); mt_keys.push_back("ARG2=" + lab + "_ref"); mt_keys.push_back("FUNC=y*log(y/(0.5*(x+y)))");
      mt_keys.push_back("PERIODIC=NO"); actions.push_back( mt_keys );
      // Compute integral
      std::vector<std::string> int_keys; int_keys.push_back( lab + ":"); int_keys.push_back("INTEGRATE_GRID");
      int_keys.push_back("ARG=" + lab + "_kl"); actions.push_back( int_keys ); 
  } else {
      plumed_merror("invalid shortcut");
  }
}

void IntegrateGrid::registerKeywords( Keywords& keys ) {
  ActionWithIntegral::registerKeywords( keys );
}

IntegrateGrid::IntegrateGrid(const ActionOptions&ao):
  Action(ao),
  ActionWithIntegral(ao)
{
}

void IntegrateGrid::performTask( const unsigned& current, MultiValue& myvals ) const {
  myvals.setValue( getPntrToOutput(0)->getPositionInStream(), getVolume()*getFunctionValue( current ) );
  if( !doNotCalculateDerivatives() ) {
    myvals.addDerivative( getPntrToOutput(0)->getPositionInStream(), current, getVolume() );
    myvals.updateIndex( getPntrToOutput(0)->getPositionInStream(), current );
  }
}

}
}
