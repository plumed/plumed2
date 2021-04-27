/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2018 The plumed team
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
#include "core/ActionSetup.h"
#include "core/ActionWithValue.h"
#include "core/ActionRegister.h"
#include "tools/IFile.h"

using namespace std;

namespace PLMD {
namespace setup {

class ReferenceValue :
public ActionSetup,
public ActionWithValue
{
public:
  static void registerKeywords( Keywords& keys );
  explicit ReferenceValue(const ActionOptions&ao);
  void clearDerivatives( const bool& force=false ) {}
  unsigned getNumberOfDerivatives() const { return 0; }
};

PLUMED_REGISTER_ACTION(ReferenceValue,"READ_MATRIX")

void ReferenceValue::registerKeywords( Keywords& keys ) {
  ActionSetup::registerKeywords(keys); ActionWithValue::registerKeywords(keys);
  keys.remove("NUMERICAL_DERIVATIVES"); keys.remove("SERIAL"); keys.remove("TIMINGS");
  keys.add( "hidden", "LABEL", "a label for the action so that its output can be referenced in the input to other actions.  Actions with scalar output are referenced using their label only.  Actions with vector output must have a separate label for every component.  Individual componets are then refered to using label.component" );
  keys.add("compulsory","FILE","an input file containing the matrix of dissimilarities");
}

ReferenceValue::ReferenceValue(const ActionOptions&ao):
  Action(ao),
  ActionSetup(ao),
  ActionWithValue(ao)
{
  std::string fname; parse("FILE",fname); checkRead();
  IFile mfile; mfile.open(fname); std::vector<unsigned> shape(2);
  std::vector<std::vector<double> > dissimilarities;
  // Read in first line
  std::vector<std::string> words; shape[1]=0;
  while( shape[1]==0 ) {
    Tools::getParsedLine( mfile, words );
    shape[1]=words.size();
  }

  std::vector<double> tmpdis( shape[1] );
  for(unsigned j=0; j<shape[1]; ++j) Tools::convert( words[j], tmpdis[j] );
  dissimilarities.push_back( tmpdis );

  while( Tools::getParsedLine( mfile, words ) ) {
    if( words.size()!=shape[1] ) error("bad formatting in matrix file");
    for(unsigned j=0; j<shape[1]; ++j) Tools::convert( words[j], tmpdis[j] );
    dissimilarities.push_back( tmpdis );
  }
  mfile.close(); shape[0] = dissimilarities.size(); 
  addValue( shape ); setNotPeriodic(); getPntrToComponent(0)->alwaysStoreValues();
  for(unsigned i=0;i<shape[0];++i) {
      for(unsigned j=0;j<shape[1];++j) getPntrToComponent(0)->set( i*shape[1] + j, dissimilarities[i][j] );
  }
}

}
}

