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

class ConstantValue :
public ActionSetup,
public ActionWithValue
{
public:
  static void registerKeywords( Keywords& keys );
  explicit ConstantValue(const ActionOptions&ao);
  void clearDerivatives( const bool& force=false ) {}
  unsigned getNumberOfDerivatives() const { return 0; }
};

PLUMED_REGISTER_ACTION(ConstantValue,"CONSTANT_VALUE")

void ConstantValue::registerKeywords( Keywords& keys ) {
  ActionSetup::registerKeywords(keys); ActionWithValue::registerKeywords(keys);
  keys.remove("NUMERICAL_DERIVATIVES"); keys.remove("SERIAL"); keys.remove("TIMINGS");
  keys.add( "hidden", "LABEL", "a label for the action so that its output can be referenced in the input to other actions.  Actions with scalar output are referenced using their label only.  Actions with vector output must have a separate label for every component.  Individual componets are then refered to using label.component" );
  keys.add("optional","FILE","an input file containing the matrix");
  keys.add("compulsory","NROWS","0","the number of rows in your input matrix");
  keys.add("compulsory","NCOLS","0","the number of columns in your matrix");
  keys.add("optional","VALUES","the elements of your matrix");
}

ConstantValue::ConstantValue(const ActionOptions&ao):
  Action(ao),
  ActionSetup(ao),
  ActionWithValue(ao)
{
  std::string fname; parse("FILE",fname); 
  std::vector<unsigned> shape; std::vector<double> vals;
  if( fname.length()>0 ) {
       IFile mfile; mfile.open(fname);
       // Read in first line
       std::vector<std::string> words; unsigned nline=0;
       while( nline==0 ) {
         Tools::getParsedLine( mfile, words );
         nline=words.size();
       }
       if( nline==1 ) shape.resize(1); 
       else { shape.resize(2); shape[1]=nline; }

       std::vector<double> tmpdis( shape[1] );
       std::vector<std::vector<double> > dissimilarities;
       for(unsigned j=0; j<shape[1]; ++j) Tools::convert( words[j], tmpdis[j] );
       dissimilarities.push_back( tmpdis );

       while( Tools::getParsedLine( mfile, words ) ) {
         if( words.size()!=nline ) error("bad formatting in matrix file");
         for(unsigned j=0; j<nline; ++j) Tools::convert( words[j], tmpdis[j] );
         dissimilarities.push_back( tmpdis );
       }
       mfile.close(); shape[0] = dissimilarities.size(); vals.resize(shape[0]);
       if( shape.size()==2 ) vals.resize( shape[0]*shape[1] );
       for(unsigned i=0;i<shape[0];++i) {
           for(unsigned j=0;j<nline;++j) vals[i*nline+j] = dissimilarities[i][j];
       }
  } else {
       unsigned nr, nc; parse("NROWS",nr); parse("NCOLS",nc); 
       if( nr>0 && nc>0 ) { 
           shape.resize(2); shape[0]=nr; shape[1]=nc; vals.resize( nr*nc );
           log.printf("  reading in %d by %d matrix \n", nr, nc ); 
       } else if( nr>0 || nc>0 ) error("makes no sense to set only one of NROWS and NCOLS to a non-zero value");
       parseVector("VALUES",vals); 
       if( getLabel().find("_ones")==std::string::npos ) {
           log.printf("  read in %d values :", vals.size() );
           for(unsigned i=0; i<vals.size(); ++i) log.printf(" %f", vals[i] );
       } else log.printf("  vector of %d ones", vals.size() );
       log.printf("\n"); if( shape.size()==0 && vals.size()>1 ) { shape.resize(1); shape[0] = vals.size(); }
  }
  // Now set the value
  addValue( shape ); setNotPeriodic(); 
  getPntrToComponent(0)->setConstant(); getPntrToComponent(0)->alwaysStoreValues();
  for(unsigned i=0; i<vals.size(); ++i) getPntrToComponent(0)->set( i, vals[i] );
}

}
}

