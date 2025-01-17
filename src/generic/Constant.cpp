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
#include "core/ActionWithValue.h"
#include "core/ActionRegister.h"
#include "tools/IFile.h"

//+PLUMEDOC COLVAR CONSTANT
/*
Create a constant value that can be passed to actions

Useful in combination with functions that
takes in input constants or parameters.

\par Examples

The following input instructs plumed to compute the distance
between atoms 1 and 2. If this distance is between 1.0 and 2.0, it is
printed. If it is lower than 1.0 (larger than 2.0), 1.0 (2.0) is printed

\plumedfile
cn: CONSTANT VALUES=1.0,2.0
dis: DISTANCE ATOMS=1,2
sss: SORT ARG=cn.v-0,dis,cn.v-1
PRINT ARG=sss.2
\endplumedfile

In case you want to pass a single value you can use VALUE:
\plumedfile
cn: CONSTANT VALUE=1.0
dis: DISTANCE ATOMS=1,2
sss: SORT ARG=cn,dis
PRINT ARG=sss.1
\endplumedfile

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace generic {

class Constant : public ActionWithValue {
public:
  static void registerKeywords( Keywords& keys );
  explicit Constant(const ActionOptions&ao);
  void clearDerivatives( const bool& force=false ) {}
  unsigned getNumberOfDerivatives() override {
    return 0;
  }
  void calculate() override {}
  void apply() override {}
};

PLUMED_REGISTER_ACTION(Constant,"CONSTANT")

void Constant::registerKeywords( Keywords& keys ) {
  Action::registerKeywords(keys);
  ActionWithValue::registerKeywords(keys);
  keys.remove("NUMERICAL_DERIVATIVES");
  keys.add("optional","FILE","an input file containing the matrix");
  keys.add("compulsory","NROWS","0","the number of rows in your input matrix");
  keys.add("compulsory","NCOLS","0","the number of columns in your matrix");
  keys.add("optional","VALUE","the single number that you would like to store");
  keys.add("optional","VALUES","the numbers that are in your constant value");
  keys.addFlag("SCALARS",false,"treat the input list of numbers as a set of scalars");
  keys.addFlag("NOLOG",false,"do not report all the read in scalars in the log");
  keys.addOutputComponent("v","SCALARS","the # value");
  keys.setValueDescription("the constant value that was read from the plumed input");
}

Constant::Constant(const ActionOptions&ao):
  Action(ao),
  ActionWithValue(ao) {
  bool nolog=false;
  parseFlag("NOLOG",nolog);
  bool scalars=false;
  std::string fname, vname;
  parse("FILE",fname);
  std::vector<unsigned> shape;
  std::vector<double> vals;
  if( fname.length()>0 ) {
    IFile mfile;
    mfile.open(fname);
    // Read in first line
    std::vector<std::string> words;
    unsigned nline=0;
    while( nline==0 ) {
      Tools::getParsedLine( mfile, words );
      nline=words.size();
    }
    std::vector<std::vector<double> > dissimilarities;
    if( nline==1 ) {
      shape.resize(1);
      error("invalid matrix in input file");
    }
    shape.resize(2);
    shape[1]=nline;
    std::vector<double> tmpdis( shape[1] );
    for(unsigned j=0; j<shape[1]; ++j) {
      Tools::convert( words[j], tmpdis[j] );
    }
    dissimilarities.push_back( tmpdis );

    while( Tools::getParsedLine( mfile, words ) ) {
      if( words.size()!=nline ) {
        error("bad formatting in matrix file");
      }
      for(unsigned j=0; j<nline; ++j) {
        Tools::convert( words[j], tmpdis[j] );
      }
      dissimilarities.push_back( tmpdis );
    }
    mfile.close();
    shape[0] = dissimilarities.size();
    vals.resize(shape[0]);
    if( shape.size()==2 ) {
      vals.resize( shape[0]*shape[1] );
    }
    for(unsigned i=0; i<shape[0]; ++i) {
      for(unsigned j=0; j<nline; ++j) {
        vals[i*nline+j] = dissimilarities[i][j];
      }
    }
  } else {
    unsigned nr, nc;
    parse("NROWS",nr);
    parse("NCOLS",nc);
    if( nr>0 && nc>0 ) {
      shape.resize(2);
      shape[0]=nr;
      shape[1]=nc;
      vals.resize( nr*nc );
      log.printf("  reading in %d by %d matrix \n", nr, nc );
    } else if( nr>0 || nc>0 ) {
      error("makes no sense to set only one of NROWS and NCOLS to a non-zero value");
    }
    parseVector("VALUES",vals);
    parseFlag("SCALARS",scalars);
    if( vals.size()==0 ) {
      parseVector("VALUE",vals);
      if( vals.size()!=1 ) {
        error("VALUE keyword should take a single scalar");
      }
    } else if( vals.size()==1 ) {
      scalars=false;
    }

    log.printf("  read in %d values :", vals.size() );
    if( !nolog ) {
      for(unsigned i=0; i<vals.size(); ++i) {
        log.printf(" %f", vals[i] );
      }
    }
    log.printf("\n");
    if( !scalars && shape.size()==0 && vals.size()>1 ) {
      shape.resize(1);
      shape[0] = vals.size();
    }
  }
  if( !scalars ) {
    // Now set the value
    addValue( shape );
    setNotPeriodic();
    getPntrToComponent(0)->setConstant();
    for(unsigned i=0; i<vals.size(); ++i) {
      getPntrToComponent(0)->set( i, vals[i] );
    }
  } else {
    for(unsigned i=0; i<vals.size(); i++) {
      std::string num;
      Tools::convert(i,num);
      addComponent("v-"+num);
      componentIsNotPeriodic("v-"+num);
      Value* comp=getPntrToComponent("v-"+num);
      comp->setConstant();
      comp->set(vals[i]);
    }
  }
}

}
}

