/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2020 The plumed team
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
#include "core/ActionWithArguments.h"
#include "core/ActionRegister.h"
#include "tools/Random.h"

namespace PLMD {
namespace generic {

//+PLUMEDOC DIMRED CREATE_MASK
/*
Create a mask vector to use for landmark selection

\par Examples

*/
//+ENDPLUMEDOC

class CreateMask :
  public ActionWithValue,
  public ActionWithArguments {
private:
  Random r;
  unsigned nzeros;
  enum {nomask,stride,random} type;
public:
  static void registerKeywords( Keywords& keys );
  CreateMask( const ActionOptions& );
  unsigned getNumberOfDerivatives() override { return 0; }
  void prepare() override ;
  void calculate() override ;
  void apply() override {}
};

PLUMED_REGISTER_ACTION(CreateMask,"CREATE_MASK")

void CreateMask::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys ); ActionWithValue::registerKeywords( keys );
  ActionWithArguments::registerKeywords( keys ); keys.use("ARG");
  keys.add("compulsory","TYPE","the way the zeros are supposed to be set");
  keys.add("compulsory","NZEROS","the number of zeros that you want to put in the mask");
  keys.add("optional","SEED","the seed to use for the random number generator");
  keys.setValueDescription("a vector of zeros and ones that is used that can be used to mask some of the elements in a time series");
}


CreateMask::CreateMask( const ActionOptions& ao ) :
  Action(ao),
  ActionWithValue(ao),
  ActionWithArguments(ao),
  nzeros(0)
{
  if( getNumberOfArguments()!=1 ) error("should only be one argument to this action");
  if( getPntrToArgument(0)->getRank()!=1 ) error("argument should be a vector");
  std::string stype; parse("TYPE",stype); if( stype!="nomask" ) parse("NZEROS",nzeros);

  if( stype=="nomask" ) {
    type=nomask; log.printf("  setting all points in output mask to zero \n");
  } else if( stype=="stride" ) {
    type=stride; log.printf("  setting every %d equally spaced points in output mask to zero \n", nzeros );
  } else if( stype=="random" ) {
    unsigned seed=230623; parse("SEED",seed); r.setSeed(-seed); getPntrToArgument(0)->buildDataStore();
    type=random; log.printf("  choosing %d points to set to non-zero in mask in accordance with input weights \n", nzeros );
  } else error( stype + " is not a valid way input for TYPE");
  std::vector<unsigned> shape(1); shape[0] = getPntrToArgument(0)->getShape()[0];
  addValue( shape ); setNotPeriodic(); getPntrToComponent(0)->buildDataStore();
  for(unsigned i=0; i<shape[0]; ++i) getPntrToComponent(0)->set( i, 1.0 );
}

void CreateMask::prepare() {
  Value* out=getPntrToComponent(0); Value* arg=getPntrToArgument(0);
  if( out->getShape()[0]!=arg->getShape()[0] ) {
    std::vector<unsigned> shape(1); shape[0] = arg->getShape()[0]; out->setShape( shape );
  }
  if( type!=nomask ) {
    for(unsigned i=nzeros; i<out->getShape()[0]; ++i) out->set( i, 1 );
  }
}

void CreateMask::calculate() {
  Value* out=getPntrToComponent(0); Value* arg=getPntrToArgument(0);
  unsigned ns = arg->getShape()[0];
  for(unsigned i=0; i<ns; ++i) out->set( i, 1.0 );

  if( type==stride ) {
    std::size_t ss = int( std::floor( ns / nzeros ) );
    for(unsigned i=0; i<nzeros; ++i) out->set( i*ss, 0.0 );
  } else if( type==random ) {
    for(unsigned i=0; i<nzeros; ++i ) {
      double totweights = 0;
      for(unsigned j=0; j<ns; ++j) {
        if( out->get(j)>0 ) totweights += arg->get(j);
      }
      double rr = r.U01()*totweights; double accum=0;
      for(unsigned j=0; j<ns; ++j) {
        if( out->get(j)>0 ) accum += arg->get(j);
        if( accum<rr ) continue;
        out->set( j, 0 ); break;
      }
    }
  } else if( type==nomask ) {
    for(unsigned i=0; i<ns; ++i) out->set( i, 0.0 );
  } else error("invalid mask creation type");
}

}
}
