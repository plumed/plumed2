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
#include "core/ActionWithValue.h"
#include "core/ActionWithArguments.h"
#include "core/ActionRegister.h"
#include "tools/Communicator.h"

//+PLUMEDOC ANALYSIS GATHER_REPLICAS
/*
Create a vector that contains the copies of the input quantities from all replicas

\par Examples


*/
//+ENDPLUMEDOC

namespace PLMD {
namespace generic {

class GatherReplicas :
  public ActionWithValue,
  public ActionWithArguments {
private:
  unsigned nreplicas;
public:
  static void registerKeywords( Keywords& keys );
  explicit GatherReplicas( const ActionOptions& );
  unsigned getNumberOfDerivatives();
  void calculate();
  void apply();
};

PLUMED_REGISTER_ACTION(GatherReplicas,"GATHER_REPLICAS")

void GatherReplicas::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys );
  ActionWithValue::registerKeywords( keys );
  ActionWithArguments::registerKeywords( keys );
  keys.addInputKeyword("compulsory","ARG","scalar/vector/matrix/grid","the argument from the various replicas that you would like to gather");
  keys.addOutputComponent("rep","default","scalar/vector/matrix/grid","the input arguments for each of the replicas");
}

GatherReplicas::GatherReplicas( const ActionOptions& ao ):
  Action(ao),
  ActionWithValue(ao),
  ActionWithArguments(ao) {
  if( getNumberOfArguments()!=1 ) {
    error("you can only gather one argument at a time with GatherReplicas");
  }

  std::vector<unsigned> shape( getPntrToArgument(0)->getShape() );
  std::string min, max;
  nreplicas=multi_sim_comm.Get_size();
  bool periodic=false;
  if( getPntrToArgument(0)->isPeriodic() ) {
    periodic=true;
    getPntrToArgument(0)->getDomain( min, max );
  }

  for(unsigned i=0; i<nreplicas; ++i) {
    std::string num;
    Tools::convert( i+1, num);
    if( getPntrToArgument(0)->hasDerivatives() ) {
      addComponentWithDerivatives( "rep-" + num, shape );
    } else {
      addComponent( "rep-" + num, shape );
    }
    if( periodic ) {
      componentIsPeriodic( "rep-" + num, min, max );
    } else {
      componentIsNotPeriodic( "rep-" + num );
    }
  }
}

unsigned GatherReplicas::getNumberOfDerivatives() {
  return getPntrToArgument(0)->getNumberOfDerivatives();
}

void GatherReplicas::calculate() {
  Value* myarg = getPntrToArgument(0);
  unsigned nvals = myarg->getNumberOfValues(), nder = myarg->getNumberOfDerivatives();
  std::vector<double> dval( nvals*(1+nder) ), datap(nreplicas*nvals*(1+nder) );
  for(unsigned i=0; i<nvals; ++i) {
    dval[i*(1+nder)] = myarg->get(i);
    if( myarg->getRank()==0 ) {
      for(unsigned j=0; j<nder; ++j) {
        dval[i*(1+nder)+1+j] = myarg->getDerivative(j);
      }
    } else if( myarg->hasDerivatives() ) {
      for(unsigned j=0; j<nder; ++j) {
        dval[i*(1+nder)+1+j] = myarg->getGridDerivative( i, j );
      }
    }
  }
  if(comm.Get_rank()==0) {
    multi_sim_comm.Allgather(dval,datap);
  }

  for(unsigned k=0; k<nreplicas; k++) {
    Value* myout = getPntrToComponent(k);
    if( myout->getNumberOfDerivatives()!=myarg->getNumberOfDerivatives() ) {
      myout->resizeDerivatives( myarg->getNumberOfDerivatives() );
    }
    unsigned sstart=k*nvals*(1+nder);
    for(unsigned i=0; i<nvals; ++i) {
      myout->set( i, datap[sstart+i*(1+nder)] );
      if( myarg->getRank()==0 ) {
        for(unsigned j=0; j<nder; ++j) {
          myout->setDerivative( j, dval[i*(1+nder)+1+j] );
        }
      } else if( myarg->hasDerivatives() ) {
        for(unsigned j=0; j<nder; ++j) {
          myout->addGridDerivatives( i, j, dval[i*(1+nder)+1+j] );
        }
      }
    }
  }
}

void GatherReplicas::apply() {
  if( doNotCalculateDerivatives() ) {
    return;
  }
  error("apply has not been implemented for GatherReplicas");
}

}
}
