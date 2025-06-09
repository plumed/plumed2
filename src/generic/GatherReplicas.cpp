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

There are three ways that you can generate indistinguishable replicas of a system for calculating
ensemble averages:

- You can take a time average.
- You can run a multiple replica simulation and average over the replicas.
- You can do a spatial average and assume that there are multiple replicas of the system in the same simulation box and average over them.

Many actions within PLUMED will perform average over more than one of these type of replica.  However, since v2.10
we have tried to separate out these three ways of generating multiple replicas for averaging.  There are thus many actions
where you take time averages by using [ACCUMULATE](ACCUMULATE.md) or [COLLECT](COLLECT.md).  You take spatial averages by passing
working the vectors of values.  You then use this action to average over replicas.  The way this works in practise is illustrated by the simple
example shown below.

```plumed
#SETTINGS NREPLICAS=2
# Calculate distance between atoms 1 and 2 on for all 2 replicas
d: DISTANCE ATOMS=1,2
# Now gather the values of this distance on all the replicas:
g: GATHER_REPLICAS ARG=d
# Now average the two distances on the two replicas
s: COMBINE ARG=g.rep-1,g.rep-2 COEFFICIENTS=0.5,0.5 PERIODIC=NO
# And print the instaneous average for the distance on the two replicas
PRINT ARG=s FILE=colvar
```

Now suppose that we wanted to calculate a time average of the distance and an average over the replicas.  We could use an input like this:

```plumed
#SETTINGS NREPLICAS=2
d: DISTANCE ATOMS=1,2
g: GATHER_REPLICAS ARG=d
s: COMBINE ARG=g.rep-1,g.rep-2 COEFFICIENTS=0.5,0.5 PERIODIC=NO
# This does the time averaging
a: AVERAGE ARG=s
# And output the final average
PRINT ARG=a FILE=colvar
```

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
  keys.remove("NUMERICAL_DERIVATIVES");
}

GatherReplicas::GatherReplicas( const ActionOptions& ao ):
  Action(ao),
  ActionWithValue(ao),
  ActionWithArguments(ao) {
  if( getNumberOfArguments()!=1 ) {
    error("you can only gather one argument at a time with GatherReplicas");
  }

  std::vector<std::size_t> shape( getPntrToArgument(0)->getShape() );
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
  std::size_t nvals = myarg->getNumberOfValues(), nder = myarg->getNumberOfDerivatives();
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
