/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2019 The plumed team
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
#include "core/Action.h"
#include "core/ActionRegister.h"
#include "core/PlumedMain.h"
#include "tools/Exception.h"
#include "core/ExchangePatterns.h"

using namespace std;

namespace PLMD {
namespace generic {

//+PLUMEDOC GENERIC RANDOM_EXCHANGES
/*
Set random pattern for exchanges.

In this way, exchanges will not be done between replicas with consecutive index, but
will be done using a random pattern.  Typically used in bias exchange \cite piana.

\par Examples

Using the following three input files one can run a bias exchange
metadynamics simulation using a different angle in each replica.
Exchanges will be randomly tried between replicas 0-1, 0-2 and 1-2

Here is plumed.0.dat
\plumedfile
RANDOM_EXCHANGES
t: TORSION ATOMS=1,2,3,4
METAD ARG=t HEIGHT=0.1 PACE=100 SIGMA=0.3
\endplumedfile

Here is plumed.1.dat
\plumedfile
RANDOM_EXCHANGES
t: TORSION ATOMS=2,3,4,5
METAD ARG=t HEIGHT=0.1 PACE=100 SIGMA=0.3
\endplumedfile

Here is plumed.2.dat
\plumedfile
RANDOM_EXCHANGES
t: TORSION ATOMS=3,4,5,6
METAD ARG=t HEIGHT=0.1 PACE=100 SIGMA=0.3
\endplumedfile

\warning Multi replica simulations are presently only working with gromacs.

\warning The directive should appear in input files for every replicas. In case SEED is specified, it
should be the same in all input files.

*/
//+ENDPLUMEDOC

class RandomExchanges:
  public Action
{
public:
  static void registerKeywords( Keywords& keys );
  explicit RandomExchanges(const ActionOptions&ao);
  void calculate() {}
  void apply() {}
};

PLUMED_REGISTER_ACTION(RandomExchanges,"RANDOM_EXCHANGES")

void RandomExchanges::registerKeywords( Keywords& keys ) {
  Action::registerKeywords(keys);
  keys.add("optional","SEED","seed for random exchanges");
}

RandomExchanges::RandomExchanges(const ActionOptions&ao):
  Action(ao)
{
  plumed.getExchangePatterns().setFlag(ExchangePatterns::RANDOM);
// I convert the seed to -seed because I think it is more general to use a positive seed in input
  int seed=-1;
  parse("SEED",seed);
  if(seed>=0) plumed.getExchangePatterns().setSeed(-seed);
}

}
}

