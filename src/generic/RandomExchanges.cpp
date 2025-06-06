/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2023 The plumed team
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

namespace PLMD {
namespace generic {

//+PLUMEDOC GENERIC RANDOM_EXCHANGES
/*
Set random pattern for exchanges.

This command is typically used if you are using the bias exchange method that is
discussed in the paper in the bibliography.  When it is present it tells PLUMED
not do do exchanges between replicas with consecutive indices and instad to use
use a random pattern.

The following three example input files show how one can run a bias exchange
metadynamics simulation using a different angle in each replica.
Exchanges are randomly tried between replicas 0-1, 0-2 and 1-2

Here is plumed.0.dat

```plumed
RANDOM_EXCHANGES
t: TORSION ATOMS=1,2,3,4
METAD ARG=t HEIGHT=0.1 PACE=100 SIGMA=0.3
```

Here is plumed.1.dat

```plumed
RANDOM_EXCHANGES
t: TORSION ATOMS=2,3,4,5
METAD ARG=t HEIGHT=0.1 PACE=100 SIGMA=0.3
```

Here is plumed.2.dat

```plumed
RANDOM_EXCHANGES
t: TORSION ATOMS=3,4,5,6
METAD ARG=t HEIGHT=0.1 PACE=100 SIGMA=0.3
```

Notice that you can perform the same calculation with the following single PLUMED input file:

```plumed
#SETTINGS NREPLICAS=3

RANDOM_EXCHANGES SEED=23
t1: TORSION ATOMS=1,2,3,4
t2: TORSION ATOMS=2,3,4,5
t3: TORSION ATOMS=3,4,5,6
METAD ARG=@replicas:t1,t2,t3 HEIGHT=0.1 PACE=100 SIGMA=0.3
```

!!! caution ""

    Multi replica simulations are presently only working with gromacs.

!!! caution ""

    The directive should appear in the input file for every replica. If SEED is specified, it
    should be the same in all input files.

*/
//+ENDPLUMEDOC

class RandomExchanges:
  public Action {
public:
  static void registerKeywords( Keywords& keys );
  explicit RandomExchanges(const ActionOptions&ao);
  void calculate() override {}
  void apply() override {}
};

PLUMED_REGISTER_ACTION(RandomExchanges,"RANDOM_EXCHANGES")

void RandomExchanges::registerKeywords( Keywords& keys ) {
  Action::registerKeywords(keys);
  keys.add("optional","SEED","seed for random exchanges");
  keys.addDOI("10.1021/jp067873l");
}

RandomExchanges::RandomExchanges(const ActionOptions&ao):
  Action(ao) {
  plumed.getExchangePatterns().setFlag(ExchangePatterns::RANDOM);
// I convert the seed to -seed because I think it is more general to use a positive seed in input
  int seed=-1;
  parse("SEED",seed);
  if(seed>=0) {
    plumed.getExchangePatterns().setSeed(-seed);
  }
}

}
}

