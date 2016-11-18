/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2016 The plumed team
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

namespace PLMD{
namespace generic{

//+PLUMEDOC GENERIC NEIGHBOR2_EXCHANGES
/*
Set exchanges pattern between replicas i and i+2 or i and i-2.

\par Examples

Here is plumed.dat
\verbatim
NEIGHBOR2_EXCHANGES
t: TORSION ATOMS=1,2,3,4
METAD ARG=t HEIGHT=0.1 PACE=100 SIGMA=0.3
\endverbatim

\warning Multi replica simulations are presently only working with gromacs.

\warning The directive should appear in input files for every replicas.

*/
//+ENDPLUMEDOC

class Neighbor2Exchanges:
  public Action
{
public:
  static void registerKeywords( Keywords& keys );
  explicit Neighbor2Exchanges(const ActionOptions&ao);
  void calculate(){}
  void apply(){}
};

PLUMED_REGISTER_ACTION(Neighbor2Exchanges,"NEIGHBOR2_EXCHANGES")

void Neighbor2Exchanges::registerKeywords( Keywords& keys ){
  Action::registerKeywords(keys);
}

Neighbor2Exchanges::Neighbor2Exchanges(const ActionOptions&ao):
Action(ao)
{
  plumed.getExchangePatterns().setFlag(ExchangePatterns::NEIGHBOR2);
}

}
}

