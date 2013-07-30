/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

   This file is part of plumed, version 2.0.

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

//+PLUMEDOC GENERIC RANDOM_EXCHANGES
/*
Set random pattern for exchanges

\par Examples

*/
//+ENDPLUMEDOC

class RandomExchanges:
  public Action
{
public:
  static void registerKeywords( Keywords& keys );
  RandomExchanges(const ActionOptions&ao);
  void calculate(){}
  void apply(){}
};

PLUMED_REGISTER_ACTION(RandomExchanges,"RANDOM_EXCHANGES")

void RandomExchanges::registerKeywords( Keywords& keys ){
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

