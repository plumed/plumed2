/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2023 The plumed team
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
#include "colvar/Colvar.h"
#include "core/ActionRegister.h"
#include "tools/NeighborList.h"

#include <string>
#include <cmath>

using namespace std;

namespace PLMD{

class NLCounter : public colvar::Colvar {
  bool pbc{true};
  bool firsttime{true};
  bool invalidateList{true};
  bool serial{false};
  std::unique_ptr<PLMD::NeighborList> nl{nullptr};
public:
  static void registerKeywords( Keywords& keys );
  NLCounter(const ActionOptions&);
// active methods:
  virtual void calculate();
  void prepare() override;
};

PLUMED_REGISTER_ACTION(NLCounter,"NLCOUNTER")

void NLCounter::registerKeywords( Keywords& keys ){
  Colvar::registerKeywords( keys );
  keys.add("optional","NL_CUTOFF","The cutoff for the neighbor list");
  keys.add("optional","NL_STRIDE","The frequency with which we are updating the atoms in the neighbor list");
  keys.add("atoms","GROUPA","First list of atoms");
  keys.add("atoms","GROUPB","Second list of atoms (if empty, N*(N-1)/2 pairs in GROUPA are counted)");
}

void NLCounter::prepare() {
  if(nl->getStride()>0) {
    if(firsttime || (getStep()%nl->getStride()==0)) {
      requestAtoms(nl->getFullAtomList());
      invalidateList=true;
      firsttime=false;
    } else {
      requestAtoms(nl->getReducedAtomList());
      invalidateList=false;
      if(getExchangeStep()) error("Neighbor lists should be updated on exchange steps - choose a NL_STRIDE which divides the exchange stride!");
    }
    if(getExchangeStep()) firsttime=true;
  }
}

NLCounter::NLCounter(const ActionOptions&ao):
PLUMED_COLVAR_INIT(ao)
{
  std::vector<AtomNumber> ga_lista,gb_lista;
  parseAtomList("GROUPA",ga_lista);
  parseAtomList("GROUPB",gb_lista);

  bool dopair=false;
  bool nopbc=!pbc;
  parseFlag("NOPBC",nopbc);
  pbc=!nopbc;

  double nl_cut=0.0;
  int nl_st=0;
  
  parse("NL_CUTOFF",nl_cut);
  if(nl_cut<=0.0)
    error("NL_CUTOFF should be explicitly specified and positive");
  parse("NL_STRIDE",nl_st);
  if(nl_st<=0)
    error("NL_STRIDE should be explicitly specified and positive");
  addValueWithDerivatives();
  setNotPeriodic();
  if(gb_lista.size()>0) {
    nl=Tools::make_unique<NeighborList>(ga_lista,gb_lista,serial,dopair,pbc,getPbc(),comm,nl_cut,nl_st);
  } else {
    nl=Tools::make_unique<NeighborList>(ga_lista,serial,pbc,getPbc(),comm,nl_cut,nl_st);
  }

  requestAtoms(nl->getFullAtomList());

  log.printf("  between two groups of %u and %u atoms\n",static_cast<unsigned>(ga_lista.size()),static_cast<unsigned>(gb_lista.size()));
  log.printf("  first group:\n");
  for(unsigned int i=0; i<ga_lista.size(); ++i) {
    if ( (i+1) % 25 == 0 ) log.printf("  \n");
    log.printf("  %d", ga_lista[i].serial());
  }
  log.printf("  \n  second group:\n");
  for(unsigned int i=0; i<gb_lista.size(); ++i) {
    if ( (i+1) % 25 == 0 ) log.printf("  \n");
    log.printf("  %d", gb_lista[i].serial());
  }
  log.printf("  \n");

  if(pbc)
    log.printf("  using periodic boundary conditions\n");
  else
    log.printf("  without periodic boundary conditions\n");
  
  log.printf("  using neighbor lists with\n");
  log.printf("  update every %d steps and cutoff %f\n",nl_st,nl_cut);
  
}

// calculator
void NLCounter::calculate(){
  setValue(nl->size());
}

}



