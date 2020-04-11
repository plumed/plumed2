/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2014-2020 The plumed team
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
#include "MetainferenceBase.h"
#include "core/ActionRegister.h"
#include "tools/NeighborList.h"
#include "tools/Pbc.h"

#include <string>
#include <cmath>
#include <memory>

using namespace std;

namespace PLMD {
namespace isdb {

//+PLUMEDOC ISDB_COLVAR NOE
/*
Calculates NOE intensities as sums of 1/r^6, also averaging over multiple equivalent atoms
 or ambiguous NOE.

Each NOE is defined by two groups containing the same number of atoms, distances are
calculated in pairs, transformed in 1/r^6, summed and saved as components.

\f[
NOE() = (\frac{1}{N_{eq}}\sum_j^{N_{eq}} (\frac{1}{r_j^6}))
\f]

NOE can be used to calculate a Metainference score over one or more replicas using the intrinsic implementation
of \ref METAINFERENCE that is activated by DOSCORE.

\par Examples
In the following examples three noes are defined, the first is calculated based on the distances
of atom 1-2 and 3-2; the second is defined by the distance 5-7 and the third by the distances
4-15,4-16,8-15,8-16. \ref METAINFERENCE is activated using DOSCORE.

\plumedfile
NOE ...
GROUPA1=1,3 GROUPB1=2,2
GROUPA2=5 GROUPB2=7
GROUPA3=4,4,8,8 GROUPB3=15,16,15,16
DOSCORE
LABEL=noes
... NOE

PRINT ARG=noes.* FILE=colvar
\endplumedfile

*/
//+ENDPLUMEDOC

class NOE :
  public MetainferenceBase
{
private:
  bool             pbc;
  vector<unsigned> nga;
  std::unique_ptr<NeighborList> nl;
  unsigned         tot_size;
public:
  static void registerKeywords( Keywords& keys );
  explicit NOE(const ActionOptions&);
  void calculate() override;
  void update() override;
};

PLUMED_REGISTER_ACTION(NOE,"NOE")

void NOE::registerKeywords( Keywords& keys ) {
  componentsAreNotOptional(keys);
  MetainferenceBase::registerKeywords(keys);
  keys.addFlag("NOPBC",false,"ignore the periodic boundary conditions when calculating distances");
  keys.add("numbered","GROUPA","the atoms involved in each of the contacts you wish to calculate. "
           "Keywords like GROUPA1, GROUPA2, GROUPA3,... should be listed and one contact will be "
           "calculated for each ATOM keyword you specify.");
  keys.add("numbered","GROUPB","the atoms involved in each of the contacts you wish to calculate. "
           "Keywords like GROUPB1, GROUPB2, GROUPB3,... should be listed and one contact will be "
           "calculated for each ATOM keyword you specify.");
  keys.reset_style("GROUPA","atoms");
  keys.reset_style("GROUPB","atoms");
  keys.add("numbered","NOEDIST","Add an experimental value for each NOE.");
  keys.addOutputComponent("noe","default","the # NOE");
  keys.addOutputComponent("exp","NOEDIST","the # NOE experimental distance");
}

NOE::NOE(const ActionOptions&ao):
  PLUMED_METAINF_INIT(ao),
  pbc(true)
{
  bool nopbc=!pbc;
  parseFlag("NOPBC",nopbc);
  pbc=!nopbc;

  // Read in the atoms
  vector<AtomNumber> t, ga_lista, gb_lista;
  for(int i=1;; ++i ) {
    parseAtomList("GROUPA", i, t );
    if( t.empty() ) break;
    for(unsigned j=0; j<t.size(); j++) ga_lista.push_back(t[j]);
    nga.push_back(t.size());
    t.resize(0);
  }
  vector<unsigned> ngb;
  for(int i=1;; ++i ) {
    parseAtomList("GROUPB", i, t );
    if( t.empty() ) break;
    for(unsigned j=0; j<t.size(); j++) gb_lista.push_back(t[j]);
    ngb.push_back(t.size());
    if(ngb[i-1]!=nga[i-1]) error("The same number of atoms is expected for the same GROUPA-GROUPB couple");
    t.resize(0);
  }
  if(nga.size()!=ngb.size()) error("There should be the same number of GROUPA and GROUPB keywords");
  // Create neighbour lists
  nl.reset( new NeighborList(ga_lista,gb_lista,true,pbc,getPbc()) );

  // Optionally add an experimental value (like with RDCs)
  vector<double> noedist;
  noedist.resize( nga.size() );
  unsigned ntarget=0;
  for(unsigned i=0; i<nga.size(); ++i) {
    if( !parseNumbered( "NOEDIST", i+1, noedist[i] ) ) break;
    ntarget++;
  }
  bool addexp=false;
  if(ntarget!=nga.size() && ntarget!=0) error("found wrong number of NOEDIST values");
  if(ntarget==nga.size()) addexp=true;
  if(getDoScore()&&!addexp) error("with DOSCORE you need to set the NOEDIST values");

  // Ouput details of all contacts
  unsigned index=0;
  for(unsigned i=0; i<nga.size(); ++i) {
    log.printf("  The %uth NOE is calculated using %u equivalent couples of atoms\n", i, nga[i]);
    for(unsigned j=0; j<nga[i]; j++) {
      log.printf("    couple %u is %d %d.\n", j, ga_lista[index].serial(), gb_lista[index].serial() );
      index++;
    }
  }
  tot_size = index;

  if(pbc)      log.printf("  using periodic boundary conditions\n");
  else         log.printf("  without periodic boundary conditions\n");

  log << " Bibliography" << plumed.cite("Bonomi, Camilloni, Bioinformatics, 33, 3999 (2017)") << "\n";

  if(!getDoScore()) {
    for(unsigned i=0; i<nga.size(); i++) {
      string num; Tools::convert(i,num);
      addComponentWithDerivatives("noe-"+num);
      componentIsNotPeriodic("noe-"+num);
    }
    if(addexp) {
      for(unsigned i=0; i<nga.size(); i++) {
        string num; Tools::convert(i,num);
        addComponent("exp-"+num);
        componentIsNotPeriodic("exp-"+num);
        Value* comp=getPntrToComponent("exp-"+num);
        comp->set(noedist[i]);
      }
    }
  } else {
    for(unsigned i=0; i<nga.size(); i++) {
      string num; Tools::convert(i,num);
      addComponent("noe-"+num);
      componentIsNotPeriodic("noe-"+num);
    }
    for(unsigned i=0; i<nga.size(); i++) {
      string num; Tools::convert(i,num);
      addComponent("exp-"+num);
      componentIsNotPeriodic("exp-"+num);
      Value* comp=getPntrToComponent("exp-"+num);
      comp->set(noedist[i]);
    }
  }

  requestAtoms(nl->getFullAtomList(), false);
  if(getDoScore()) {
    setParameters(noedist);
    Initialise(nga.size());
  }
  setDerivatives();
  checkRead();
}

void NOE::calculate()
{
  const unsigned ngasz=nga.size();
  vector<Vector> deriv(tot_size, Vector{0,0,0});

  #pragma omp parallel for num_threads(OpenMP::getNumThreads())
  for(unsigned i=0; i<ngasz; i++) {
    Tensor dervir;
    double noe=0;
    unsigned index=0;
    for(unsigned k=0; k<i; k++) index+=nga[k];
    string num; Tools::convert(i,num);
    Value* val=getPntrToComponent("noe-"+num);
    // cycle over equivalent atoms
    for(unsigned j=0; j<nga[i]; j++) {
      const unsigned i0=nl->getClosePair(index+j).first;
      const unsigned i1=nl->getClosePair(index+j).second;

      Vector distance;
      if(pbc) distance=pbcDistance(getPosition(i0),getPosition(i1));
      else    distance=delta(getPosition(i0),getPosition(i1));

      const double ir2=1./distance.modulo2();
      const double ir6=ir2*ir2*ir2;
      const double ir8=6*ir6*ir2;

      noe += ir6;
      deriv[index+j] = ir8*distance;
      if(!getDoScore()) {
        dervir += Tensor(distance, deriv[index+j]);
        setAtomsDerivatives(val, i0,  deriv[index+j]);
        setAtomsDerivatives(val, i1, -deriv[index+j]);
      }
    }
    val->set(noe);
    if(!getDoScore()) {
      setBoxDerivatives(val, dervir);
    } else setCalcData(i, noe);
  }

  if(getDoScore()) {
    /* Metainference */
    Tensor dervir;
    double score = getScore();
    setScore(score);

    /* calculate final derivatives */
    Value* val=getPntrToComponent("score");
    for(unsigned i=0; i<ngasz; i++) {
      unsigned index=0;
      for(unsigned k=0; k<i; k++) index+=nga[k];
      // cycle over equivalent atoms
      for(unsigned j=0; j<nga[i]; j++) {
        const unsigned i0=nl->getClosePair(index+j).first;
        const unsigned i1=nl->getClosePair(index+j).second;

        Vector distance;
        if(pbc) distance=pbcDistance(getPosition(i0),getPosition(i1));
        else    distance=delta(getPosition(i0),getPosition(i1));

        dervir += Tensor(distance,deriv[index+j]*getMetaDer(i));
        setAtomsDerivatives(val, i0,  deriv[index+j]*getMetaDer(i));
        setAtomsDerivatives(val, i1, -deriv[index+j]*getMetaDer(i));
      }
    }
    setBoxDerivatives(val, dervir);
  }
}

void NOE::update() {
  // write status file
  if(getWstride()>0&& (getStep()%getWstride()==0 || getCPT()) ) writeStatus();
}

}
}
