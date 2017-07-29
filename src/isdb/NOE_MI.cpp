/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2014-2017 The plumed team
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
#include "Meta.h"
#include "colvar/ActionRegister.h"
#include "tools/NeighborList.h"
#include "tools/OpenMP.h"
#include "tools/Pbc.h"

#include <string>
#include <cmath>

using namespace std;

namespace PLMD {
namespace isdb {

//+PLUMEDOC COLVAR NOEMI
/*
Calculates NOEMI intensities as sums of 1/r^6, also averaging over multiple equivalent atoms
or ambiguous NOEMI.

Each NOEMI is defined by two groups containing the same number of atoms, distances are
calculated in pairs, transformed in 1/r^6, summed and saved as components.

\f[
NOEMI() = (\frac{1}{N_{eq}}\sum_j^{N_{eq}} (\frac{1}{r_j^6}))
\f]

Intensities can then in principle ensemble averaged using \ref ENSEMBLE and used to
calculate a scoring function for example with \ref METAINFERENCE.

\par Examples
In the following examples three noes are defined, the first is calculated based on the distances
of atom 1-2 and 3-2; the second is defined by the distance 5-7 and the third by the distances
4-15,4-16,8-15,8-16.

\plumedfile
NOEMI ...
GROUPA1=1,3 GROUPB1=2,2
GROUPA2=5 GROUPB2=7
GROUPA3=4,4,8,8 GROUPB3=15,16,15,16
LABEL=noes
... NOEMI

PRINT ARG=noes.* FILE=colvar
\endplumedfile

*/
//+ENDPLUMEDOC

class NOEMI :
  public Meta
{
private:
  bool             pbc;
  vector<unsigned> nga;
  NeighborList     *nl;
public:
  static void registerKeywords( Keywords& keys );
  explicit NOEMI(const ActionOptions&);
  ~NOEMI();
  void calculate_simple();
  void calculate();
};

PLUMED_REGISTER_ACTION(NOEMI,"NOEMI")

void NOEMI::registerKeywords( Keywords& keys ) {
  componentsAreNotOptional(keys);
  useCustomisableComponents(keys);
  Meta::registerKeywords(keys);
  keys.addFlag("NOPBC",false,"ignore the periodic boundary conditions when calculating distances");
  keys.add("numbered","GROUPA","the atoms involved in each of the contacts you wish to calculate. "
           "Keywords like GROUPA1, GROUPA2, GROUPA3,... should be listed and one contact will be "
           "calculated for each ATOM keyword you specify.");
  keys.add("numbered","GROUPB","the atoms involved in each of the contacts you wish to calculate. "
           "Keywords like GROUPB1, GROUPB2, GROUPB3,... should be listed and one contact will be "
           "calculated for each ATOM keyword you specify.");
  keys.reset_style("GROUPA","atoms");
  keys.reset_style("GROUPB","atoms");
  keys.addFlag("ADDEXP",false,"Set to TRUE if you want to have fixed components with the experimental reference values.");
  keys.add("numbered","NOEMIDIST","Add an experimental value for each NOEMI.");
  keys.addOutputComponent("noe","default","the # NOEMI");
  keys.addOutputComponent("exp","ADDEXP","the # NOEMI experimental distance");
}

NOEMI::NOEMI(const ActionOptions&ao):
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
  nl= new NeighborList(ga_lista,gb_lista,true,pbc,getPbc());

  bool addexp=false;
  parseFlag("ADDEXP",addexp);

  vector<double> noedist;
  if(addexp) {
    noedist.resize( nga.size() );
    unsigned ntarget=0;

    for(unsigned i=0; i<nga.size(); ++i) {
      if( !parseNumbered( "NOEMIDIST", i+1, noedist[i] ) ) break;
      ntarget++;
    }
    if( ntarget!=nga.size() ) error("found wrong number of NOEMIDIST values");
  }

  // Ouput details of all contacts
  unsigned index=0;
  for(unsigned i=0; i<nga.size(); ++i) {
    log.printf("  The %uth NOEMI is calculated using %u equivalent couples of atoms\n", i, nga[i]);
    for(unsigned j=0; j<nga[i]; j++) {
      log.printf("    couple %u is %d %d.\n", j, ga_lista[index].serial(), gb_lista[index].serial() );
      index++;
    }
  }

  if(pbc)      log.printf("  using periodic boundary conditions\n");
  else         log.printf("  without periodic boundary conditions\n");

  if(!doscore_) {
    for(unsigned i=0; i<nga.size(); i++) {
      string num; Tools::convert(i,num);
      addComponentWithDerivatives("noe_"+num);
      componentIsNotPeriodic("noe_"+num);
    }
    if(addexp) {
      for(unsigned i=0; i<nga.size(); i++) {
        string num; Tools::convert(i,num);
        addComponent("exp_"+num);
        componentIsNotPeriodic("exp_"+num);
        Value* comp=getPntrToComponent("exp_"+num);
        comp->set(noedist[i]);
      }
    }
  } else {
    setParameters(noedist);
    for(unsigned i=0; i<nga.size(); i++) {
      string num; Tools::convert(i,num);
      addComponent("noe_"+num);
      componentIsNotPeriodic("noe_"+num);
    }
    for(unsigned i=0; i<nga.size(); i++) {
      string num; Tools::convert(i,num);
      addComponent("exp_"+num);
      componentIsNotPeriodic("exp_"+num);
      Value* comp=getPntrToComponent("exp_"+num);
      comp->set(noedist[i]);
    }
  }

  requestAtoms(nl->getFullAtomList());
  setDerivatives();
  checkRead();
}

NOEMI::~NOEMI() {
  delete nl;
}

void NOEMI::calculate_simple()
{
  const unsigned ngasz=nga.size();

  #pragma omp parallel for num_threads(OpenMP::getNumThreads())
  for(unsigned i=0; i<ngasz; i++) {
    Tensor dervir;
    double noe=0;
    unsigned index=0;
    for(unsigned k=0; k<i; k++) index+=nga[k];
    const double c_aver=1./static_cast<double>(nga[i]);
    string num; Tools::convert(i,num);
    Value* val=getPntrToComponent("noe_"+num);
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
      const double tmpir6=c_aver*ir6;
      const double tmpir8=c_aver*ir8;

      noe += tmpir6;
      Vector deriv = tmpir8*distance;
      dervir += Tensor(distance,deriv);
      setAtomsDerivatives(val, i0,  deriv);
      setAtomsDerivatives(val, i1, -deriv);
    }
    val->set(noe);
    setBoxDerivatives(val, dervir);
  }
}


void NOEMI::calculate()
{
  if(!doscore_) calculate_simple();
  else {
    vector<Vector> deriv(getNumberOfAtoms());
    const unsigned ngasz=nga.size();

    #pragma omp parallel for num_threads(OpenMP::getNumThreads())
    for(unsigned i=0; i<ngasz; i++) {
      double noe=0;
      unsigned index=0;
      for(unsigned k=0; k<i; k++) index+=nga[k];
      const double c_aver=1./static_cast<double>(nga[i]);
      string num; Tools::convert(i,num);
      Value* val=getPntrToComponent("noe_"+num);
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
        const double tmpir6=c_aver*ir6;
        const double tmpir8=c_aver*ir8;

        noe += tmpir6;
        Vector tmpder = tmpir8*distance;
        deriv[i0] += tmpder;
        deriv[i1] -= tmpder;
      }
      val->set(noe);
      setCalcData(i, noe);
    }
    /* Metainference */
    /* 1) collect weights */
    double fact = 0.;
    double var_fact = 0.;
    get_weights(fact, var_fact);

    /* 2) calculate average */
    vector<double> mean(getNarg(),0);
    // this is the derivative of the mean with respect to the argument
    vector<double> dmean_x(getNarg(),fact);
    // this is the derivative of the mean with respect to the bias
    vector<double> dmean_b(getNarg(),0);
    // calculate it
    replica_averaging(fact, mean, dmean_b);

    /* 3) calculates parameters */
    get_sigma_mean(fact, var_fact, mean);

    /* 4) run monte carlo */
    doMonteCarlo(mean);

    /* 5) calculate score */
    double score = getScore(mean, dmean_x, dmean_b);
    setScore(score);

    /* calculate final derivatives */
    Tensor dervir;
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

        dervir += Tensor(distance,deriv[i0]*getMetaDer(i));
        setAtomsDerivatives(val, i0, deriv[i0]*getMetaDer(i));
        setAtomsDerivatives(val, i1, deriv[i1]*getMetaDer(i));
      }
    }
    setBoxDerivatives(val, dervir);
  }
}

}
}
