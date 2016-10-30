/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2014-2016 The plumed team
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
#include "Colvar.h"
#include "ActionRegister.h"
#include "tools/NeighborList.h"
#include "tools/OpenMP.h"

#include <string>
#include <cmath>

using namespace std;

namespace PLMD {
namespace colvar {

//+PLUMEDOC COLVAR NOE 
/*
Calculates NOE intensities as sums of 1/r^6, also averaging over multiple equivalent atoms
or ambiguous NOE.

Each NOE is defined by two groups containing the same number of atoms, distances are
calculated in pairs, transformed in 1/r^6, summed and saved as components.

\f[
NOE() = (\frac{1}{N_{eq}}\sum_j^{N_{eq}} (\frac{1}{r_j^6}))^{\frac{-1}{6}} 
\f]

Intensities can then in principle ensemble averaged using \ref ENSEMBLE and used to
calculate a scoring function for example with \ref METAINFERENCE.

\par Examples
In the following examples three noes are defined, the first is calculated based on the distances
of atom 1-2 and 3-2; the second is defined by the distance 5-7 and the third by the distances
4-15,4-16,8-15,8-16.

\verbatim
NOE ...
GROUPA1=1,3 GROUPB1=2,2 
GROUPA2=5 GROUPB2=7
GROUPA3=4,4,8,8 GROUPB3=15,16,15,16
LABEL=noes
... NOE

PRINT ARG=noes.* FILE=colvar
\endverbatim
(See also \ref PRINT) 

*/
//+ENDPLUMEDOC

class NOE : public Colvar {   
private:
  bool             pbc;
  vector<unsigned> nga;
  NeighborList     *nl;
public:
  static void registerKeywords( Keywords& keys );
  explicit NOE(const ActionOptions&);
  ~NOE();
  virtual void calculate();
};

PLUMED_REGISTER_ACTION(NOE,"NOE")

void NOE::registerKeywords( Keywords& keys ){
  Colvar::registerKeywords( keys );
  componentsAreNotOptional(keys);
  useCustomisableComponents(keys);
  keys.add("numbered","GROUPA","the atoms involved in each of the contacts you wish to calculate. "
                   "Keywords like GROUPA1, GROUPA2, GROUPA3,... should be listed and one contact will be "
                   "calculated for each ATOM keyword you specify.");
  keys.add("numbered","GROUPB","the atoms involved in each of the contacts you wish to calculate. "
                   "Keywords like GROUPB1, GROUPB2, GROUPB3,... should be listed and one contact will be "
                   "calculated for each ATOM keyword you specify.");
  keys.reset_style("GROUPA","atoms");
  keys.reset_style("GROUPB","atoms");
  keys.addFlag("ADDEXP",false,"Set to TRUE if you want to have fixed components with the experimental reference values.");  
  keys.add("numbered","NOEDIST","Add an experimental value for each NOE.");
  keys.addOutputComponent("noe","default","the # NOE");
  keys.addOutputComponent("exp","ADDEXP","the # NOE experimental distance");
}

NOE::NOE(const ActionOptions&ao):
PLUMED_COLVAR_INIT(ao),
pbc(true)
{
  bool nopbc=!pbc;
  parseFlag("NOPBC",nopbc);
  pbc=!nopbc;  

  // Read in the atoms
  vector<AtomNumber> t, ga_lista, gb_lista;
  for(int i=1;;++i ){
     parseAtomList("GROUPA", i, t );
     if( t.empty() ) break;
     for(unsigned j=0;j<t.size();j++) ga_lista.push_back(t[j]);
     nga.push_back(t.size());
     t.resize(0); 
  }
  vector<unsigned> ngb;
  for(int i=1;;++i ){
     parseAtomList("GROUPB", i, t );
     if( t.empty() ) break;
     for(unsigned j=0;j<t.size();j++) gb_lista.push_back(t[j]);
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

    for(unsigned i=0;i<nga.size();++i){
       if( !parseNumbered( "NOEDIST", i+1, noedist[i] ) ) break;
       ntarget++; 
    }
    if( ntarget!=nga.size() ) error("found wrong number of NOEDIST values");
  }

  // Ouput details of all contacts
  unsigned index=0; 
  for(unsigned i=0;i<nga.size();++i){
    log.printf("  The %uth NOE is calculated using %u equivalent couples of atoms\n", i, nga[i]);
    for(unsigned j=0;j<nga[i];j++) {
      log.printf("    couple %u is %d %d.\n", j, ga_lista[index].serial(), gb_lista[index].serial() );
      index++;
    }
  }

  if(pbc)      log.printf("  using periodic boundary conditions\n");
  else         log.printf("  without periodic boundary conditions\n");

  for(unsigned i=0;i<nga.size();i++) {
    string num; Tools::convert(i,num);
    addComponentWithDerivatives("noe_"+num);
    componentIsNotPeriodic("noe_"+num);
  }

  if(addexp) {
    for(unsigned i=0;i<nga.size();i++) {
      string num; Tools::convert(i,num);
      addComponent("exp_"+num);
      componentIsNotPeriodic("exp_"+num);
      Value* comp=getPntrToComponent("exp_"+num);
      comp->set(noedist[i]);
    }
  }

  requestAtoms(nl->getFullAtomList());
  checkRead();
}

NOE::~NOE(){
  delete nl;
} 

void NOE::calculate()
{
  const unsigned ngasz=nga.size();

  #pragma omp parallel for num_threads(OpenMP::getNumThreads()) 
  for(unsigned i=0;i<ngasz;i++) {
    Tensor dervir;
    double noe=0;
    unsigned index=0;
    for(unsigned k=0; k<i; k++) index+=nga[k]; 
    const double c_aver=1./static_cast<double>(nga[i]);
    Value* val=getPntrToComponent(i);
    // cycle over equivalent atoms 
    for(unsigned j=0;j<nga[i];j++) {
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

}
}
