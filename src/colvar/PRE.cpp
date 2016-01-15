/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015,2016 The plumed team
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

#include <string>
#include <cmath>

using namespace std;

namespace PLMD {
namespace colvar {

//+PLUMEDOC COLVAR PRE 
/*
Calculates the PRE intensity ratio.

*/
//+ENDPLUMEDOC

class PRE : public Colvar {   
private:
  bool             pbc;
  bool             serial;
  double           constant, inept;
  vector<double>   rtwo;
  vector<unsigned> nga;
  NeighborList     *nl;
public:
  static void registerKeywords( Keywords& keys );
  explicit PRE(const ActionOptions&);
  ~PRE();
  virtual void calculate();
};

PLUMED_REGISTER_ACTION(PRE,"PRE")

void PRE::registerKeywords( Keywords& keys ){
  Colvar::registerKeywords( keys );
  keys.add("compulsory","INEPT","is the INEPT time (in ms).");
  keys.add("compulsory","TAUC","is the correlation time (in ns) for this electron-nuclear interaction.");
  keys.add("compulsory","OMEGA","is the Larmor frequency of the nuclear spin (in MHz).");
  keys.add("atoms","SPINLABEL","The atom to be used as the paramagnetic center.");
  keys.add("numbered","GROUPA","the atoms involved in each of the contacts you wish to calculate. "
                   "Keywords like GROUPA1, GROUPA2, GROUPA3,... should be listed and one contact will be "
                   "calculated for each ATOM keyword you specify.");
  keys.reset_style("GROUPA","atoms");
  keys.add("numbered","RTWO","The relaxation of the atom/atoms in the corresponding GROUPA of atoms. "
                   "Keywords like RTWO1, RTWO2, RTWO3,... should be listed.");
  keys.addFlag("SERIAL",false,"Perform the calculation in serial - for debug purpose");
  keys.addOutputComponent("pre","default","the # PRE");
}

PRE::PRE(const ActionOptions&ao):
PLUMED_COLVAR_INIT(ao),
pbc(true),
serial(false)
{
  parseFlag("SERIAL",serial);

  bool nopbc=!pbc;
  parseFlag("NOPBC",nopbc);
  pbc=!nopbc;  

  vector<AtomNumber> atom;
  parseAtomList("SPINLABEL",atom);
  if(atom.size()!=1)
    error("Number of specified atom should be 1");

  // Read in the atoms
  vector<AtomNumber> t, ga_lista, gb_lista;
  for(int i=1;;++i ){
     parseAtomList("GROUPA", i, t );
     if( t.empty() ) break;
     for(unsigned j=0;j<t.size();j++) {ga_lista.push_back(t[j]); gb_lista.push_back(atom[0]);}
     nga.push_back(t.size());
     t.resize(0); 
  }

  // Read in reference values
  rtwo.resize( nga.size() ); 
  unsigned ntarget=0;
  for(unsigned i=0;i<nga.size();++i){
     if( !parseNumbered( "RTWO", i+1, rtwo[i] ) ) break;
     ntarget++; 
  }
  if( ntarget==0 ){
      parse("RTWO",rtwo[0]);
      for(unsigned i=1;i<nga.size();++i) rtwo[i]=rtwo[0];
  } else if( ntarget!=nga.size() ) error("found wrong number of RTWO values");

  double tauc=0.;
  parse("TAUC",tauc);
  if(tauc==0.) error("TAUC must be set");

  double omega=0.;
  parse("OMEGA",omega);
  if(omega==0.) error("OMEGA must be set");

  inept=0.;
  parse("INEPT",inept);
  if(inept==0.) error("INEPT must be set");
  inept *= 0.001; // ms2s

  const double ns2s   = 0.000000001;
  const double MHz2Hz = 1000000.;
  const double Kappa  = 12300000000.00; // this is 1/15*S*(S+1)*gamma^2*g^2*beta^2 
                                        // where gamma is the nuclear gyromagnetic ratio, 
                                        // g is the electronic g factor, and beta is the Bohr magneton
                                        // in nm^6/s^2
  constant = (4.*tauc*ns2s+(3.*tauc*ns2s)/(1+omega*omega*MHz2Hz*MHz2Hz*tauc*tauc*ns2s*ns2s))*Kappa;

  // Create neighbour lists
  nl= new NeighborList(gb_lista,ga_lista,true,pbc,getPbc());

  // Ouput details of all contacts
  unsigned index=0; 
  for(unsigned i=0;i<nga.size();++i){
    log.printf("  The %uth PRE is calculated using %u equivalent atoms:\n", i, nga[i]);
    log.printf("    %d", ga_lista[index].serial());
    index++;
    for(unsigned j=1;j<nga[i];j++) {
      log.printf(" %d", ga_lista[index].serial());
      index++;
    }
    log.printf("\n");
  }

  if(!serial)  log.printf("  The PREs are calculated in parallel\n");
  else         log.printf("  The PREs are calculated in serial\n");
  if(pbc)      log.printf("  using periodic boundary conditions\n");
  else         log.printf("  without periodic boundary conditions\n");

  for(unsigned i=0;i<nga.size();i++) {
    std::string num; Tools::convert(i,num);
    addComponentWithDerivatives("pre_"+num);
    componentIsNotPeriodic("pre_"+num);
  }

  requestAtoms(nl->getFullAtomList());
  checkRead();
}

PRE::~PRE(){
  delete nl;
} 

void PRE::calculate(){ 
  unsigned sga = nga.size();
  std::vector<Vector> deriv(getNumberOfAtoms());
  std::vector<Tensor> dervir(sga);
  vector<double> pre(sga), ratio(sga);

  // internal parallelisation
  unsigned stride=comm.Get_size();
  unsigned rank=comm.Get_rank();
  if(serial){
    stride=1;
    rank=0;
  }
 
  for(unsigned i=0;i<sga;i++) ratio[i]=pre[i]=0.;

  unsigned index=0; for(unsigned k=0;k<rank;k++) index += nga[k];

  for(unsigned i=rank;i<sga;i+=stride) { //cycle over the number of pre 
    for(unsigned j=0;j<nga[i];j++) {
      double aver=1./((double)nga[i]);
      unsigned i0=nl->getClosePair(index).first;
      unsigned i1=nl->getClosePair(index).second;
      Vector distance;
      if(pbc){
        distance=pbcDistance(getPosition(i0),getPosition(i1));
      } else {
        distance=delta(getPosition(i0),getPosition(i1));
      }
      double d=distance.modulo();
      double r2=d*d;
      double r6=r2*r2*r2;
      double r8=r6*r2;
      double tmpir6=constant*aver/r6;
      double tmpir8=-6.*constant*aver/r8;

      pre[i] += tmpir6;

      deriv[i0] = -tmpir8*distance;
      deriv[i1] = +tmpir8*distance;

      dervir[i] += Tensor(distance,deriv[i0]);

      index++;
    }
    ratio[i] = rtwo[i]*exp(-pre[i]*inept) / (rtwo[i]+pre[i]);
    for(unsigned k=i+1;k<i+stride;k++) index += nga[k];
  }

  if(!serial){
    comm.Sum(&pre[0],sga);
    comm.Sum(&ratio[0],sga);
    comm.Sum(&deriv[0][0],3*deriv.size());
    comm.Sum(&dervir[0][0][0],9*sga);
  }

  index = 0;
  for(unsigned i=0;i<sga;i++) { //cycle over the number of pre 
    Value* val=getPntrToComponent(i);
    val->set(ratio[i]);
    double fact = -ratio[i]*(inept+1./(rtwo[i]+pre[i]));
    Tensor virial = fact*dervir[i];
    setBoxDerivatives(val, virial);
    for(unsigned j=0;j<nga[i];j++) {
      unsigned i0=nl->getClosePair(index).first;
      unsigned i1=nl->getClosePair(index).second;
      setAtomsDerivatives(val, i0, fact*deriv[i0]);
      setAtomsDerivatives(val, i1, fact*deriv[i1]); 
      index++;
    }
  }

}

}
}
