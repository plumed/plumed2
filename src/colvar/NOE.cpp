/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2014,2015 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

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

//+PLUMEDOC COLVAR NOE 
/*
Calculates the deviation of current distances from experimental NOE derived distances.

NOE distances are calculated between couple of atoms, averaging over equivalent couples, and compared with a set of reference distances. 
Distances can also be averaged over multiple replicas to perform replica-averaged simulations.
Each NOE is defined by two groups containing the same number of atoms and by a reference distance, distances
are calculated in pairs.

\f[
NOE() = \sum_i^{noes}((\frac{1}{N_{eq}}\sum_j^{N_{eq}} (\frac{1}{r_j^6}))^{\frac{-1}{6}} - d_i^{exp})^2 
\f]

Reference distances can also be considered as upper limits only, in this case the sum is over a half
parabola. 

\par Examples
In the following examples three noes are defined, the first is calculated based on the distances
of atom 1-2 and 3-2; the second is defined by the distance 5-7 and the third by the distances
4-15,4-16,8-15,8-16.

\verbatim
NOE ...
GROUPA1=1,3 GROUPB1=2,2 NOEDIST1=0.5
GROUPA2=5 GROUPB2=7 NOEDIST2=0.4
GROUPA3=4,4,8,8 GROUPB3=15,16,15,16 NOEDIST3=0.3
LABEL=noes
... NOE

PRINT ARG=noes FILE=colvar
\endverbatim
(See also \ref PRINT) 

*/
//+ENDPLUMEDOC

class NOE : public Colvar {   
private:
  bool             pbc;
  vector<double>   noedist;
  vector<unsigned> nga, ngb;
  NeighborList     *nl;
  unsigned         ens_dim;
  unsigned         pperiod;
  bool             isupper;
  bool             ensemble;
  bool             serial;
public:
  static void registerKeywords( Keywords& keys );
  explicit NOE(const ActionOptions&);
  ~NOE();
  virtual void calculate();
};

PLUMED_REGISTER_ACTION(NOE,"NOE")

void NOE::registerKeywords( Keywords& keys ){
  Colvar::registerKeywords( keys );
  keys.add("numbered","GROUPA","the atoms involved in each of the contacts you wish to calculate. "
                   "Keywords like GROUPA1, GROUPA2, GROUPA3,... should be listed and one contact will be "
                   "calculated for each ATOM keyword you specify.");
  keys.add("numbered","GROUPB","the atoms involved in each of the contacts you wish to calculate. "
                   "Keywords like GROUPB1, GROUPB2, GROUPB3,... should be listed and one contact will be "
                   "calculated for each ATOM keyword you specify.");
  keys.reset_style("GROUPA","atoms");
  keys.reset_style("GROUPB","atoms");
  keys.add("numbered","NOEDIST","A compulsory reference distance for a given NOE"
                                "You can either specify a global reference value using NOEDIST or one "
                                "reference value for each contact.");
  keys.addFlag("UPPER_LIMITS",false,"Set to TRUE if you want to consider the reference distances as upper limits.");
  keys.add("compulsory","WRITE_NOE","0","Write the back-calculated chemical shifts every # steps.");
  keys.addFlag("ENSEMBLE",false,"Set to TRUE if you want to average over multiple replicas.");
  keys.addFlag("SERIAL",false,"Perform the calculation in serial - for debug purpose");
}

NOE::NOE(const ActionOptions&ao):
PLUMED_COLVAR_INIT(ao),
pbc(true),
isupper(false),
ensemble(false),
serial(false)
{
  parseFlag("SERIAL",serial);

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

  // Read in reference values
  noedist.resize( nga.size() ); 
  unsigned ntarget=0;
  for(unsigned i=0;i<nga.size();++i){
     if( !parseNumbered( "NOEDIST", i+1, noedist[i] ) ) break;
     ntarget++; 
  }
  if( ntarget==0 ){
      parse("NOEDIST",noedist[0]);
      for(unsigned i=1;i<nga.size();++i) noedist[i]=noedist[0];
  } else if( ntarget!=nga.size() ) error("found wrong number of NOEDIST values");

  parseFlag("UPPER_LIMITS",isupper);

  unsigned w_period=0;
  parse("WRITE_NOE", w_period);
  pperiod=w_period;

  ensemble=false;
  parseFlag("ENSEMBLE",ensemble);
  if(ensemble){
    if(comm.Get_rank()==0) {
      if(multi_sim_comm.Get_size()<2) error("You CANNOT run Replica-Averaged simulations without running multiple replicas!\n");
      ens_dim=multi_sim_comm.Get_size();
    } else ens_dim=0;
    comm.Sum(&ens_dim, 1);
  } else ens_dim=1;

  // Ouput details of all contacts
  unsigned index=0; 
  for(unsigned i=0;i<nga.size();++i){
    log.printf("  The %uth NOE is calculated using %u equivalent couples of atoms and compared with a %f reference distance\n", i, nga[i], noedist[i]);
    for(unsigned j=0;j<nga[i];j++) {
      log.printf("    couple %u is %d %d.\n", j, ga_lista[index].serial(), gb_lista[index].serial() );
      index++;
    }
  }

  if(isupper)  log.printf("  NOEs reference distances are considered as upper limits only\n");
  if(!serial)  log.printf("  The NOEs are calculated in parallel\n");
  else         log.printf("  The NOEs are calculated in serial\n");
  if(ensemble) log.printf("  ENSEMBLE averaging over %u replicas\n", ens_dim);
  if(pbc)      log.printf("  using periodic boundary conditions\n");
  else         log.printf("  without periodic boundary conditions\n");

  addValueWithDerivatives(); 
  setNotPeriodic(); 
  requestAtoms(nl->getFullAtomList());
  checkRead();
}

NOE::~NOE(){
  delete nl;
} 

void NOE::calculate(){ 
  Tensor virial;
  double score=0.;
  std::vector<Vector> deriv(getNumberOfAtoms());
  unsigned sga = nga.size();
  vector<double> noe(sga);
  vector<double> dnoe(sga);
 
  // internal parallelisation
  unsigned stride=comm.Get_size();
  unsigned rank=comm.Get_rank();
  if(serial){
    stride=1;
    rank=0;
  }
 
  for(unsigned i=0;i<sga;i++) { noe[i]=0.; dnoe[i]=0.;}

  unsigned index=0; for(unsigned k=0;k<rank;k++) index += nga[k];

  for(unsigned i=rank;i<sga;i+=stride) { //cycle over the number of noe 
    for(unsigned j=0;j<nga[i];j++) {
      Vector distance;
      double aver=1./((double)nga[i]);
      unsigned i0=nl->getClosePair(index).first;
      unsigned i1=nl->getClosePair(index).second;
      if(pbc){
        distance=pbcDistance(getPosition(i0),getPosition(i1));
      } else {
        distance=delta(getPosition(i0),getPosition(i1));
      }
      double d=distance.modulo();
      double r2=d*d;
      double r3=r2*d;
      double r6=r3*r3;
      double r8=r6*r2;
      double tmpir6=aver/r6;
      double tmpir8=aver/r8;
      noe[i] += tmpir6;
      deriv[i0] = -tmpir8*distance;
      deriv[i1] = +tmpir8*distance;
      index++;
    }
    for(unsigned k=i+1;k<i+stride;k++) index += nga[k];
    if(!ensemble) {
      double diff = pow(noe[i],(-1./6.)) - noedist[i];
      bool doscore = (isupper&&diff>0.) || (!isupper); 
      if(doscore) {
        score += diff*diff;
        dnoe[i] = 2.*diff/pow(noe[i],(7./6.));
      }
    }
  }

  bool printout=false;
  if(pperiod>0) printout = (!(getStep()%pperiod));
  if(printout) {
    // share the calculated noe
    if(!serial) comm.Sum(&noe[0],noe.size());
    // print only if master
    if(comm.Get_rank()==0) {
      char tmp1[21]; sprintf(tmp1, "%ld", getStep()); 
      string csfile = string("noe")+"-"+getLabel()+"-"+tmp1+string(".dat");
      FILE *outfile = fopen(csfile.c_str(), "w");
      fprintf(outfile, "#index calc exp\n");
      for(unsigned i=0;i<sga;i++) { 
        fprintf(outfile," %4u %10.6f %10.6f\n", i, pow(noe[i],(-1./6.)), noedist[i]);
      }
      fclose(outfile);
    }
  }

  // Ensemble averaging
  double fact=1.0;
  if(ensemble) {
    fact = 1./((double) ens_dim);
    // share the calculated noe unless they have been already shared by printout
    if(!serial&&!printout) comm.Sum(&noe[0],noe.size());
    // I am the master of my replica
    if(comm.Get_rank()==0) {
      // among replicas
      multi_sim_comm.Sum(&noe[0], sga );
      for(unsigned i=0;i<sga;i++) noe[i] *= fact; 
    } else for(unsigned i=0;i<sga;i++) noe[i] = 0.;
    // inside each replica
    comm.Sum(&noe[0], sga );
    for(unsigned i=rank;i<sga;i+=stride) {
      double diff = pow(noe[i],(-1./6.)) - noedist[i];
      bool doscore = (isupper&&diff>0.) || (!isupper); 
      if(doscore) {
        score += diff*diff;
        dnoe[i] = 2.*diff/pow(noe[i],(7./6.));
      }
    }
  }

  index=0; for(unsigned k=0;k<rank;k++) index += nga[k];

  for(unsigned i=rank;i<sga;i+=stride) { //cycle over the number of groups
    for(unsigned j=0;j<nga[i];j++) {
      unsigned i0=nl->getClosePair(index).first;
      unsigned i1=nl->getClosePair(index).second;
      deriv[i0] = fact*dnoe[i]*deriv[i0];
      deriv[i1] = fact*dnoe[i]*deriv[i1];
      virial=virial+(-1.*Tensor(getPosition(i0),deriv[i0]));
      virial=virial+(-1.*Tensor(getPosition(i1),deriv[i1]));
      index++;
    }
    for(unsigned k=i+1;k<i+stride;k++) index += nga[k];
  }

  if(!serial){
    comm.Sum(score);
    comm.Sum(&deriv[0][0],3*deriv.size());
    comm.Sum(virial);
  }

  for(unsigned i=0;i<deriv.size();++i) setAtomsDerivatives(i,deriv[i]);
  setValue           (score);
  setBoxDerivatives  (virial);
}

}
}
