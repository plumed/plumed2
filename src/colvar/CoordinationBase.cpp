/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013 The plumed team
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
#include "CoordinationBase.h"
#include "tools/NeighborList.h"
#include "tools/Communicator.h"

#include <string>

using namespace std;

namespace PLMD{
namespace colvar{

void CoordinationBase::registerKeywords( Keywords& keys ){
  Colvar::registerKeywords(keys);
  keys.addFlag("SERIAL",false,"Perform the calculation in serial - for debug purpose");
  keys.addFlag("PAIR",false,"Pair only 1st element of the 1st group with 1st element in the second, etc");
  keys.addFlag("NLIST",false,"Use a neighbour list to speed up the calculation");
  keys.add("optional","NL_CUTOFF","The cutoff for the neighbour list");
  keys.add("optional","NL_STRIDE","The frequency with which we are updating the atoms in the neighbour list");
  keys.add("atoms","GROUPA","First list of atoms");
  keys.add("atoms","GROUPB","Second list of atoms");
}

CoordinationBase::CoordinationBase(const ActionOptions&ao):
PLUMED_COLVAR_INIT(ao),
pbc(true),
serial(false),
invalidateList(true),
firsttime(true)
{

  parseFlag("SERIAL",serial);

  vector<AtomNumber> ga_lista,gb_lista;
  parseAtomList("GROUPA",ga_lista);
  parseAtomList("GROUPB",gb_lista);

  bool nopbc=!pbc;
  parseFlag("NOPBC",nopbc);
  pbc=!nopbc;

// pair stuff
  bool dopair=false;
  parseFlag("PAIR",dopair);

// neighbor list stuff
  bool doneigh=false;
  double nl_cut=0.0;
  int nl_st=0;
  parseFlag("NLIST",doneigh);
  if(doneigh){
   parse("NL_CUTOFF",nl_cut);
   if(nl_cut<=0.0) error("NL_CUTOFF should be explicitly specified and positive");
   parse("NL_STRIDE",nl_st);
   if(nl_st<=0) error("NL_STRIDE should be explicitly specified and positive");
  }
  
  addValueWithDerivatives(); setNotPeriodic();
  if(doneigh)  nl= new NeighborList(ga_lista,gb_lista,dopair,pbc,getPbc(),nl_cut,nl_st);
  else         nl= new NeighborList(ga_lista,gb_lista,dopair,pbc,getPbc());
  
  requestAtoms(nl->getFullAtomList());
 
  log.printf("  between two groups of %d and %d atoms\n",ga_lista.size(),gb_lista.size());
  log.printf("  first group:\n");
  for(unsigned int i=0;i<ga_lista.size();++i){
   if ( (i+1) % 25 == 0 ) log.printf("  \n");
   log.printf("  %d", ga_lista[i].serial());
  }
  log.printf("  \n  second group:\n");
  for(unsigned int i=0;i<gb_lista.size();++i){
   if ( (i+1) % 25 == 0 ) log.printf("  \n");
   log.printf("  %d", gb_lista[i].serial());
  }
  log.printf("  \n");
  if(pbc) log.printf("  using periodic boundary conditions\n");
  else    log.printf("  without periodic boundary conditions\n");
  if(dopair) log.printf("  with PAIR option\n");
  if(doneigh){
   log.printf("  using neighbor lists with\n");
   log.printf("  update every %d steps and cutoff %lf\n",nl_st,nl_cut);
  }
}

CoordinationBase::~CoordinationBase(){
  delete nl;
}

void CoordinationBase::prepare(){
  if(nl->getStride()>0){
    if(firsttime || (getStep()%nl->getStride()==0)){
      requestAtoms(nl->getFullAtomList());
      invalidateList=true;
      firsttime=false;
    }else{
      requestAtoms(nl->getReducedAtomList());
      invalidateList=false;
      if(getExchangeStep()) error("Neighbor lists should be updated on exchange steps - choose a NL_STRIDE which divides the exchange stride!");
    }
    if(getExchangeStep()) firsttime=true;
  }
}

// calculator
void CoordinationBase::calculate()
{

 double ncoord=0.;
 Tensor virial;
 vector<Vector> deriv(getNumberOfAtoms());
// deriv.resize(getPositions().size());

 if(nl->getStride()>0 && invalidateList){
   nl->update(getPositions());
 }

 unsigned stride=comm.Get_size();
 unsigned rank=comm.Get_rank();
 if(serial){
   stride=1;
   rank=0;
 }else{
   stride=comm.Get_size();
   rank=comm.Get_rank();
 }

 for(unsigned int i=rank;i<nl->size();i+=stride) {                   // sum over close pairs
 
  Vector distance;
  unsigned i0=nl->getClosePair(i).first;
  unsigned i1=nl->getClosePair(i).second;
  if(pbc){
   distance=pbcDistance(getPosition(i0),getPosition(i1));
  } else {
   distance=delta(getPosition(i0),getPosition(i1));
  }

  double dfunc=0.;
  ncoord += pairing(distance.modulo2(), dfunc,i0,i1);

  deriv[i0] = deriv[i0] + (-dfunc)*distance ;
  deriv[i1] = deriv[i1] + dfunc*distance ;
  virial=virial+(-dfunc)*Tensor(distance,distance);
 }

 if(!serial){
   comm.Sum(ncoord);
   if(!deriv.empty()) comm.Sum(&deriv[0][0],3*deriv.size());
   comm.Sum(virial);
 }

 for(unsigned i=0;i<deriv.size();++i) setAtomsDerivatives(i,deriv[i]);
 setValue           (ncoord);
 setBoxDerivatives  (virial);

}
}
}
