#include "SwitchingFunction.h"
#include "Colvar.h"
#include "ActionRegister.h"
#include "NeighborList.h"

#include <string>
#include <cmath>
#include <cassert>

using namespace std;

namespace PLMD{

//+PLUMEDOC COLVAR COORDINATION
/**
This is just a template variable

*/
//+ENDPLUMEDOC
   
class ColvarCoordination : public Colvar {
  bool pbc;
  bool serial;
  NeighborList *nl;
  SwitchingFunction switchingFunction;
  
public:
  ColvarCoordination(const ActionOptions&);
  ~ColvarCoordination();
// active methods:
  virtual void calculate();
  virtual void prepare();
};

PLUMED_REGISTER_ACTION(ColvarCoordination,"COORDINATION")

ColvarCoordination::ColvarCoordination(const ActionOptions&ao):
PLUMED_COLVAR_INIT(ao),
pbc(true),
serial(false)
{

  parseFlag("SERIAL",serial);

  vector<AtomNumber> ga_lista,gb_lista;
  parseAtomList("GROUPA",ga_lista);
  parseAtomList("GROUPB",gb_lista);
  
  int nn,mm;
  double r_0,d_0;
  nn=6;
  mm=12;
  r_0=-1;
  d_0=0.0;
  parse("NN",nn);
  parse("MM",mm);
  parse("R_0",r_0);
  parse("D_0",d_0);
  assert(r_0>0); // this is the only compulsory option
  switchingFunction.set(nn,mm,r_0,d_0);

  bool nopbc=!pbc;
  parseFlag("NOPBC",nopbc);
  pbc=!nopbc;
  parseFlag("PBC",pbc);

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
   assert(nl_cut>0.);
   parse("NL_STRIDE",nl_st);
   assert(nl_st>0);
  }
  
  checkRead();

  addValueWithDerivatives("");
  getValue("")->setPeriodicity(false);

  if(doneigh)  nl= new NeighborList(ga_lista,gb_lista,dopair,pbc,getPbc(),nl_cut,nl_st);
  else         nl= new NeighborList(ga_lista,gb_lista,dopair,pbc,getPbc());
  
  requestAtoms(nl->getFullAtomList());
 
  log.printf("  between two groups of %d and %d atoms\n",ga_lista.size(),gb_lista.size());
  log.printf("  first group:\n");
  for(unsigned int i=0;i<ga_lista.size();++i){
   log.printf("  %d", ga_lista[i].serial());
  }
  log.printf("  \n  second group:\n");
  for(unsigned int i=0;i<gb_lista.size();++i){
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

ColvarCoordination::~ColvarCoordination(){
  delete nl;
}

void ColvarCoordination::prepare(){
 if(nl->getStride()>0 && (getStep()-nl->getLastUpdate())>=nl->getStride()){
  requestAtoms(nl->getFullAtomList());
 }
}

// calculator
void ColvarCoordination::calculate()
{

 double ncoord=0.;
 Tensor virial;
 vector<Vector> deriv;
 deriv.resize(getPositions().size());

 if(nl->getStride()>0 && (getStep()-nl->getLastUpdate())>=nl->getStride()){
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
   distance=pbcDistance(getPositions(i0),getPositions(i1));
  } else {
   distance=delta(getPositions(i0),getPositions(i1));
  }

  double dfunc=0.;
  ncoord += switchingFunction.calculate(distance.modulo(), dfunc);

  deriv[i0] = deriv[i0] + (-dfunc)*distance ;
  deriv[i1] = deriv[i1] + dfunc*distance ;
  virial=virial+(-dfunc)*Tensor(distance,distance);
 }

 if(!serial){
   comm.Sum(&ncoord,1);
   comm.Sum(&deriv[0][0],3*deriv.size());
 }

 for(unsigned i=0;i<deriv.size();++i) setAtomsDerivatives(i,deriv[i]);
 setValue           (ncoord);
 setBoxDerivatives  (virial);

 if(nl->getStride()>0 && (getStep()-nl->getLastUpdate())>=nl->getStride()){
  requestAtoms(nl->getReducedAtomList());
  nl->setLastUpdate(getStep());
 }

}
}
