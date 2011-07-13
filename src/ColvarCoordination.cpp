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
  int nn;
  int mm;
  double r_0;
  double d_0;
  vector<Vector> deriv;
  NeighborList *nl;
  
public:
  ColvarCoordination(const ActionOptions&);
// active methods:
  virtual void calculate();
  virtual void prepare();
};

PLUMED_REGISTER_ACTION(ColvarCoordination,"COORDINATION")

ColvarCoordination::ColvarCoordination(const ActionOptions&ao):
PLUMED_COLVAR_INIT(ao),
pbc(true)
{
  vector<AtomNumber> ga_lista,gb_lista;
  parseAtomList("GROUPA",ga_lista);
  parseAtomList("GROUPB",gb_lista);
  
  vector<int> tnn,tmm;
  vector<double> tr_0,td_0;
  parseVector("NN",tnn);
  assert(tnn.size()==1);
  parseVector("MM",tmm);
  assert(tmm.size()==1);
  parseVector("R_0",tr_0);
  assert(tr_0.size()==1);
  parseVector("D_0",td_0);
  assert(td_0.size()==1);
  nn=tnn[0];
  mm=tmm[0];
  r_0=tr_0[0];
  d_0=td_0[0];

  bool nopbc=!pbc;
  parseFlag("NOPBC",nopbc);
  pbc=!nopbc;
  parseFlag("PBC",pbc);

// pair stuff
  bool dopair;
  parseFlag("PAIR",dopair);

// neighbor list stuff
  bool doneigh;
  const Pbc & pbc_c=getPbc();
  vector<double> nl_cut;
  vector<int> nl_st;
  parseFlag("NLIST",doneigh);
  if(doneigh){
   parseVector("NL_CUTOFF",nl_cut);
   assert(nl_cut.size()==1);
   parseVector("NL_STRIDE",nl_st);
   assert(nl_st.size()==1);
  }
  
  checkRead();

  addValueWithDerivatives("");

  if(doneigh)  *nl= NeighborList(ga_lista,gb_lista,dopair,pbc,pbc_c,nl_cut[0],nl_st[0]);
  else         *nl= NeighborList(ga_lista,gb_lista,dopair,pbc,pbc_c);
  
  requestAtoms(nl->getFullAtomList());
  deriv.resize(nl->getFullAtomList().size());

  log.printf("  between two groups of %d and %d atoms\n",ga_lista.size(),gb_lista.size());
  log.printf("  first group:\n");
  for(unsigned int i=0;i<ga_lista.size();++i){
   log.printf("  %d", ga_lista[i].serial());
  }
  log.printf("  \nsecond group:\n");
  for(unsigned int i=0;i<gb_lista.size();++i){
   log.printf("  %d", gb_lista[i].serial());
  }
  log.printf("  \n");
  if(pbc) log.printf("  using periodic boundary conditions\n");
  else    log.printf("  without periodic boundary conditions\n");
  if(dopair) log.printf("  with PAIR option\n");
  if(doneigh){
   log.printf("  using neighbor lists with\n");
   log.printf("  update every %d steps and cutoff %lf\n",nl_st[0],nl_cut[0]);
  }
}

void ColvarCoordination::prepare(){
 if(nl->doUpdate(getStep())) {requestAtoms(nl->getFullAtomList());}
}

// calculator
void ColvarCoordination::calculate()
{

 double ncoord=0.;
 //double threshold=pow(0.00001,1./(nn-mm));


 Tensor virial;
 deriv.resize(getPositions().size());
 for(unsigned int i=0;i<deriv.size();++i) deriv[i].clear();

 if(nl->doUpdate(getStep())) {nl->update(getPositions(),getStep());}

 for(unsigned int i=0;i<nl->size();++i) {                   // sum over close pairs
  Vector distance;
  unsigned i0=nl->getClosePair(i).first;
  unsigned i1=nl->getClosePair(i).second;
  if(pbc){
   distance=pbcDistance(getPositions(i0),getPositions(i1));
  } else {
   distance=delta(getPositions(i0),getPositions(i1));
  }

// we should define switching functions and derivatives in Tools 
// Carlo: I put an implementation of the switching function in Tools
// It can be called as
  double dfunc=0.;
  ncoord += Tools::switchingFunc(distance.modulo(), nn, mm, r_0, d_0, &dfunc);
/*
  const double rdist = (distance.modulo()-d_0)/r_0;
  double dfunc=0.;
  /// analitic limit of the switching function 
  if(rdist<=0.){
     ncoord+=1.;
     dfunc=0.;
  }else if(rdist>0.999999 && rdist<1.000001){
     ncoord+=nn/mm;
     dfunc=0.5*nn*(nn-mm)/mm;
  }else if(rdist>threshold){
     dfunc=0.;
  }else{
     double rNdist = pow(rdist, nn-1);
     double rMdist = pow(rdist, mm-1);
     double num = 1.-rNdist*rdist;
     double iden = 1./(1.-rMdist*rdist);
     double func = num*iden;
     ncoord += func;
     dfunc = ((-nn*rNdist*iden)+(func*(iden*mm)*rMdist))/(distance.modulo()*r_0);
  }
  */
  deriv[i0] = deriv[i0] + (-dfunc)*distance ;
  deriv[i1] = deriv[i1] + dfunc*distance ;
  virial=virial+(-dfunc)*Tensor(distance,distance);
 }

 for(int i=0;i<deriv.size();++i) setAtomsDerivatives(i,deriv[i]);
 setValue           (ncoord);
 setBoxDerivatives  (virial);

 if(nl->doUpdate(getStep())){requestAtoms(nl->getReducedAtomList());}

}

}



