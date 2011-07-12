#include "Colvar.h"
#include "ActionRegister.h"

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
  int gasize;
  int gbsize;
  int nn;
  int mm;
  double r_0;
  double d_0;

public:
  ColvarCoordination(const ActionOptions&);
// active methods:
  virtual void calculate();
};

PLUMED_REGISTER_ACTION(ColvarCoordination,"COORDINATION")

ColvarCoordination::ColvarCoordination(const ActionOptions&ao):
PLUMED_COLVAR_INIT(ao),
pbc(true)
{
  vector<string> ga;
  vector<string> gb;
  vector<int> lista;
  parseVector("GROUPA",ga); // this is just parsing "commas" ("1-4","14-16")
  Tools::interpretRanges(ga); // this is expanding the array ("1","2","3","4","14","15","16")
  parseVector("GROUPB",gb); // this is just parsing "commas" ("1-4","14-16")
  Tools::interpretRanges(gb); // this is expanding the array ("1","2","3","4","14","15","16")
  gasize=ga.size();
  gbsize=gb.size();
  lista.resize((gasize+gbsize));
  for(int i=0;i<gasize;i++) Tools::convert(ga[i],lista[i]); // this is converting strings to int
  for(int i=gasize;i<(gasize+gbsize);i++) Tools::convert(gb[i-gasize],lista[i]); // this is converting strings to int

  vector<int> tnn;
  vector<int> tmm;
  vector<double> tr_0;
  vector<double> td_0;
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

  checkRead();

  log.printf("  between two groups of %d and %d atoms\n",gasize,gbsize);
  log.printf("  first group:\n");
  for(int i=0;i<gasize;i++) log.printf("  %d", lista[i]);
  log.printf("  \nsecond group:\n");
  for(int i=gasize;i<(gasize+gbsize);i++) log.printf("  %d", lista[i]);
  log.printf("  \n");
  if(pbc) log.printf("  using periodic boundary conditions\n");
  else    log.printf("  without periodic boundary conditions\n");

  addValueWithDerivatives("");

  requestAtoms(lista);
}


// calculator
void ColvarCoordination::calculate(){

 double ncoord=0.;
 double threshold=pow(0.00001,1./(nn-mm));

 for(int i=0;i<gasize;i++) {                                           // sum over CoordNumber(i)
    for(int j=gasize;j<(gasize+gbsize);j++) {
      Vector distance;
      if(pbc){
        distance=pbcDistance(getPositions(i),getPositions(j));
      } else {
        distance=delta(getPositions(i),getPositions(j));
      }
      const double rdist = (distance.modulo()-d_0)/r_0;
      double dfunc=0.;
      /* analitic limit of the switching function */
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
      setAtomsDerivatives(i,dfunc*distance);
      setAtomsDerivatives(j,dfunc*distance);
    }
  }
  setValue           (ncoord);

  //setBoxDerivatives  (-invvalue*Tensor(distance,distance));

}

}



