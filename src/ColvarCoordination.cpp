#include "SwitchingFunction.h"
#include "Colvar.h"
#include "ActionRegister.h"
#include "NeighborList.h"
#include "PlumedCommunicator.h"

#include <string>
#include <cmath>
#include <cassert>

using namespace std;

namespace PLMD{

//+PLUMEDOC COLVAR COORDINATION
/**
This keyword can be used to calculate the coordination numbers for atoms in your system. 
We use the following switching function to make the coordination number differentiable:

\f[
s = \frac{ 1 - \left(\frac{r-d_0}{r_0}\right)^n } { 1 - \left(\frac{r-d_0}{r_0}\right)^m }
\f]

\par Examples

The following example instructs plumed to the average coordination number of the atoms in group 1-10 with the atoms in group 20-100.  These coordination numbers count the number of atoms within the group that are within 0.3 nm of the central atom.  A neighbour list is used to make this calculation faster, this neighbour list is updated every 100 steps.
\verbatim
COORDINATION GROUPA=1-10 GROUPB=20-100 R_0=0.3 NL_CUTOFF=0.5 UPDATE=100 
\endverbatim

*/
//+ENDPLUMEDOC
   
class ColvarCoordination : public Colvar {
  bool pbc;
  bool serial;
  NeighborList *nl;
  SwitchingFunction switchingFunction;

  bool reduceListAtNextStep;
  
public:
  ColvarCoordination(const ActionOptions&);
  ~ColvarCoordination();
// active methods:
  virtual void calculate();
  virtual void prepare();
  static void registerKeywords( Keywords& keys );
};

PLUMED_REGISTER_ACTION(ColvarCoordination,"COORDINATION")

void ColvarCoordination::registerKeywords( Keywords& keys ){
  Colvar::registerKeywords(keys);
  keys.addFlag("SERIAL",false,"perform the calcualtion of the coordination number in serial");
  keys.addFlag("PAIR",false,"Evaulate the switching functions for only the 1st element of the 1st group with the first element in the second group etc only");
  keys.addFlag("NLIST",false,"Use a neighbour list to speed up the calculation");
  keys.add("atoms","GROUPA","The list of central atoms for which we are calculating our coordination numbers");
  keys.add("atoms","GROUPB","The list of neighbourhood atoms for which we are using to calculate coordination numbers");
  keys.add("optional","NN","The n parameter of the switching function ");
  keys.add("optional","MM","The m parameter of the switching function ");
  keys.add("optional","D_0","The d_0 parameter of the switching function");
  keys.add("optional","R_0","The r_0 parameter of the switching function");
  keys.add("optional","SWITCH","A generic switching function");
  keys.add("optional","NL_CUTOFF","The cutoff for the neighbour list");
  keys.add("optional","NL_STRIDE","The frequency with which we are updating the atoms in the neighbour list");
}

ColvarCoordination::ColvarCoordination(const ActionOptions&ao):
PLUMED_COLVAR_INIT(ao),
pbc(true),
serial(false),
reduceListAtNextStep(false)
{

  parseFlag("SERIAL",serial);

  vector<AtomNumber> ga_lista,gb_lista;
  parseAtomList("GROUPA",ga_lista);
  parseAtomList("GROUPB",gb_lista);

  string sw,errors;
  parse("SWITCH",sw);
  if(sw.length()>0) switchingFunction.set(sw,errors);
  else {
    int nn=6;
    int mm=12;
    double d0=0.0;
    double r0=0.0;
    parse("R_0",r0);
    plumed_massert(r0>0.0,"R_0 is compulsory");
    parse("D_0",r0);
    parse("NN",nn);
    parse("MM",mm);
    switchingFunction.set(nn,mm,r0,d0);
  }
  
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
   plumed_assert(nl_cut>0.);
   parse("NL_STRIDE",nl_st);
   plumed_assert(nl_st>0);
  }
  
  checkRead();

  addValueWithDerivatives(); setNotPeriodic();
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
 if(reduceListAtNextStep){
   requestAtoms(nl->getReducedAtomList());
   reduceListAtNextStep=false;
 }
 if(nl->getStride()>0 && (getStep()-nl->getLastUpdate())>=nl->getStride()){
  requestAtoms(nl->getFullAtomList());
 }
}

// calculator
void ColvarCoordination::calculate()
{

 double ncoord=0.;
 Tensor virial;
 vector<Vector> deriv(getNumberOfAtoms());
// deriv.resize(getPositions().size());

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
   distance=pbcDistance(getPosition(i0),getPosition(i1));
  } else {
   distance=delta(getPosition(i0),getPosition(i1));
  }

  double dfunc=0.;
  ncoord += switchingFunction.calculate(distance.modulo(), dfunc);

  deriv[i0] = deriv[i0] + (-dfunc)*distance ;
  deriv[i1] = deriv[i1] + dfunc*distance ;
  virial=virial+(-dfunc)*Tensor(distance,distance);
 }

 if(!serial){
   comm.Sum(&ncoord,1);
   if(deriv.size()>0) comm.Sum(&deriv[0][0],3*deriv.size());
   comm.Sum(&virial[0][0],9);
 }

 for(unsigned i=0;i<deriv.size();++i) setAtomsDerivatives(i,deriv[i]);
 setValue           (ncoord);
 setBoxDerivatives  (virial);

 if(nl->getStride()>0 && (getStep()-nl->getLastUpdate())>=nl->getStride()){
  reduceListAtNextStep=true;
  nl->setLastUpdate(getStep());
 }

}
}
