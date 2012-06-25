#include "Colvar.h"
#include "NeighborList.h"
#include "ActionRegister.h"
#include "SwitchingFunction.h"

#include <string>
#include <cmath>
#include <cassert>

using namespace std;

namespace PLMD {

//+PLUMEDOC MCOLVAR CONTACTMAP
/**

*/
//+ENDPLUMEDOC

class ContactMap : public Colvar {   
private:
  bool pbc, dosum;
  NeighborList *nl;
  std::vector<SwitchingFunction> sfs;
  bool reduceListAtNextStep;
public:
  static void registerKeywords( Keywords& keys );
  ContactMap(const ActionOptions&);
  ~ContactMap();
// active methods:
  virtual void calculate();
  void checkFieldsAllowed(){};
};

PLUMED_REGISTER_ACTION(ContactMap,"CONTACTMAP")

void ContactMap::registerKeywords( Keywords& keys ){
  Colvar::registerKeywords( keys );
  keys.add("numbered","ATOMS","the atoms involved in each of the contacts you wish to calculate. "
                   "Keywords like ATOMS1, ATOMS2, ATOMS3,... should be listed and one contact will be "
                   "calculated for each ATOM keyword you specify.");
  keys.reset_style("ATOMS","atoms");
  keys.add("numbered","SWITCH","The switching functions to use for each of the contacts in your map.  You can either specify a global switching function using SWITCH or one switching function for each contact"); keys.reset_style("SWITCH","compulsory"); 
  keys.addFlag("SUM",false,"calculate the sum of all the contacts in the input");
}

ContactMap::ContactMap(const ActionOptions&ao):
PLUMED_COLVAR_INIT(ao),
pbc(true),
dosum(false),
reduceListAtNextStep(false)
{
  parseFlag("SUM",dosum);
  bool nopbc=!pbc;
  parseFlag("NOPBC",nopbc);
  pbc=!nopbc;  

  // Read in the atoms
  std::vector<AtomNumber> t, ga_lista, gb_lista;
  for(int i=1;;++i ){
     parseAtomList("ATOMS", i, t );
     if( t.empty() ) break;

     log.printf("  The %dth contact is calculated from atoms : ", i);
     for(unsigned j=0;j<t.size();++j) log.printf("%d ",t[j].serial() );
     log.printf("\n");

     if( t.size()!=2 ){
         std::string ss; Tools::convert(i,ss);
         error("ATOMS" + ss + " keyword has the wrong number of atoms");
     }
     ga_lista.push_back(t[0]); gb_lista.push_back(t[1]);
     t.resize(0); 

     // Add a value for this contact
     std::string num; Tools::convert(i,num);
     if(!dosum) addComponentWithDerivatives("contact"+num); componentIsNotPeriodic("contact"+num);
  }
  // Create neighbour lists
  nl= new NeighborList(ga_lista,gb_lista,true,pbc,getPbc());
  requestAtoms(nl->getFullAtomList());

  // Set up if it is just a list of contacts
  if(dosum){
     addValueWithDerivatives(); setNotPeriodic();
     log.printf("  colvar is sum of all contacts in contact map");
  }

  // Read in switching functions
  std::string sw, errors; parse("SWITCH",sw);
  for(unsigned i=0;i<ga_lista.size();++i) sfs.push_back( SwitchingFunction() );
  if(sw.length()>0){
     for(unsigned i=0;i<ga_lista.size();++i){
        sfs[i].set(sw,errors);
        if( errors.length()!=0 ) error("problem reading SWITCH keyword : " + errors );
     }
  } else {
     for(unsigned i=0;i<ga_lista.size();++i){
         std::string num, sw1; Tools::convert(i+1, num);
         if( !parseNumbered( "SWITCH", i+1, sw1 ) ) error("missing SWITCH" + num + " keyword");
         sfs[i].set(sw1,errors);
         if( errors.length()!=0 ) error("problem reading SWITCH" + num + " keyword : " + errors );
     }
  }
  checkRead();
}

ContactMap::~ContactMap(){
  delete nl;
}

void ContactMap::calculate(){ 
     
 double ncoord=0., coord;
 Tensor virial;
 std::vector<Vector> deriv(getNumberOfAtoms());

   for(unsigned int i=0;i<nl->size();i++) {                   // sum over close pairs
      Vector distance;
      unsigned i0=nl->getClosePair(i).first;
      unsigned i1=nl->getClosePair(i).second;
      if(pbc){
       distance=pbcDistance(getPosition(i0),getPosition(i1));
      } else {
       distance=delta(getPosition(i0),getPosition(i1));
      }

      double dfunc=0.;
      coord = sfs[i].calculate(distance.modulo(), dfunc);
      if( !dosum ) {
         Value* val=getPntrToComponent( i );
         setAtomsDerivatives( val, i0, (-dfunc)*distance );
         setAtomsDerivatives( val, i1, dfunc*distance ); 
         setBoxDerivatives( val, (-dfunc)*Tensor(distance,distance) );
         val->set(coord);
      } else {  
         deriv[i0] = deriv[i0] + (-dfunc)*distance ;
         deriv[i1] = deriv[i1] + dfunc*distance ;
         virial=virial+(-dfunc)*Tensor(distance,distance);
         ncoord += coord;
      }
   }

  if( dosum ){
    for(unsigned i=0;i<deriv.size();++i) setAtomsDerivatives(i,deriv[i]);
    setValue           (ncoord);
    setBoxDerivatives  (virial);
  }
}

}
