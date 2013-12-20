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
#include "Colvar.h"
#include "tools/NeighborList.h"
#include "ActionRegister.h"
#include "tools/SwitchingFunction.h"

#include <string>
#include <cmath>

using namespace std;

namespace PLMD {
namespace colvar {

//+PLUMEDOC COLVAR CONTACTMAP
/*
Calculate the distances between a number of pairs of atoms and transform each distance by a switching function.

\par Examples

The following example calculates switching functions based on the distances between atoms
1 and 2, 3 and 4 and 4 and 5. The values of these three switching functions are then output
to a file named colvar.

\verbatim
CONTACTMAP ATOMS1=1,2 ATOMS2=3,4 ATOMS3=4,5 ATOMS4=5,6 SWITCH=(RATIONAL R_0=1.5) LABEL=f1
PRINT ARG=f1.* FILE=colvar
\endverbatim
(See also \ref PRINT)

*/
//+ENDPLUMEDOC

class ContactMap : public Colvar {   
private:
  bool pbc, dosum;
  NeighborList *nl;
  std::vector<SwitchingFunction> sfs;
public:
  static void registerKeywords( Keywords& keys );
  ContactMap(const ActionOptions&);
  ~ContactMap();
// active methods:
  virtual void calculate();
  void checkFieldsAllowed(){}
};

PLUMED_REGISTER_ACTION(ContactMap,"CONTACTMAP")

void ContactMap::registerKeywords( Keywords& keys ){
  Colvar::registerKeywords( keys );
  keys.add("numbered","ATOMS","the atoms involved in each of the contacts you wish to calculate. "
                   "Keywords like ATOMS1, ATOMS2, ATOMS3,... should be listed and one contact will be "
                   "calculated for each ATOM keyword you specify.");
  keys.reset_style("ATOMS","atoms");
  keys.add("numbered","SWITCH","The switching functions to use for each of the contacts in your map. "
                               "You can either specify a global switching function using SWITCH or one "
                               "switching function for each contact. Details of the various switching "
                               "functions you can use are provided on \\ref switchingfunction."); 
  keys.reset_style("SWITCH","compulsory"); 
  keys.addFlag("SUM",false,"calculate the sum of all the contacts in the input");
}

ContactMap::ContactMap(const ActionOptions&ao):
PLUMED_COLVAR_INIT(ao),
pbc(true),
dosum(false)
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

  // Read in switching functions
  std::string errors; sfs.resize( ga_lista.size() ); unsigned nswitch=0;
  for(unsigned i=0;i<ga_lista.size();++i){
      std::string num, sw1; Tools::convert(i+1, num);
      if( !parseNumbered( "SWITCH", i+1, sw1 ) ) break;
      nswitch++; sfs[i].set(sw1,errors);
      if( errors.length()!=0 ) error("problem reading SWITCH" + num + " keyword : " + errors );
  }
  if( nswitch==0 ){
     std::string sw; parse("SWITCH",sw);
     if(sw.length()==0) error("no switching function specified use SWITCH keyword");
     for(unsigned i=0;i<ga_lista.size();++i){
        sfs[i].set(sw,errors);
        if( errors.length()!=0 ) error("problem reading SWITCH keyword : " + errors );
     }
  } else if( nswitch!=sfs.size()  ){
     std::string num; Tools::convert(nswitch+1, num);
     error("missing SWITCH" + num + " keyword");
  }

  // Ouput details of all contacts 
  for(unsigned i=0;i<sfs.size();++i){
     log.printf("  The %dth contact is calculated from atoms : %d %d. Inflection point of switching function is at %s\n", i+1, ga_lista[i].serial(), gb_lista[i].serial() , ( sfs[i].description() ).c_str() );
  }

  // Set up if it is just a list of contacts
  if(dosum){
     addValueWithDerivatives(); setNotPeriodic();
     log.printf("  colvar is sum of all contacts in contact map");
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
      coord = sfs[i].calculateSqr(distance.modulo2(), dfunc);
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
}
