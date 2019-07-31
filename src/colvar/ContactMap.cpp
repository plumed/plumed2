/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2019 The plumed team
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
#include "tools/NeighborList.h"
#include "ActionRegister.h"
#include "tools/SwitchingFunction.h"

#include <string>
#include <cmath>
#include <memory>

using namespace std;

namespace PLMD {
namespace colvar {

//+PLUMEDOC COLVAR CONTACTMAP
/*
Calculate the distances between a number of pairs of atoms and transform each distance by a switching function.

The transformed distance can be compared with a reference value in order to calculate the squared distance
between two contact maps. Each distance can also be weighted for a given value. CONTACTMAP can be used together
with \ref FUNCPATHMSD to define a path in the contactmap space.

The individual contact map distances related to each contact can be accessed as components
named `cm.contact-1`, `cm.contact-2`, etc, assuming that the label of the CONTACTMAP is `cm`.

\par Examples

The following example calculates switching functions based on the distances between atoms
1 and 2, 3 and 4 and 4 and 5. The values of these three switching functions are then output
to a file named colvar.

\plumedfile
CONTACTMAP ATOMS1=1,2 ATOMS2=3,4 ATOMS3=4,5 ATOMS4=5,6 SWITCH={RATIONAL R_0=1.5} LABEL=f1
PRINT ARG=f1.* FILE=colvar
\endplumedfile

The following example calculates the difference of the current contact map with respect
to a reference provided. In this case REFERENCE is the fraction of contact that is formed
(i.e. the distance between two atoms transformed with the SWITCH), while R_0 is the contact
distance. WEIGHT gives the relative weight of each contact to the final distance measure.

\plumedfile
CONTACTMAP ...
ATOMS1=1,2 REFERENCE1=0.1 WEIGHT1=0.5
ATOMS2=3,4 REFERENCE2=0.5 WEIGHT2=1.0
ATOMS3=4,5 REFERENCE3=0.25 WEIGHT3=1.0
ATOMS4=5,6 REFERENCE4=0.0 WEIGHT4=0.5
SWITCH={RATIONAL R_0=1.5}
LABEL=cmap
CMDIST
... CONTACTMAP

PRINT ARG=cmap FILE=colvar
\endplumedfile

The next example calculates calculates fraction of native contacts (Q)
for Trp-cage mini-protein. R_0 is the distance at which the switch function is guaranteed to
be 1.0 – it doesn't really matter for Q and  should be something very small, like 1 A.
REF is the reference distance for the contact, e.g. the distance from a crystal structure.
LAMBDA is the tolerance for the distance – if set to 1.0, the contact would have to have exactly
the reference value to be formed; instead for lambda values of 1.5–1.8 are usually used to allow some slack.
BETA is the softness of the switch function, default is 50nm.
WEIGHT is the 1/(number of contacts) giving equal weight to each contact.

When using native contact Q switch function, please cite \cite best2013

\plumedfile
# The full (much-longer) example available in regtest/basic/rt72/

CONTACTMAP ...
ATOMS1=1,67 SWITCH1={Q R_0=0.01 BETA=50.0 LAMBDA=1.5 REF=0.4059} WEIGHT1=0.003597
ATOMS2=1,68 SWITCH2={Q R_0=0.01 BETA=50.0 LAMBDA=1.5 REF=0.4039} WEIGHT2=0.003597
ATOMS3=1,69 SWITCH3={Q R_0=0.01 BETA=50.0 LAMBDA=1.5 REF=0.3215} WEIGHT3=0.003597
ATOMS4=5,61 SWITCH4={Q R_0=0.01 BETA=50.0 LAMBDA=1.5 REF=0.4277} WEIGHT4=0.003597
ATOMS5=5,67 SWITCH5={Q R_0=0.01 BETA=50.0 LAMBDA=1.5 REF=0.3851} WEIGHT5=0.003597
ATOMS6=5,68 SWITCH6={Q R_0=0.01 BETA=50.0 LAMBDA=1.5 REF=0.3811} WEIGHT6=0.003597
ATOMS7=5,69 SWITCH7={Q R_0=0.01 BETA=50.0 LAMBDA=1.5 REF=0.3133} WEIGHT7=0.003597
LABEL=cmap
SUM
... CONTACTMAP

PRINT ARG=cmap FILE=colvar
\endplumedfile
(See also \ref switchingfunction)

*/
//+ENDPLUMEDOC

class ContactMap : public Colvar {
private:
  bool pbc, serial, docomp, dosum, docmdist;
  std::unique_ptr<NeighborList> nl;
  std::vector<SwitchingFunction> sfs;
  vector<double> reference, weight;
public:
  static void registerKeywords( Keywords& keys );
  explicit ContactMap(const ActionOptions&);
// active methods:
  void calculate() override;
  void checkFieldsAllowed() override {}
};

PLUMED_REGISTER_ACTION(ContactMap,"CONTACTMAP")

void ContactMap::registerKeywords( Keywords& keys ) {
  Colvar::registerKeywords( keys );
  keys.add("numbered","ATOMS","the atoms involved in each of the contacts you wish to calculate. "
           "Keywords like ATOMS1, ATOMS2, ATOMS3,... should be listed and one contact will be "
           "calculated for each ATOM keyword you specify.");
  keys.reset_style("ATOMS","atoms");
  keys.add("numbered","SWITCH","The switching functions to use for each of the contacts in your map. "
           "You can either specify a global switching function using SWITCH or one "
           "switching function for each contact. Details of the various switching "
           "functions you can use are provided on \\ref switchingfunction.");
  keys.add("numbered","REFERENCE","A reference value for a given contact, by default is 0.0 "
           "You can either specify a global reference value using REFERENCE or one "
           "reference value for each contact.");
  keys.add("numbered","WEIGHT","A weight value for a given contact, by default is 1.0 "
           "You can either specify a global weight value using WEIGHT or one "
           "weight value for each contact.");
  keys.reset_style("SWITCH","compulsory");
  keys.addFlag("SUM",false,"calculate the sum of all the contacts in the input");
  keys.addFlag("CMDIST",false,"calculate the distance with respect to the provided reference contact map");
  keys.addFlag("SERIAL",false,"Perform the calculation in serial - for debug purpose");
  keys.addOutputComponent("contact","default","By not using SUM or CMDIST each contact will be stored in a component");
}

ContactMap::ContactMap(const ActionOptions&ao):
  PLUMED_COLVAR_INIT(ao),
  pbc(true),
  serial(false),
  docomp(true),
  dosum(false),
  docmdist(false)
{
  parseFlag("SERIAL",serial);
  parseFlag("SUM",dosum);
  parseFlag("CMDIST",docmdist);
  if(docmdist==true&&dosum==true) error("You cannot use SUM and CMDIST together");
  bool nopbc=!pbc;
  parseFlag("NOPBC",nopbc);
  pbc=!nopbc;

  // Read in the atoms
  std::vector<AtomNumber> t, ga_lista, gb_lista;
  for(int i=1;; ++i ) {
    parseAtomList("ATOMS", i, t );
    if( t.empty() ) break;

    if( t.size()!=2 ) {
      std::string ss; Tools::convert(i,ss);
      error("ATOMS" + ss + " keyword has the wrong number of atoms");
    }
    ga_lista.push_back(t[0]); gb_lista.push_back(t[1]);
    t.resize(0);

    // Add a value for this contact
    std::string num; Tools::convert(i,num);
    if(!dosum&&!docmdist) {addComponentWithDerivatives("contact-"+num); componentIsNotPeriodic("contact-"+num);}
  }
  // Create neighbour lists
  nl.reset(new NeighborList(ga_lista,gb_lista,true,pbc,getPbc()));

  // Read in switching functions
  std::string errors; sfs.resize( ga_lista.size() ); unsigned nswitch=0;
  for(unsigned i=0; i<ga_lista.size(); ++i) {
    std::string num, sw1; Tools::convert(i+1, num);
    if( !parseNumbered( "SWITCH", i+1, sw1 ) ) break;
    nswitch++; sfs[i].set(sw1,errors);
    if( errors.length()!=0 ) error("problem reading SWITCH" + num + " keyword : " + errors );
  }
  if( nswitch==0 ) {
    std::string sw; parse("SWITCH",sw);
    if(sw.length()==0) error("no switching function specified use SWITCH keyword");
    for(unsigned i=0; i<ga_lista.size(); ++i) {
      sfs[i].set(sw,errors);
      if( errors.length()!=0 ) error("problem reading SWITCH keyword : " + errors );
    }
  } else if( nswitch!=sfs.size()  ) {
    std::string num; Tools::convert(nswitch+1, num);
    error("missing SWITCH" + num + " keyword");
  }

  // Read in reference values
  nswitch=0;
  reference.resize(ga_lista.size(), 0.);
  for(unsigned i=0; i<ga_lista.size(); ++i) {
    if( !parseNumbered( "REFERENCE", i+1, reference[i] ) ) break;
    nswitch++;
  }
  if( nswitch==0 ) {
    parse("REFERENCE",reference[0]);
    for(unsigned i=1; i<ga_lista.size(); ++i) {
      reference[i]=reference[0];
      nswitch++;
    }
  }
  if(nswitch == 0 && docmdist) error("with CMDIST one must use REFERENCE to setup the reference contact map");

  // Read in weight values
  nswitch=0;
  weight.resize(ga_lista.size(), 1.0);
  for(unsigned i=0; i<ga_lista.size(); ++i) {
    if( !parseNumbered( "WEIGHT", i+1, weight[i] ) ) break;
    nswitch++;
  }
  if( nswitch==0 ) {
    parse("WEIGHT",weight[0]);
    for(unsigned i=1; i<ga_lista.size(); ++i) {
      weight[i]=weight[0];
    }
    nswitch = ga_lista.size();
  }

  // Ouput details of all contacts
  for(unsigned i=0; i<sfs.size(); ++i) {
    log.printf("  The %uth contact is calculated from atoms : %d %d. Inflection point of switching function is at %s. Reference contact value is %f\n",
               i+1, ga_lista[i].serial(), gb_lista[i].serial(), ( sfs[i].description() ).c_str(), reference[i] );
  }

  if(dosum) {
    addValueWithDerivatives(); setNotPeriodic();
    log.printf("  colvar is sum of all contacts in contact map\n");
  }
  if(docmdist) {
    addValueWithDerivatives(); setNotPeriodic();
    log.printf("  colvar is distance between the contact map matrix and the provided reference matrix\n");
  }

  if(dosum || docmdist) {
    docomp=false;
  } else {
    serial=true;
    docomp=true;
  }

  // Set up if it is just a list of contacts
  requestAtoms(nl->getFullAtomList());
  checkRead();
}

void ContactMap::calculate() {

  double ncoord=0.;
  Tensor virial;
  std::vector<Vector> deriv(getNumberOfAtoms());

  unsigned stride;
  unsigned rank;
  if(serial) {
    // when using components the parallelisation do not work
    stride=1;
    rank=0;
  } else {
    stride=comm.Get_size();
    rank=comm.Get_rank();
  }

// sum over close pairs
  for(unsigned i=rank; i<nl->size(); i+=stride) {
    Vector distance;
    unsigned i0=nl->getClosePair(i).first;
    unsigned i1=nl->getClosePair(i).second;
    if(pbc) {
      distance=pbcDistance(getPosition(i0),getPosition(i1));
    } else {
      distance=delta(getPosition(i0),getPosition(i1));
    }

    double dfunc=0.;
    double coord = weight[i]*(sfs[i].calculate(distance.modulo(), dfunc) - reference[i]);
    Vector tmpder = weight[i]*dfunc*distance;
    Tensor tmpvir = weight[i]*dfunc*Tensor(distance,distance);
    if(!docmdist) {
      deriv[i0] -= tmpder;
      deriv[i1] += tmpder;
      virial    -= tmpvir;
      ncoord    += coord;
    } else {
      tmpder *= 2.*coord;
      tmpvir *= 2.*coord;
      deriv[i0] -= tmpder;
      deriv[i1] += tmpder;
      virial    -= tmpvir;
      ncoord    += coord*coord;
    }

    if(docomp) {
      Value* val=getPntrToComponent( i );
      setAtomsDerivatives( val, i0, deriv[i0] );
      setAtomsDerivatives( val, i1, deriv[i1] );
      setBoxDerivatives( val, -tmpvir );
      val->set(coord);
    }
  }

  if(!serial) {
    comm.Sum(&ncoord,1);
    if(!deriv.empty()) comm.Sum(&deriv[0][0],3*deriv.size());
    comm.Sum(&virial[0][0],9);
  }

  if( !docomp ) {
    for(unsigned i=0; i<deriv.size(); ++i) setAtomsDerivatives(i,deriv[i]);
    setValue           (ncoord);
    setBoxDerivatives  (virial);
  }
}

}
}
