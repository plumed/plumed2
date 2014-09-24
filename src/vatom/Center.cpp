/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2014 The plumed team
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
#include "ActionWithVirtualAtom.h"
#include "ActionRegister.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"

using namespace std;

namespace PLMD{
namespace vatom{

//+PLUMEDOC VATOM CENTER
/*
Calculate the center for a group of atoms, with arbitrary weights.

The computed
center is stored as a virtual atom that can be accessed in
an atom list through the label for the CENTER action that creates it.
Notice that the generated virtual atom has charge equal to the sum of the
charges and mass equal to the sum of the masses. If used with the MASS flag,
then it provides a result identical to \ref COM.

\par Examples

\verbatim
# a point which is on the line connecting atoms 1 and 10, so that its distance
# from 10 is twice its distance from 1:
c1: CENTER ATOMS=1,1,10
# this is another way of stating the same:
c1bis: CENTER ATOMS=1,10 WEIGHTS=2,1

# center of mass among these atoms:
c2: CENTER ATOMS=2,3,4,5 MASS

d1: DISTANCE ATOMS=c1,c2

PRINT ARG=d1
\endverbatim
(See also \ref DISTANCE, \ref COM and \ref PRINT).

*/
//+ENDPLUMEDOC


class Center:
  public ActionWithVirtualAtom
{
  std::vector<double> weights;
  bool weight_mass;
public:
  Center(const ActionOptions&ao);
  void calculate();
  static void registerKeywords( Keywords& keys );
};

PLUMED_REGISTER_ACTION(Center,"CENTER")

void Center::registerKeywords(Keywords& keys){
  ActionWithVirtualAtom::registerKeywords(keys);
  keys.add("optional","WEIGHTS","Center is computed as a weighted average.");
  keys.addFlag("MASS",false,"If set center is mass weighted");
}

Center::Center(const ActionOptions&ao):
  Action(ao),
  ActionWithVirtualAtom(ao),
  weight_mass(false)
{
  vector<AtomNumber> atoms;
  parseAtomList("ATOMS",atoms);
  if(atoms.size()==0) error("at least one atom should be specified");
  parseVector("WEIGHTS",weights);
  parseFlag("MASS",weight_mass);
  checkRead();
  log.printf("  of atoms");
  for(unsigned i=0;i<atoms.size();++i) log.printf(" %d",atoms[i].serial());
  if(weight_mass){
    log<<"  mass weighted\n";
    if(weights.size()!=0) error("WEIGHTS and MASS keywords should not be used simultaneously");
  } else {
    if( weights.size()==0) {
      weights.resize( atoms.size() );
      for(unsigned i=0;i<atoms.size();i++) weights[i] = 1.;
    }
    log<<" with weights";
    if( weights.size()!=atoms.size() ) error("number of elements in weight vector does not match the number of atoms");
    for(unsigned i=0;i<weights.size();++i) log.printf(" %f",weights[i]);
    log.printf("\n");
  }
  requestAtoms(atoms);
}

void Center::calculate(){
  Vector pos;
  double mass(0.0);
  vector<Tensor> deriv(getNumberOfAtoms());
  for(unsigned i=0;i<getNumberOfAtoms();i++) mass+=getMass(i);
  if( plumed.getAtoms().chargesWereSet() ){
     double charge(0.0);
     for(unsigned i=0;i<getNumberOfAtoms();i++) charge+=getCharge(i);
     setCharge(charge);
  } else {
     setCharge(0.0);
  }
  double wtot=0.0;
  for(unsigned i=0;i<weights.size();i++) wtot+=weights[i];
  for(unsigned i=0;i<getNumberOfAtoms();i++){
    double w=0;
    if(weight_mass) w=getMass(i)/mass;
    else w=weights[i]/wtot;
    pos+=w*getPosition(i);
    deriv[i]=w*Tensor::identity();
  }
  setPosition(pos);
  setMass(mass);
  setAtomsDerivatives(deriv);
}

}
}
