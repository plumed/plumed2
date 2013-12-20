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
#include "tools/SwitchingFunction.h"
#include "ActionRegister.h"

#include <string>

using namespace std;

namespace PLMD{
namespace colvar{

//+PLUMEDOC COLVAR COORDINATION
/*
Calculate coordination numbers.

This keyword can be used to calculate the coordination numbers for atoms in your system. 
We use the following switching function to make the coordination number differentiable:

\f[
s = \frac{ 1 - \left(\frac{r-d_0}{r_0}\right)^n } { 1 - \left(\frac{r-d_0}{r_0}\right)^m }
\f]

To make your calculation faster you can use a neighbor list, which makes it that only a
relevant subset of the pairwise distance are calculated at every step.

\par Examples

The following example instructs plumed to calculate the total coordination number of the atoms in group 1-10 with the atoms in group 20-100.  For atoms 1-10 coordination numbers are calculated that count the number of atoms from the second group that are within 0.3 nm of the central atom.  A neighbour list is used to make this calculation faster, this neighbour list is updated every 100 steps.
\verbatim
COORDINATION GROUPA=1-10 GROUPB=20-100 R_0=0.3 NLIST NL_CUTOFF=0.5 NL_STRIDE=100 
\endverbatim

*/
//+ENDPLUMEDOC
   
class Coordination : public CoordinationBase{
  SwitchingFunction switchingFunction;

public:
  Coordination(const ActionOptions&);
// active methods:
  static void registerKeywords( Keywords& keys );
  virtual double pairing(double distance,double&dfunc,unsigned i,unsigned j)const;
};

PLUMED_REGISTER_ACTION(Coordination,"COORDINATION")

void Coordination::registerKeywords( Keywords& keys ){
  CoordinationBase::registerKeywords(keys);
  keys.add("compulsory","NN","6","The n parameter of the switching function ");
  keys.add("compulsory","MM","12","The m parameter of the switching function ");
  keys.add("compulsory","D_0","0.0","The d_0 parameter of the switching function");
  keys.add("compulsory","R_0","The r_0 parameter of the switching function");
  keys.add("optional","SWITCH","This keyword is used if you want to employ an alternative to the continuous swiching function defined above. "
                               "The following provides information on the \\ref switchingfunction that are available. " 
                               "When this keyword is present you no longer need the NN, MM, D_0 and R_0 keywords."); 
}

Coordination::Coordination(const ActionOptions&ao):
Action(ao),
CoordinationBase(ao)
{

  string sw,errors;
  parse("SWITCH",sw);
  if(sw.length()>0){
    switchingFunction.set(sw,errors);
    if( errors.length()!=0 ) error("problem reading SWITCH keyword : " + errors );
  } else {
    int nn=6;
    int mm=12;
    double d0=0.0;
    double r0=0.0;
    parse("R_0",r0);
    if(r0<=0.0) error("R_0 should be explicitly specified and positive");
    parse("D_0",d0);
    parse("NN",nn);
    parse("MM",mm);
    switchingFunction.set(nn,mm,r0,d0);
  }
  
  checkRead();
}

double Coordination::pairing(double distance,double&dfunc,unsigned i,unsigned j)const{
  (void) i; // avoid warnings
  (void) j; // avoid warnings
  return switchingFunction.calculateSqr(distance,dfunc);
}

}

}
