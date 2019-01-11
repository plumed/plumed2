/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2019 The plumed team
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
#include "core/ActionWithVirtualAtom.h"
#include "core/ActionRegister.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "MultiColvarBase.h"
#include "CatomPack.h"
#include "BridgedMultiColvarFunction.h"
#include "vesselbase/StoreDataVessel.h"

using namespace std;

//+PLUMEDOC VATOM CENTER_OF_MULTICOLVAR
/*
Calculate a a weighted average position based on the value of some multicolvar.

This action calculates the position of a new virtual atom using the following formula:

\f[
x_\alpha = \frac{1}{2\pi} \arctan \left[ \frac{ \sum_i w_i f_i \sin\left( 2\pi x_{i,\alpha} \right) }{ \sum_i w_i f_i \cos\left( 2\pi x_{i,\alpha} \right) } \right]
\f]

Where in this expression the \f$w_i\f$ values are a set of weights calculated within a multicolvar
action and the \f$f_i\f$ are the values of the multicolvar functions. The \f$x_{i,\alpha}\f$ values are
the positions (in scaled coordinates) associated with each of the multicolvars calculated.

\bug The virial contribution for this type of virtual atom is not currently evaluated so do not use in bias functions unless the volume of the cell is fixed

\par Examples

Lets suppose that you are examining the formation of liquid droplets from gas.  You may want to
determine the center of mass of any of the droplets formed.  In doing this calculation you recognize that
the atoms in the liquid droplets will have a higher coordination number than those in the surrounding gas.
As you want to calculate the position of the droplets you thus recognize that these atoms with high coordination
numbers should have a high weight in the weighted average you are using to calculate the position of the droplet.
You can thus calculate the position of the droplet using an input like the one shown below:

\plumedfile
c1: COORDINATIONNUMBER SPECIES=1-512 SWITCH={EXP D_0=4.0 R_0=0.5}
cc: CENTER_OF_MULTICOLVAR DATA=c1
\endplumedfile

The first line here calculates the coordination numbers of all the atoms in the system.  The virtual atom then uses the values
of the coordination numbers calculated by the action labelled c1 when it calculates the Berry Phase average described above.
(N.B. the \f$w_i\f$ in the above expression are all set equal to 1 in this case)

The above input is fine we can, however, refine this somewhat by making use of a multicolvar transform action as shown below:

\plumedfile
c1: COORDINATIONNUMBER SPECIES=1-512 SWITCH={EXP D_0=4.0 R_0=0.5}
cf: MTRANSFORM_MORE DATA=c1 SWITCH={RATIONAL D_0=2.0 R_0=0.1} LOWMEM
cc: CENTER_OF_MULTICOLVAR DATA=cf
\endplumedfile

This input once again calculates the coordination numbers of all the atoms in the system.  The middle line then transforms these
coordination numbers to numbers between 0 and 1.  Essentially any atom with a coordination number larger than 2.0 is given a weight
of one and below this value the transformed value decays to zero.  It is these transformed coordination numbers that are used to calculate
the Berry phase average described in the previous section.

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace multicolvar {


class CenterOfMultiColvar : public ActionWithVirtualAtom {
private:
  unsigned comp;
  vesselbase::StoreDataVessel* mystash;
  MultiColvarBase* mycolv;
public:
  static void registerKeywords( Keywords& keys );
  explicit CenterOfMultiColvar(const ActionOptions&ao);
  void calculate();
};

PLUMED_REGISTER_ACTION(CenterOfMultiColvar,"CENTER_OF_MULTICOLVAR")

void CenterOfMultiColvar::registerKeywords(Keywords& keys) {
  ActionWithVirtualAtom::registerKeywords(keys);
  keys.add("compulsory","DATA","find the average value for a multicolvar");
  keys.add("optional","COMPONENT","if your input multicolvar is a vector then specify which component you would like to use in calculating the weight");
}

CenterOfMultiColvar::CenterOfMultiColvar(const ActionOptions&ao):
  Action(ao),
  ActionWithVirtualAtom(ao)
{
  std::string mlab; parse("DATA",mlab);
  mycolv= plumed.getActionSet().selectWithLabel<MultiColvarBase*>(mlab);
  if(!mycolv) error("action labelled " +  mlab + " does not exist or does not have vessels");
  // Copy the atoms from the input multicolvar
  BridgedMultiColvarFunction* mybr=dynamic_cast<BridgedMultiColvarFunction*>( mycolv );
  if( mybr ) {
    requestAtoms( (mybr->getPntrToMultiColvar())->getAbsoluteIndexes() ); comp=1;
  } else {
    if( mycolv->getNumberOfQuantities()>5 ) {
      int incomp=-1; parse("COMPONENT",incomp);
      if( incomp<0 ) error("vector input but component was not specified");
      comp=incomp;
    } else {
      comp=1;
    }
    requestAtoms( mycolv->getAbsoluteIndexes () );
  }
  // We need the derivatives
  mycolv->turnOnDerivatives(); addDependency(mycolv);
  mystash = mycolv->buildDataStashes( NULL );
  log.printf("  building center of mass based on weights calculated in multicolvar action named %s \n",mycolv->getLabel().c_str() );
}

void CenterOfMultiColvar::calculate() {
  // Retrieve the periodic boundary conditions
  const Pbc& pbc=mycolv->getPbc();
  if( !pbc.isOrthorombic() ) error("Berry phase does not work for non orthorhombic cells");

  // Create a multivalue to store the derivatives
  MultiValue myvals( 7, mycolv->getNumberOfDerivatives() ); myvals.clearAll();
  MultiValue tvals( mycolv->getNumberOfQuantities(), mycolv->getNumberOfDerivatives() );
  tvals.clearAll();

  // Now loop over all active multicolvars
  Vector stmp, ctmp, scom, ccom, sder, cder;
  scom.zero(); ccom.zero(); double norm=0;
  std::vector<double> cvals( mycolv->getNumberOfQuantities() );
  for(unsigned i=0; i<mystash->getNumberOfStoredValues(); ++i) {
    // Retrieve value and derivatives
    mystash->retrieveSequentialValue( i, false, cvals );
    mystash->retrieveDerivatives( mycolv->getPositionInFullTaskList(i), false, tvals );
    // Convert position into fractionals
    Vector fpos = pbc.realToScaled( mycolv->getCentralAtomPos( mycolv->getPositionInFullTaskList(i) ) );
    // Now accumulate Berry phase averages
    for(unsigned j=0; j<3; ++j) {
      stmp[j] = sin( 2*pi*fpos[j] ); ctmp[j] = cos( 2*pi*fpos[j] );
      scom[j] += cvals[0]*cvals[comp]*stmp[j]; ccom[j] += cvals[0]*cvals[comp]*ctmp[j];
      double icell = 1.0 / getPbc().getBox().getRow(j).modulo();
      sder[j] = 2*pi*icell*cvals[0]*cvals[comp]*cos( 2*pi*fpos[j] );
      cder[j]=-2*pi*icell*cvals[0]*cvals[comp]*sin( 2*pi*fpos[j] );
    }
    // Now accumulate derivatives
    for(unsigned k=0; k<tvals.getNumberActive(); ++k) {
      unsigned icomp=tvals.getActiveIndex(k);
      myvals.addDerivative( 0, icomp, cvals[0]*tvals.getDerivative( comp, icomp ) + cvals[comp]*tvals.getDerivative( 0, icomp ) );
      for(unsigned k=0; k<3; ++k) {
        myvals.addDerivative( 1+k, icomp, stmp[k]*( cvals[0]*tvals.getDerivative( comp, icomp ) +
                              cvals[comp]*tvals.getDerivative( 0, icomp ) ) );
        myvals.addDerivative( 4+k, icomp, ctmp[k]*( cvals[0]*tvals.getDerivative( comp, icomp ) +
                              cvals[comp]*tvals.getDerivative( 0, icomp ) ) );
      }
    }
    // Get the central atom pack
    CatomPack mypack; mycolv->getCentralAtomPack( 0, mycolv->getPositionInFullTaskList(i), mypack );
    for(unsigned j=0; j<mypack.getNumberOfAtomsWithDerivatives(); ++j) {
      unsigned jder=3*mypack.getIndex(j);
      // Derivatives of sine
      myvals.addDerivative( 1, jder+0, mypack.getDerivative(j, 0, sder) );
      myvals.addDerivative( 2, jder+1, mypack.getDerivative(j, 1, sder) );
      myvals.addDerivative( 3, jder+2, mypack.getDerivative(j, 2, sder) );
      // Derivatives of cosine
      myvals.addDerivative( 4, jder+0, mypack.getDerivative(j, 0, cder) );
      myvals.addDerivative( 5, jder+1, mypack.getDerivative(j, 1, cder) );
      myvals.addDerivative( 6, jder+2, mypack.getDerivative(j, 2, cder) );
    }
    norm += cvals[0]*cvals[comp]; tvals.clearAll();
  }

  // And now finish Berry phase average
  scom /= norm; ccom /=norm; Vector cpos;
  for(unsigned j=0; j<3; ++j) cpos[j] = atan2( scom[j], ccom[j] ) / (2*pi);
  Vector cart_pos = pbc.scaledToReal( cpos );
  setPosition(cart_pos); setMass(1.0);   // This could be much cleverer but not without changing many things in PLMED

  // And derivatives
  Vector tander; myvals.updateDynamicList(); double inv_weight = 1.0 / norm;
  for(unsigned j=0; j<3; ++j) {
    double tmp = scom[j] / ccom[j];
    tander[j] = getPbc().getBox().getRow(j).modulo() / (2*pi*( 1 + tmp*tmp ));
  }
  for(unsigned i=0; i<myvals.getNumberActive(); ++i) {
    unsigned ider=myvals.getActiveIndex(i);
    for(unsigned j=0; j<3; ++j) {
      double sderv = inv_weight*myvals.getDerivative(1+j,ider) - inv_weight*scom[j]*myvals.getDerivative(0,ider);
      double cderv = inv_weight*myvals.getDerivative(4+j,ider) - inv_weight*ccom[j]*myvals.getDerivative(0,ider);
      myvals.setDerivative( 1+j, ider, tander[j]*(sderv/ccom[j]  - scom[j]*cderv/(ccom[j]*ccom[j])) );
      //if( j==2 ) printf("DERIV %d %10.4f %10.4f %10.4f %10.4f \n",i,myvals.getDerivative(0,ider),sderv,cderv,myvals.getDerivative(1+j,ider ) );
    }
  }

  // Atom derivatives
  std::vector<Tensor> fderiv( getNumberOfAtoms() );
  for(unsigned j=0; j<getNumberOfAtoms(); ++j) {
    for(unsigned k=0; k<3; ++k) {
      if( myvals.isActive(3*j+k) ) for(unsigned n=0; n<3; ++n) fderiv[j](k,n) = myvals.getDerivative( 1+n, 3*j+k );
      else for(unsigned n=0; n<3; ++n) fderiv[j](k,n) = 0;
    }
  }
  setAtomsDerivatives( fderiv );
  // Box derivatives?
}

}
}
