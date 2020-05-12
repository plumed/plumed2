/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2020 The plumed team
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
#include "OrientationSphere.h"
#include "core/PlumedMain.h"
#include "VectorMultiColvar.h"

using namespace std;

namespace PLMD {
namespace crystallization {

void OrientationSphere::registerKeywords( Keywords& keys ) {
  multicolvar::MultiColvarBase::registerKeywords( keys );
  keys.add("compulsory","NN","6","The n parameter of the switching function ");
  keys.add("compulsory","MM","0","The m parameter of the switching function; 0 implies 2*NN");
  keys.add("compulsory","D_0","0.0","The d_0 parameter of the switching function");
  keys.add("compulsory","R_0","The r_0 parameter of the switching function");
  keys.add("optional","SWITCH","This keyword is used if you want to employ an alternative to the continuous switching function defined above. "
           "The following provides information on the \\ref switchingfunction that are available. "
           "When this keyword is present you no longer need the NN, MM, D_0 and R_0 keywords.");
  // Use actionWithDistributionKeywords
  keys.use("SPECIES"); keys.use("SPECIESA"); keys.use("SPECIESB");
  keys.use("MEAN"); keys.use("MORE_THAN"); keys.use("LESS_THAN");
  keys.use("MIN"); keys.use("BETWEEN"); keys.use("HISTOGRAM"); keys.use("MOMENTS");
  keys.use("LOWEST"); keys.use("HIGHEST");
}

OrientationSphere::OrientationSphere(const ActionOptions&ao):
  Action(ao),
  MultiColvarBase(ao)
{
  if( getNumberOfBaseMultiColvars()>1 ) warning("not sure if orientation sphere works with more than one base multicolvar - check numerical derivatives");
  // Read in the switching function
  std::string sw, errors; parse("SWITCH",sw);
  if(sw.length()>0) {
    switchingFunction.set(sw,errors);
  } else {
    double r_0=-1.0, d_0; int nn, mm;
    parse("NN",nn); parse("MM",mm);
    parse("R_0",r_0); parse("D_0",d_0);
    if( r_0<0.0 ) error("you must set a value for R_0");
    switchingFunction.set(nn,mm,r_0,d_0);
  }
  log.printf("  degree of overlap in orientation between central molecule and those within %s\n",( switchingFunction.description() ).c_str() );
  log<<"  Bibliography "<<plumed.cite("Tribello, Giberti, Sosso, Salvalaglio and Parrinello, J. Chem. Theory Comput. 13, 1317 (2017)")<<"\n";
  // Set the link cell cutoff
  rcut2 = switchingFunction.get_dmax()*switchingFunction.get_dmax();
  setLinkCellCutoff( switchingFunction.get_dmax() );
  std::vector<AtomNumber> all_atoms; setupMultiColvarBase( all_atoms );
}

double OrientationSphere::compute( const unsigned& tindex, multicolvar::AtomValuePack& myatoms ) const {
  double sw, value=0, denom=0, dfunc; Vector ddistance;
  unsigned ncomponents=getBaseMultiColvar(0)->getNumberOfQuantities();
  std::vector<double> catom_orient( ncomponents ), this_orient( ncomponents );
  std::vector<double> this_der( ncomponents ), catom_der( ncomponents );

  getInputData( 0, true, myatoms, catom_orient );
  MultiValue& myder0=getInputDerivatives( 0, true, myatoms );

  for(unsigned i=1; i<myatoms.getNumberOfAtoms(); ++i) {
    Vector& distance=myatoms.getPosition(i);
    double d2;
    if ( (d2=distance[0]*distance[0])<rcut2 &&
         (d2+=distance[1]*distance[1])<rcut2 &&
         (d2+=distance[2]*distance[2])<rcut2 &&
         d2>epsilon ) {

      sw = switchingFunction.calculateSqr( d2, dfunc );

      getInputData( i, true, myatoms, this_orient );
      // Calculate the dot product wrt to this position
      double f_dot = computeVectorFunction( distance, catom_orient, this_orient, ddistance, catom_der, this_der );

      if( !doNotCalculateDerivatives() ) {
        for(unsigned k=2; k<catom_orient.size(); ++k) { this_der[k]*=sw; catom_der[k]*=sw; }
        MultiValue& myder1=getInputDerivatives( i, true, myatoms );
        mergeInputDerivatives( 1, 2, this_orient.size(), 0, catom_der, myder0, myatoms );
        mergeInputDerivatives( 1, 2, catom_der.size(), i, this_der, myder1, myatoms );
        addAtomDerivatives( 1, 0, f_dot*(-dfunc)*distance - sw*ddistance, myatoms );
        addAtomDerivatives( 1, i, f_dot*(dfunc)*distance + sw*ddistance, myatoms );
        myatoms.addBoxDerivatives( 1, (-dfunc)*f_dot*Tensor(distance,distance) - sw*extProduct(distance,ddistance) );
        myder1.clearAll();

        addAtomDerivatives( -1, 0, (-dfunc)*distance, myatoms );
        addAtomDerivatives( -1, i, (dfunc)*distance, myatoms );
        myatoms.addTemporyBoxDerivatives( (-dfunc)*Tensor(distance,distance) );

      }
      value += sw*f_dot;
      denom += sw;
    }
  }
  double rdenom, df2, pref=calculateCoordinationPrefactor( denom, df2 );
  if( fabs(denom)>epsilon ) { rdenom = 1.0 / denom; }
  else { plumed_assert(fabs(value)<epsilon); rdenom=1.0; }

  // Now divide everything
  double rdenom2=rdenom*rdenom;
  updateActiveAtoms( myatoms ); MultiValue& myvals=myatoms.getUnderlyingMultiValue();
  for(unsigned i=0; i<myvals.getNumberActive(); ++i) {
    unsigned ider=myvals.getActiveIndex(i);
    double  dgd=myvals.getTemporyDerivative(ider);
    myvals.setDerivative( 1, ider, rdenom*(pref*myvals.getDerivative(1,ider)+value*df2*dgd) - (value*pref*dgd)*rdenom2 );
  }

  return pref*rdenom*value;
}

}
}

