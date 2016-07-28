/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2016 The plumed team
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
#include "multicolvar/MultiColvarFunction.h"
#include "core/ActionRegister.h"
#include "tools/SwitchingFunction.h"
#include "tools/Torsion.h" 
    
#include <string>
#include <cmath>

//+PLUMEDOC MCOLVARF INTERMOLECULARTORSIONS
/*
Calculate torsions between vectors on adjacent molecules

This variable can be used to calculate the average torsional angles between vectors.  In other words,
it can be used to compute quantities like this:

\f[
s = \frac{ \sum_{i \ne j} \sigma(r_{ij}) \theta_{ij} }{ \sum_{i \ne j} \sigma(r_{ij}) }
f\]

Here the sums run over all pairs of molecules. \f$\sigma(r_{ij})\f$ is a \ref switchingfunction that
action on the distance between the centers of molecules \f$i\f$ and \f$j\f$.  \f$\theta_{ij}\f$ is then
the torsional angle between an orientation vector for molecule \f$i\f$ and molecule \f$j\f$.

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace crystallization {

class InterMolecularTorsions : public multicolvar::MultiColvarFunction {
private:
/// The switching function that tells us if atoms are close enough together
  SwitchingFunction switchingFunction;
public:
  static void registerKeywords( Keywords& keys );
  explicit InterMolecularTorsions(const ActionOptions&);
/// Do the stuff with the switching functions
  double calculateWeight( const unsigned& taskCode, const double& weight, multicolvar::AtomValuePack& myatoms ) const ;
/// Actually do the calculation
  double compute( const unsigned& tindex, multicolvar::AtomValuePack& myatoms ) const ;
/// Is the variable periodic
  bool isPeriodic(){ return true; }
  void retrieveDomain( std::string& min, std::string& max ){ min="-pi"; max="+pi"; }
};

PLUMED_REGISTER_ACTION(InterMolecularTorsions,"INTERMOLECULARTORSIONS")

void InterMolecularTorsions::registerKeywords( Keywords& keys ){
  MultiColvarFunction::registerKeywords( keys );
  keys.add("atoms","MOLS","The molecules you would like to calculate the torsional angles between. This should be the label/s of \\ref MOLECULES or \\ref PLANES actions");
  keys.add("atoms-1","MOLSA","In this version of the input the torsional angles between all pairs of atoms including one atom from MOLA one atom from MOLB will be computed. " 
                             "This should be the label/s of \\ref MOLECULES or \\ref PLANES actions");
  keys.add("atoms-1","MOLSB","In this version of the input the torsional angles between all pairs of atoms including one atom from MOLA one atom from MOLB will be computed. "
                             "This should be the label/s of \\ref MOLECULES or \\ref PLANES actions");  
  keys.add("compulsory","NN","6","The n parameter of the switching function ");
  keys.add("compulsory","MM","0","The m parameter of the switching function; 0 implies 2*NN");
  keys.add("compulsory","D_0","0.0","The d_0 parameter of the switching function");
  keys.add("compulsory","R_0","The r_0 parameter of the switching function");
  keys.add("optional","SWITCH","This keyword is used if you want to employ an alternative to the continuous swiching function defined above. "
                               "The following provides information on the \\ref switchingfunction that are available. "
                               "When this keyword is present you no longer need the NN, MM, D_0 and R_0 keywords.");
  // Use actionWithDistributionKeywords
  keys.remove("LOWMEM"); keys.remove("DATA"); 
  keys.addFlag("LOWMEM",false,"lower the memory requirements");
}

InterMolecularTorsions::InterMolecularTorsions(const ActionOptions& ao):
Action(ao),
MultiColvarFunction(ao)
{
  for(unsigned i=0;i<getNumberOfBaseMultiColvars();++i){
      if( getBaseMultiColvar(i)->getNumberOfQuantities()!=5 ) error("input multicolvar does not calculate molecular orientations");
  }
  // The weight of this does have derivatives
  weightHasDerivatives=true;

  // Read in the switching function
  std::string sw, errors; parse("SWITCH",sw);
  if(sw.length()>0){
     switchingFunction.set(sw,errors);
  } else {
     double r_0=-1.0, d_0; int nn, mm;
     parse("NN",nn); parse("MM",mm);
     parse("R_0",r_0); parse("D_0",d_0);
     if( r_0<0.0 ) error("you must set a value for R_0");
     switchingFunction.set(nn,mm,r_0,d_0);
  }
  log.printf("  calculating number of links with atoms separation of %s\n",( switchingFunction.description() ).c_str() );
  std::vector<AtomNumber> all_atoms; readTwoGroups( "MOLS", "MOLSA", "MOLSB", all_atoms );
  setupMultiColvarBase( all_atoms ); setLinkCellCutoff( switchingFunction.get_dmax() );

  for(unsigned i=0;i<getNumberOfBaseMultiColvars();++i){
    if( !getBaseMultiColvar(i)->hasDifferentiableOrientation() ) error("cannot use multicolvar of type " + getBaseMultiColvar(i)->getName() );
  }

  // Create holders for the collective variable
  readVesselKeywords();
  plumed_assert( getNumberOfVessels()==0 );
  std::string input; addVessel( "SUM", input, -1 );
  readVesselKeywords();
}

double InterMolecularTorsions::calculateWeight( const unsigned& taskCode, const double& weight, multicolvar::AtomValuePack& myatoms ) const {
  Vector distance = getSeparation( myatoms.getPosition(0), myatoms.getPosition(1) );
  double dfunc, sw = switchingFunction.calculateSqr( distance.modulo2(), dfunc );

  if( !doNotCalculateDerivatives() ){
      addAtomDerivatives( 0, 0, (-dfunc)*weight*distance, myatoms );
      addAtomDerivatives( 0, 1, (dfunc)*weight*distance, myatoms );
      myatoms.addBoxDerivatives( 0, (-dfunc)*weight*Tensor(distance,distance) );
  }
  return sw;
}

double InterMolecularTorsions::compute( const unsigned& tindex, multicolvar::AtomValuePack& myatoms ) const {
   Vector v1, v2, dv1, dv2, dconn, conn = getSeparation( myatoms.getPosition(0), myatoms.getPosition(1) );

   // Retrieve vectors
   std::vector<double> orient0( 5 ), orient1( 5 ); 
   getInputData( 0, true, myatoms, orient0 );
   getInputData( 1, true, myatoms, orient1 );
   for(unsigned i=0;i<3;++i){ v1[i]=orient0[2+i]; v2[i]=orient1[2+i]; }
   if( getBaseMultiColvar(0)->getNumberOfQuantities()<3 ) return 1.0;

   // Evaluate angle
   Torsion t; double angle = t.compute( v1, conn, v2, dv1, dconn, dv2 );
   for(unsigned i=0;i<3;++i){ orient0[i+2]=dv1[i]; orient1[i+2]=dv2[i]; }

   // And accumulate derivatives
   if( !doNotCalculateDerivatives() ){
      MultiValue& myder0=getInputDerivatives( 0, true, myatoms );
      mergeInputDerivatives( 1, 2, orient1.size(), 0, orient0, myder0, myatoms );
      MultiValue& myder1=getInputDerivatives( 1, true, myatoms );
      mergeInputDerivatives( 1, 2, orient0.size(), 1, orient1, myder1, myatoms );
      addAtomDerivatives( 1, 0, -dconn, myatoms ); addAtomDerivatives( 1, 1, dconn, myatoms ); 
      myatoms.addBoxDerivatives( 1, -extProduct( conn, dconn ) );
   }

   return angle;
}

}
}
