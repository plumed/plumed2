/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2019 The plumed team
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
#include "MultiColvarBase.h"
#include "AtomValuePack.h"
#include "tools/Angle.h"
#include "tools/SwitchingFunction.h"
#include "core/ActionRegister.h"

#include <string>
#include <cmath>

using namespace std;

namespace PLMD {
namespace multicolvar {

//+PLUMEDOC MCOLVAR ANGLES
/*
Calculate functions of the distribution of angles .

You can use this command to calculate functions such as:

\f[
 f(x) = \sum_{ijk} g( \theta_{ijk} )
\f]

Alternatively you can use this command to calculate functions such as:

\f[
f(x) = \sum_{ijk} s(r_{ij})s(r_{jk}) g(\theta_{ijk})
\f]

where \f$s(r)\f$ is a \ref switchingfunction.  This second form means that you can
use this to calculate functions of the angles in the first coordination sphere of
an atom / molecule \cite lj-recon.

\par Examples

The following example instructs plumed to find the average of two angles and to
print it to a file

\plumedfile
ANGLES ATOMS1=1,2,3 ATOMS2=4,5,6 MEAN LABEL=a1
PRINT ARG=a1.mean FILE=colvar
\endplumedfile

The following example tells plumed to calculate all angles involving
at least one atom from GROUPA and two atoms from GROUPB in which the distances
are less than 1.0. The number of angles between \f$\frac{\pi}{4}\f$ and
\f$\frac{3\pi}{4}\f$ is then output

\plumedfile
ANGLES GROUPA=1-10 GROUPB=11-100 BETWEEN={GAUSSIAN LOWER=0.25pi UPPER=0.75pi} SWITCH={GAUSSIAN R_0=1.0} LABEL=a1
PRINT ARG=a1.between FILE=colvar
\endplumedfile

This final example instructs plumed to calculate all the angles in the first coordination
spheres of the atoms. The bins for a normalized histogram of the distribution is then output

\plumedfile
ANGLES GROUP=1-38 HISTOGRAM={GAUSSIAN LOWER=0.0 UPPER=pi NBINS=20} SWITCH={GAUSSIAN R_0=1.0} LABEL=a1
PRINT ARG=a1.* FILE=colvar
\endplumedfile

*/
//+ENDPLUMEDOC

class Angles : public MultiColvarBase {
private:
  bool use_sf;
  double rcut2_1, rcut2_2;
  SwitchingFunction sf1;
  SwitchingFunction sf2;
public:
  static void registerKeywords( Keywords& keys );
  explicit Angles(const ActionOptions&);
/// Updates neighbor list
  virtual double compute( const unsigned& tindex, AtomValuePack& ) const ;
/// Returns the number of coordinates of the field
  double calculateWeight( const unsigned& taskCode, const double& weight, AtomValuePack& ) const ;
  bool isPeriodic() { return false; }
};

PLUMED_REGISTER_ACTION(Angles,"ANGLES")

void Angles::registerKeywords( Keywords& keys ) {
  MultiColvarBase::registerKeywords( keys );
  keys.use("MEAN"); keys.use("LESS_THAN");
  keys.use("BETWEEN"); keys.use("HISTOGRAM"); keys.use("MORE_THAN");
  // Could also add Region here in theory
  keys.add("numbered","ATOMS","the atoms involved in each of the angles you wish to calculate. "
           "Keywords like ATOMS1, ATOMS2, ATOMS3,... should be listed and one angle will be "
           "calculated for each ATOM keyword you specify (all ATOM keywords should "
           "provide the indices of three atoms).  The eventual number of quantities calculated by this "
           "action will depend on what functions of the distribution you choose to calculate.");
  keys.reset_style("ATOMS","atoms");
  keys.add("atoms-1","GROUP","Calculate angles for each distinct set of three atoms in the group");
  keys.add("atoms-2","GROUPA","A group of central atoms about which angles should be calculated");
  keys.add("atoms-2","GROUPB","When used in conjunction with GROUPA this keyword instructs plumed "
           "to calculate all distinct angles involving one atom from GROUPA "
           "and two atoms from GROUPB. The atom from GROUPA is the central atom.");
  keys.add("atoms-3","GROUPC","This must be used in conjunction with GROUPA and GROUPB.  All angles "
           "involving one atom from GROUPA, one atom from GROUPB and one atom from "
           "GROUPC are calculated. The GROUPA atoms are assumed to be the central "
           "atoms");
  keys.add("optional","SWITCH","A switching function that ensures that only angles between atoms that "
           "are within a certain fixed cutoff are calculated. The following provides "
           "information on the \\ref switchingfunction that are available.");
  keys.add("optional","SWITCHA","A switching function on the distance between the atoms in group A and the atoms in "
           "group B");
  keys.add("optional","SWITCHB","A switching function on the distance between the atoms in group A and the atoms in "
           "group B");
}

Angles::Angles(const ActionOptions&ao):
  Action(ao),
  MultiColvarBase(ao),
  use_sf(false)
{
  std::string sfinput,errors; parse("SWITCH",sfinput);
  if( sfinput.length()>0 ) {
    use_sf=true;
    weightHasDerivatives=true;
    sf1.set(sfinput,errors);
    if( errors.length()!=0 ) error("problem reading SWITCH keyword : " + errors );
    sf2.set(sfinput,errors);
    if( errors.length()!=0 ) error("problem reading SWITCH keyword : " + errors );
    log.printf("  only calculating angles for atoms separated by less than %s\n", sf1.description().c_str() );
  } else {
    parse("SWITCHA",sfinput);
    if(sfinput.length()>0) {
      use_sf=true;
      weightHasDerivatives=true;
      sf1.set(sfinput,errors);
      if( errors.length()!=0 ) error("problem reading SWITCHA keyword : " + errors );
      sfinput.clear(); parse("SWITCHB",sfinput);
      if(sfinput.length()==0) error("found SWITCHA keyword without SWITCHB");
      sf2.set(sfinput,errors);
      if( errors.length()!=0 ) error("problem reading SWITCHB keyword : " + errors );
      log.printf("  only calculating angles when the distance between GROUPA and GROUPB atoms is less than %s\n", sf1.description().c_str() );
      log.printf("  only calculating angles when the distance between GROUPA and GROUPC atoms is less than %s\n", sf2.description().c_str() );
    }
  }
  // Read in the atoms
  std::vector<AtomNumber> all_atoms;
  readGroupKeywords( "GROUP", "GROUPA", "GROUPB", "GROUPC", false, true, all_atoms );
  if( atom_lab.size()==0 ) readAtomsLikeKeyword( "ATOMS", 3, all_atoms );
  setupMultiColvarBase( all_atoms );
  // Set cutoff for link cells
  if( use_sf ) {
    setLinkCellCutoff( sf1.get_dmax() );
    rcut2_1=sf1.get_dmax()*sf1.get_dmax();
    rcut2_2=sf2.get_dmax()*sf2.get_dmax();
  }

  // And check everything has been read in correctly
  checkRead();
  // Setup stuff for central atom
  std::vector<bool> catom_ind(3, false); catom_ind[0]=true;
  setAtomsForCentralAtom( catom_ind );
}

double Angles::calculateWeight( const unsigned& taskCode, const double& weight, AtomValuePack& myatoms ) const {
  if(!use_sf) return 1.0;
  Vector dij=getSeparation( myatoms.getPosition(0), myatoms.getPosition(2) );
  Vector dik=getSeparation( myatoms.getPosition(0), myatoms.getPosition(1) );

  double w1, w2, dw1, dw2, wtot;
  double ldij = dij.modulo2(), ldik = dik.modulo2();

  if( use_sf ) {
    if( ldij>rcut2_1 || ldik>rcut2_2 ) return 0.0;
  }

  w1=sf1.calculateSqr( ldij, dw1 );
  w2=sf2.calculateSqr( ldik, dw2 );
  wtot=w1*w2; dw1*=weight*w2; dw2*=weight*w1;

  addAtomDerivatives( 0, 1, dw2*dik, myatoms );
  addAtomDerivatives( 0, 0, -dw1*dij - dw2*dik, myatoms );
  addAtomDerivatives( 0, 2, dw1*dij, myatoms );
  myatoms.addBoxDerivatives( 0, (-dw1)*Tensor(dij,dij) + (-dw2)*Tensor(dik,dik) );
  return wtot;
}

double Angles::compute( const unsigned& tindex, AtomValuePack& myatoms ) const {
  Vector dij=getSeparation( myatoms.getPosition(0), myatoms.getPosition(2) );
  Vector dik=getSeparation( myatoms.getPosition(0), myatoms.getPosition(1) );

  Vector ddij,ddik; PLMD::Angle a;
  double angle=a.compute(dij,dik,ddij,ddik);

  // And finish the calculation
  addAtomDerivatives( 1, 1, ddik, myatoms );
  addAtomDerivatives( 1, 0, - ddik - ddij, myatoms );
  addAtomDerivatives( 1, 2, ddij, myatoms );
  myatoms.addBoxDerivatives( 1, -(Tensor(dij,ddij)+Tensor(dik,ddik)) );

  return angle;
}

}
}
