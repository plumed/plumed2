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
#include "MultiColvarBase.h"
#include "AtomValuePack.h"
#include "core/ActionRegister.h"
#include "tools/Angle.h"
#include "tools/SwitchingFunction.h"

#include <string>
#include <cmath>

using namespace std;

namespace PLMD {
namespace multicolvar {

//+PLUMEDOC MCOLVAR XANGLES
/*
Calculate the angles between the vector connecting two atoms and the x axis.

\par Examples

The following input tells plumed to calculate the angles between the x-axis and the vector connecting atom 3 to atom 5 and between the x-axis
and the vector connecting atom 1 to atom 2.  The minimum of these two quantities is then
\plumedfile
XANGLES ATOMS1=3,5 ATOMS2=1,2 MIN={BETA=0.1} LABEL=d1
PRINT ARG=d1.min
\endplumedfile
(See also \ref PRINT).
*/
//+ENDPLUMEDOC

//+PLUMEDOC MCOLVAR YANGLES
/*
Calculate the angles between the vector connecting two atoms and the y axis.

\par Examples

The following input tells plumed to calculate the angles between the y-axis and the vector connecting atom 3 to atom 5 and between the y-axis
and the vector connecting atom 1 to atom 2.  The minimum of these two quantities is then
\plumedfile
YANGLES ATOMS1=3,5 ATOMS2=1,2 MIN={BETA=0.1} LABEL=d1
PRINT ARG=d1.min
\endplumedfile
(See also \ref PRINT).
*/
//+ENDPLUMEDOC

//+PLUMEDOC MCOLVAR ZANGLES
/*
Calculate the angles between the vector connecting two atoms and the z axis.

\par Examples

The following input tells plumed to calculate the angles between the z-axis and the vector connecting atom 3 to atom 5 and between the z-axis
and the vector connecting atom 1 to atom 2.  The minimum of these two quantities is then
\plumedfile
ZANGLES ATOMS1=3,5 ATOMS2=1,2 MIN={BETA=0.1} LABEL=d1
PRINT ARG=d1.min
\endplumedfile
(See also \ref PRINT).
*/
//+ENDPLUMEDOC



class XAngles : public MultiColvarBase {
private:
  bool use_sf;
  unsigned myc;
  SwitchingFunction sf1;
public:
  static void registerKeywords( Keywords& keys );
  explicit XAngles(const ActionOptions&);
// active methods:
  double compute( const unsigned& tindex, AtomValuePack& myatoms ) const override;
  double calculateWeight( const unsigned& taskCode, const double& weight, AtomValuePack& ) const override;
/// Returns the number of coordinates of the field
  bool isPeriodic() override { return false; }
};

PLUMED_REGISTER_ACTION(XAngles,"XANGLES")
PLUMED_REGISTER_ACTION(XAngles,"YANGLES")
PLUMED_REGISTER_ACTION(XAngles,"ZANGLES")

void XAngles::registerKeywords( Keywords& keys ) {
  MultiColvarBase::registerKeywords( keys );
  keys.use("MAX"); keys.use("ALT_MIN");
  keys.use("MEAN"); keys.use("MIN"); keys.use("LESS_THAN");
  keys.use("LOWEST"); keys.use("HIGHEST");
  keys.use("MORE_THAN"); keys.use("BETWEEN"); keys.use("HISTOGRAM"); keys.use("MOMENTS");
  keys.add("numbered","ATOMS","the atoms involved in each of the angles you wish to calculate. "
           "Keywords like ATOMS1, ATOMS2, ATOMS3,... should be listed and one angle will be "
           "calculated for each ATOM keyword you specify (all ATOM keywords should "
           "specify the indices of two atoms).  The eventual number of quantities calculated by this "
           "action will depend on what functions of the distribution you choose to calculate.");
  keys.reset_style("ATOMS","atoms");
  keys.add("atoms-1","GROUP","Calculate the distance between each distinct pair of atoms in the group");
  keys.add("atoms-2","GROUPA","Calculate the distances between all the atoms in GROUPA and all "
           "the atoms in GROUPB. This must be used in conjunction with GROUPB.");
  keys.add("atoms-2","GROUPB","Calculate the distances between all the atoms in GROUPA and all the atoms "
           "in GROUPB. This must be used in conjunction with GROUPA.");
  keys.add("optional","SWITCH","A switching function that ensures that only angles are only computed when atoms are within "
           "are within a certain fixed cutoff. The following provides information on the \\ref switchingfunction that are available.");
}

XAngles::XAngles(const ActionOptions&ao):
  Action(ao),
  MultiColvarBase(ao),
  use_sf(false)
{
  if( getName().find("X")!=std::string::npos) myc=0;
  else if( getName().find("Y")!=std::string::npos) myc=1;
  else if( getName().find("Z")!=std::string::npos) myc=2;
  else plumed_error();

  // Read in switching function
  std::string sfinput, errors; parse("SWITCH",sfinput);
  if( sfinput.length()>0 ) {
    use_sf=true; weightHasDerivatives=true;
    sf1.set(sfinput,errors);
    if( errors.length()!=0 ) error("problem reading SWITCH keyword : " + errors );
    log.printf("  only calculating angles for atoms separated by less than %s\n", sf1.description().c_str() );
    setLinkCellCutoff( sf1.get_dmax() );
  }

  // Read in the atoms
  std::vector<AtomNumber> all_atoms;
  readTwoGroups( "GROUP", "GROUPA", "GROUPB", all_atoms );
  if( atom_lab.size()==0 ) readAtomsLikeKeyword( "ATOMS", 2, all_atoms );
  setupMultiColvarBase( all_atoms );
  // And check everything has been read in correctly
  checkRead();
}

double XAngles::calculateWeight( const unsigned& taskCode, const double& weight, AtomValuePack& myatoms ) const {
  if(!use_sf) return 1.0;

  Vector distance=getSeparation( myatoms.getPosition(0), myatoms.getPosition(1) );
  double dw, w = sf1.calculateSqr( distance.modulo2(), dw );
  addAtomDerivatives( 0, 0, (-dw)*distance, myatoms );
  addAtomDerivatives( 0, 1, (+dw)*distance, myatoms );
  myatoms.addBoxDerivatives( 0, (-dw)*Tensor(distance,distance) );
  return w;
}

double XAngles::compute( const unsigned& tindex, AtomValuePack& myatoms ) const {
  Vector ddij, ddik, axis, distance; axis.zero(); axis[myc]=1;
  distance=getSeparation( myatoms.getPosition(0), myatoms.getPosition(1) );
  PLMD::Angle a; double angle=a.compute( distance, axis, ddij, ddik );

  addAtomDerivatives( 1, 0, -ddij, myatoms );
  addAtomDerivatives( 1, 1, ddij, myatoms );
  myatoms.addBoxDerivatives( 1, -Tensor( distance,ddij ) );
  return angle;
}

}
}

