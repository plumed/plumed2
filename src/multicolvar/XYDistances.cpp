/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2016 The plumed team
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
#include "MultiColvar.h"
#include "core/ActionRegister.h"

#include <string>
#include <cmath>

using namespace std;

namespace PLMD{
namespace multicolvar{

//+PLUMEDOC MCOLVAR XYDISTANCES
/*
Calculate distance between a pair of atoms neglecting the z-component.  
You can then calculate functions of the distribution of
values such as the minimum, the number less than a certain quantity and so on. 

\par Examples

The following input tells plumed to calculate the projection of the length of the vector connecting atom 3 
to atom 5 and that is projected in the xy-plane and the projection of the length of the vector 
the vector connecting atom 1 to atom 2 in the xy-plane.  The minimum of these two quantities is then 
printed 
\verbatim
XYDISTANCES ATOMS1=3,5 ATOMS2=1,2 MIN={BETA=0.1} LABEL=d1
PRINT ARG=d1.min
\endverbatim
(See also \ref PRINT).

*/
//+ENDPLUMEDOC

//+PLUMEDOC MCOLVAR XZDISTANCES
/*
Calculate distance between a pair of atoms neglecting the y-component.  
You can then calculate functions of the distribution of
values such as the minimum, the number less than a certain quantity and so on. 

\par Examples

See documentation for \ref XYDISTANCES for examples of how to use this command.
You just need to substitute XZDISTANCES for XYDISTANCES to investigate the xz component
rather than the xy component.

*/
//+ENDPLUMEDOC

//+PLUMEDOC MCOLVAR YZDISTANCES
/*
Calculate distance between a pair of atoms neglecting the x-component.  
You can then calculate functions of the distribution of
values such as the minimum, the number less than a certain quantity and so on. 

\par Examples

See documentation for \ref XYDISTANCES for examples of how to use this command.
You just need to substitute YZDISTANCES for XYDISTANCES to investigate the yz component
rather than the xy component.

*/
//+ENDPLUMEDOC


class XYDistances : public MultiColvar {
private:
  unsigned myc1, myc2;
public:
  static void registerKeywords( Keywords& keys );
  explicit XYDistances(const ActionOptions&);
// active methods:
  virtual double compute( const unsigned& tindex, AtomValuePack& myatoms ) const ;
/// Returns the number of coordinates of the field
  bool isPeriodic(){ return false; }
};

PLUMED_REGISTER_ACTION(XYDistances,"XYDISTANCES")
PLUMED_REGISTER_ACTION(XYDistances,"XZDISTANCES")
PLUMED_REGISTER_ACTION(XYDistances,"YZDISTANCES")

void XYDistances::registerKeywords( Keywords& keys ){
  MultiColvar::registerKeywords( keys );
  keys.use("ATOMS"); keys.use("MAX"); keys.use("ALT_MIN"); 
  keys.use("MEAN"); keys.use("MIN"); keys.use("LESS_THAN");
  keys.use("LOWEST"); keys.use("HIGHEST"); 
  keys.use("MORE_THAN"); keys.use("BETWEEN"); keys.use("HISTOGRAM"); keys.use("MOMENTS");
  keys.add("atoms-1","GROUP","Calculate the distance between each distinct pair of atoms in the group");
  keys.add("atoms-2","GROUPA","Calculate the distances between all the atoms in GROUPA and all "
                              "the atoms in GROUPB. This must be used in conjuction with GROUPB.");
  keys.add("atoms-2","GROUPB","Calculate the distances between all the atoms in GROUPA and all the atoms "
                              "in GROUPB. This must be used in conjuction with GROUPA.");
}

XYDistances::XYDistances(const ActionOptions&ao):
PLUMED_MULTICOLVAR_INIT(ao)
{
  if( getName().find("XY")!=std::string::npos){
      myc1=0; myc2=1;
  } else if( getName().find("XZ")!=std::string::npos){
      myc1=0; myc2=2;
  } else if( getName().find("YZ")!=std::string::npos){
      myc1=1; myc2=2;
  } else plumed_error();

  // Read in the atoms
  std::vector<AtomNumber> all_atoms;
  readTwoGroups( "GROUP", "GROUPA", "GROUPB", all_atoms );
  int natoms=2; readAtoms( natoms, all_atoms );
  // And check everything has been read in correctly
  checkRead();
}

double XYDistances::compute( const unsigned& tindex, AtomValuePack& myatoms ) const {
   Vector distance; 
   distance=getSeparation( myatoms.getPosition(0), myatoms.getPosition(1) );
   const double value=sqrt(distance[myc1]*distance[myc1] + distance[myc2]*distance[myc2] );
   const double invvalue=1.0/value;

   Vector myvec; myvec.zero(); 
   // And finish the calculation
   myvec[myc1]=+invvalue*distance[myc1]; myvec[myc2]=+invvalue*distance[myc2]; addAtomDerivatives( 1, 1, myvec, myatoms  );
   myvec[myc1]=-invvalue*distance[myc1]; myvec[myc2]=-invvalue*distance[myc2]; addAtomDerivatives( 1, 0, myvec, myatoms );
   myatoms.addBoxDerivatives( 1, Tensor(distance,myvec) );
   return value;
}

}
}

