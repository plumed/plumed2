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

//+PLUMEDOC MCOLVAR XDISTANCES
/*
Calculate the x components of the vectors connecting one or many pairs of atoms.  
You can then calculate functions of the distribution of
values such as the minimum, the number less than a certain quantity and so on. 

\par Examples

The following input tells plumed to calculate the x-component of the vector connecting atom 3 to atom 5 and 
the x-component of the vector connecting atom 1 to atom 2.  The minimum of these two quantities is then 
printed 
\verbatim
XDISTANCES ATOMS1=3,5 ATOMS2=1,2 MIN={BETA=0.1} LABEL=d1
PRINT ARG=d1.min
\endverbatim
(See also \ref PRINT).


The following input tells plumed to calculate the x-component of the vector connecting atom 3 to atom 5 and 
the x-component of the vector connecting atom 1 to atom 2.  The number of values that are 
less than 0.1nm is then printed to a file.
\verbatim
XDISTANCES ATOMS1=3,5 ATOMS2=1,2 LABEL=d1 LESS_THAN={RATIONAL R_0=0.1}
PRINT ARG=d1.lt0.1
\endverbatim
(See also \ref PRINT \ref switchingfunction).

The following input tells plumed to calculate the x-components of all the distinct vectors that can be created 
between atoms 1, 2 and 3 (i.e. the vectors between atoms 1 and 2, atoms 1 and 3 and atoms 2 and 3).  
The average of these quantities is then calculated.
\verbatim
XDISTANCES GROUP=1-3 AVERAGE LABEL=d1
PRINT ARG=d1.average
\endverbatim
(See also \ref PRINT)

The following input tells plumed to calculate all the vectors connecting the the atoms in GROUPA to the atoms in GROUPB.
In other words the vector between atoms 1 and 2 and the vector between atoms 1 and 3.  The number of values 
more than 0.1 is then printed to a file.
\verbatim
XDISTANCES GROUPA=1 GROUPB=2,3 MORE_THAN={RATIONAL R_0=0.1}
PRINT ARG=d1.gt0.1 
\endverbatim
(See also \ref PRINT \ref switchingfunction)
*/
//+ENDPLUMEDOC

//+PLUMEDOC MCOLVAR YDISTANCES
/*
Calculate the y components of the vectors connecting one or many pairs of atoms.  
You can then calculate functions of the distribution of
values such as the minimum, the number less than a certain quantity and so on.

\par Examples

See documentation for \ref XDISTANCES for examples of how to use this command.
You just need to substitute YDISTANCES for XDISTANCES to investigate the y component
rather than the x component.

*/
//+ENDPLUMEDOC

//+PLUMEDOC MCOLVAR ZDISTANCES
/*
Calculate the z components of the vectors connecting one or many pairs of atoms.  
You can then calculate functions of the distribution of
values such as the minimum, the number less than a certain quantity and so on.

\par Examples

See documentation for \ref XDISTANCES for examples of how to use this command.
You just need to substitute ZDISTANCES for XDISTANCES to investigate the z component
rather than the x component.

*/
//+ENDPLUMEDOC


class XDistances : public MultiColvar {
private:
  unsigned myc;
public:
  static void registerKeywords( Keywords& keys );
  explicit XDistances(const ActionOptions&);
// active methods:
  virtual double compute( const unsigned& tindex, AtomValuePack& myatoms ) const ;
/// Returns the number of coordinates of the field
  bool isPeriodic(){ return false; }
};

PLUMED_REGISTER_ACTION(XDistances,"XDISTANCES")
PLUMED_REGISTER_ACTION(XDistances,"YDISTANCES")
PLUMED_REGISTER_ACTION(XDistances,"ZDISTANCES")

void XDistances::registerKeywords( Keywords& keys ){
  MultiColvar::registerKeywords( keys );
  keys.use("ATOMS");  keys.use("MAX"); keys.use("ALT_MIN"); 
  keys.use("MEAN"); keys.use("MIN"); keys.use("LESS_THAN");
  keys.use("LOWEST"); keys.use("HIGHEST"); 
  keys.use("MORE_THAN"); keys.use("BETWEEN"); keys.use("HISTOGRAM"); keys.use("MOMENTS");
  keys.add("atoms-1","GROUP","Calculate the distance between each distinct pair of atoms in the group");
  keys.add("atoms-2","GROUPA","Calculate the distances between all the atoms in GROUPA and all "
                              "the atoms in GROUPB. This must be used in conjuction with GROUPB.");
  keys.add("atoms-2","GROUPB","Calculate the distances between all the atoms in GROUPA and all the atoms "
                              "in GROUPB. This must be used in conjuction with GROUPA.");
}

XDistances::XDistances(const ActionOptions&ao):
PLUMED_MULTICOLVAR_INIT(ao)
{
  if( getName().find("X")!=std::string::npos) myc=0;
  else if( getName().find("Y")!=std::string::npos) myc=1;
  else if( getName().find("Z")!=std::string::npos) myc=2;
  else plumed_error();

  // Read in the atoms
  std::vector<AtomNumber> all_atoms;
  readTwoGroups( "GROUP", "GROUPA", "GROUPB", all_atoms );
  int natoms=2; readAtoms( natoms, all_atoms );
  // And check everything has been read in correctly
  checkRead();
}

double XDistances::compute( const unsigned& tindex, AtomValuePack& myatoms ) const {
   Vector distance; 
   distance=getSeparation( myatoms.getPosition(0), myatoms.getPosition(1) );
   const double value=distance[myc];

   Vector myvec; myvec.zero(); 
   // And finish the calculation
   myvec[myc]=+1; addAtomDerivatives( 1, 1, myvec, myatoms );
   myvec[myc]=-1; addAtomDerivatives( 1, 0, myvec, myatoms );
   myatoms.addBoxDerivatives( 1, Tensor(distance,myvec) );
   return value;
}

}
}

