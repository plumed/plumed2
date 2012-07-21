/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

   This file is part of plumed, version 2.0.

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
#include "ActionRegister.h"

#include <string>
#include <cmath>
#include <cassert>

using namespace std;

namespace PLMD{

//+PLUMEDOC COLVAR RGYR
/*
Calculate the radius of gyration for a chain of atoms.

The radius of gyration is calculated using:

\f[
s_{\rm Gyr}=\Big ( \frac{\sum_i^{n}  
 \vert {r}_i -{r}_{\rm COM} \vert ^2 }{\sum_i^{n} m_i} \Big)^{1/2} 
\f]

with the position of the center of mass \f${r}_{\rm COM}\f$ given by:

\f[
{r}_{\rm COM}=\frac{\sum_i^{n} {r}_i\ m_i }{\sum_i^{n} m_i}
\f]

\bug This was a very quick implementation of RGYR for a project that I am working on. It has very little of the functionality that is available in plumed 1.0. 

\par Examples

The following input tells plumed to print the radius of gyration of the 
chain containing atoms 10 to 20.
\verbatim
RGYR ATOMS=10-20 LABEL=rg
PRINT ARG=rg STRIDE=1 FILE=colvar 
\endverbatim
(See also \ref PRINT)

*/
//+ENDPLUMEDOC

class ColvarRGYR : public Colvar {
private:
  bool use_masses;
public:
  static void registerKeywords( Keywords& keys );
  ColvarRGYR(const ActionOptions&);
// active methods:
  virtual void calculate();
};

PLUMED_REGISTER_ACTION(ColvarRGYR,"RGYR")

void ColvarRGYR::registerKeywords( Keywords& keys ){
  Colvar::registerKeywords( keys );
  keys.add("atoms","ATOMS","the pair of atom that we are calculating the distance between");
  keys.addFlag("NOT_MASS_WEIGHTED",false,"set the masses of all the atoms equal to one");
}

ColvarRGYR::ColvarRGYR(const ActionOptions&ao):
PLUMED_COLVAR_INIT(ao)
{
  std::vector<AtomNumber> atoms;
  parseAtomList("ATOMS",atoms);
  if(atoms.size()==0) error("no atoms specified");
  parseFlag("NOT_MASS_WEIGHTED",use_masses);
  checkRead();

  log.printf("  atoms involved : ");
  for(unsigned i=0;i<atoms.size();++i) log.printf("%d ",atoms[i].serial());
  log.printf("\n");

  addValueWithDerivatives(); setNotPeriodic();
  requestAtoms(atoms);
}

void ColvarRGYR::calculate(){
  std::vector<Vector> derivatives( getNumberOfAtoms() );

  // Find the center of mass
  double totmass = 0; Vector pos0, com, diff;
  if( use_masses ) totmass += getMass(0); 
  else totmass += 1.0;
  pos0=getPosition(0); com.zero();
  for(unsigned i=1;i<getNumberOfAtoms();++i){
     diff=delta( pos0, getPosition(i) );
     if( use_masses ){
         totmass += getMass(i);
         com += getMass(i)*diff;
     } else {
         totmass += 1.0;
         com += diff;
     }
  }
  com = com / totmass + pos0;

  // Now compute radius of gyration
  double d, rgyr=0;
  for(unsigned i=0;i<getNumberOfAtoms();++i){
     diff=delta( com, getPosition(i) );
     d=diff.modulo();
     if( use_masses ){
        rgyr += getMass(i)*d*d;
        derivatives[i]=diff*getMass(i);
     } else {
        rgyr += d*d;
        derivatives[i]=diff;
     } 
  } 
  
  rgyr=sqrt(rgyr/totmass); setValue(rgyr);
  Tensor virial; virial.zero();
  for(unsigned i=0;i<getNumberOfAtoms();++i){
     derivatives[i] /= rgyr*totmass;
     setAtomsDerivatives(i,derivatives[i]);
     virial=virial+(-1.0*Tensor(getPosition(i),derivatives[i]));
  }
  setBoxDerivatives(virial);
}

}
