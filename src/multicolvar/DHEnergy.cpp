/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013,2014 The plumed team
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
#include "vesselbase/VesselRegister.h"
#include "vesselbase/FunctionVessel.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"
#include "MultiColvar.h"

namespace PLMD {
namespace multicolvar {

class DHEnergy : public vesselbase::FunctionVessel {
private:
  MultiColvar* mycolv;
  double I, T;
  double k; // Inverse Debye screening length
  double constant;
  double epsilon;
  std::vector<double> df;
public:
  static void registerKeywords( Keywords& keys );
  static void reserveKeyword( Keywords& keys );
  DHEnergy( const vesselbase::VesselOptions& da );
  std::string function_description();
  bool calculate();
  void finish();
};

PLUMED_REGISTER_VESSEL(DHEnergy,"DHENERGY")

void DHEnergy::registerKeywords( Keywords& keys ){
  FunctionVessel::registerKeywords( keys );
  keys.add("compulsory","I","1.0","Ionic strength (M)");
  keys.add("compulsory","TEMP","300.0","Simulation temperature (K)");
  keys.add("compulsory","EPSILON","80.0","Dielectric constant of solvent");
}

void DHEnergy::reserveKeyword( Keywords& keys ){
  keys.reserve("numbered","DHENERGY","calculate the Debye-Huckel interaction energy. This is a alternative "
                                     "implementation of \\ref DHENERGY that is particularly useful if you "
                                     "want to calculate the Debye-Huckel interaction energy and some other "
                                     "function of set of distances between the atoms in the two groups. "
                                     "The input for this keyword should read "
                                     "DHENERGY={I=\\f$I\\f$ TEMP=\\f$T\\f$ EPSILON=\\f$\\epsilon\\f$}.");
  keys.addOutputComponent("dhenergy","DHENERGY","the Debye-Huckel interaction energy. You can calculate "
                                                "this quantity multiple times using different parameters");
}

DHEnergy::DHEnergy( const vesselbase::VesselOptions& da ) :
FunctionVessel(da),
df(2)
{
  mycolv=dynamic_cast<MultiColvar*>( getAction() );
  plumed_massert( mycolv, "DHENERGY can only be used with MultiColvars and should only be used with DISTANCES");

  parse("I",I); parse("TEMP",T); parse("EPSILON",epsilon);
  
  Atoms& catoms( getAction()->plumed.getAtoms() ); 
  if( catoms.usingNaturalUnits() ) error("DHENERGY cannot be used for calculations performed with natural units");
  constant=138.935458111/catoms.getUnits().getEnergy()/catoms.getUnits().getLength();
  k=sqrt(I/(epsilon*T))*502.903741125*catoms.getUnits().getLength();
}

std::string DHEnergy::function_description(){
  std::ostringstream ostr;
  ostr<<"the Debye-Huckel interaction energy "<<getAction()->plumed.cite("Do, Carloni, Varani and Bussi, J. Chem. Theory Comput. 9, 1720 (2013)")<<"."; 
  ostr<<" Parameters : temperature "<<T<<" K, ionic strength "<<I<<" M, ";
  ostr<<"solvent dielectric constant "<<epsilon;
  return ostr.str();
}

bool DHEnergy::calculate(){
  if( mycolv->getAbsoluteIndex(0)==mycolv->getAbsoluteIndex(1) ) return false;

  double val=getAction()->getElementValue(0);
  double invdistance = 1.0 / val;
  double f=exp(-k*val)*invdistance*constant*mycolv->getCharge(0)*mycolv->getCharge(1)/epsilon;
  double dval=-(k+invdistance)*f; 
  addValueIgnoringTolerance(0,f);
  getAction()->chainRuleForElementDerivatives( 0, 0, dval, this ); 
  return true;
}

void DHEnergy::finish(){
  setOutputValue( getFinalValue(0) ); 
  df[0]=1.0; df[1]=0.0;
  mergeFinalDerivatives( df );
}

} 
}
