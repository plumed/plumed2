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
#include "vesselbase/VesselRegister.h"
#include "vesselbase/SumVessel.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"
#include "MultiColvar.h"

namespace PLMD {
namespace multicolvar {

class VesselDHEnergy : public vesselbase::SumVessel {
private:
  MultiColvar* mycolv;
  double k; // Inverse Debye screening length
  double constant;
  double epsilon;
public:
  static void reserveKeyword( Keywords& keys );
  VesselDHEnergy( const vesselbase::VesselOptions& da );
  void printKeywords();
  double compute( const unsigned& i, const double& val, double& df );
};

PLUMED_REGISTER_VESSEL(VesselDHEnergy,"DHENERGY")

void VesselDHEnergy::reserveKeyword( Keywords& keys ){
  keys.reserve("numbered","DHENERGY","calculate the Debye-Huckel interaction energy. This is a alternative "
                                     "implementation of \\ref DHENERGY that is particularly useful if you "
                                     "want to calculate the Debye-Huckel interaction energy and some other "
                                     "function of set of distances between the atoms in the two groups. "
                                     "The input for this keyword should read "
                                     "DHENERGY={I=\\f$I\\f$ TEMP=\\f$T\\f$ EPSILON=\\f$\\epsilon\\f$}. "
                                     "The final value can be referenced using \\e label.dhenergy.");
}

VesselDHEnergy::VesselDHEnergy( const vesselbase::VesselOptions& da ) :
SumVessel(da)
{
  mycolv=dynamic_cast<MultiColvar*>( getAction() );
  plumed_massert( mycolv, "DHENERGY can only be used with MultiColvars and should only be used with DISTANCES");
  std::vector<std::string> data=Tools::getWords(da.parameters);

  double I,T;
  if( !Tools::parse(data,"I",I) ) I=1.0; 
  if( !Tools::parse(data,"TEMP",T) ) T=300.0; 
  if( !Tools::parse(data,"EPSILON",epsilon) ) epsilon=80.0;
  if( !data.empty() ){
      std::string errormsg="found the following rogue keywords in input to DHEnergy function input : ";
      for(unsigned i=0;i<data.size();++i) errormsg = errormsg + data[i] + " ";
      error( errormsg );
  }
  
  Atoms& catoms( getAction()->plumed.getAtoms() ); 
  if( catoms.usingNaturalUnits() ) error("DHENERGY cannot be used for calculations performed with natural units");
  constant=138.935458111/catoms.getUnits().getEnergy()/catoms.getUnits().getLength();
  k=sqrt(I/(epsilon*T))*502.903741125*catoms.getUnits().getLength();

  std::string label,name; unsigned numlab=getNumericalLabel();
  if( numlab==0 ){
     name="dhenergy"; 
  } else {
     Tools::convert(numlab,label);
     name="dhenergy-" + label;
  }
  addOutput(name);
  log.printf("  value %s.%s contains the Debye-Huckel interaction energy\n",(getAction()->getLabel()).c_str(), name.c_str());
  log.printf("  parameters temperature %f K ionic strength %f M solvent dielectric constant %f \n",T, I, epsilon); 
  log<<"  Bibliography "<<getAction()->plumed.cite("Do, Carloni, Varani and Bussi, submitted (2013)")<<" \n";
}

void VesselDHEnergy::printKeywords(){
  Keywords dhkeys;
  dhkeys.add("compulsory","I","1.0","Ionic strength (M)");
  dhkeys.add("compulsory","TEMP","300.0","Simulation temperature (K)");
  dhkeys.add("compulsory","EPSILON","80.0","Dielectric constant of solvent");
  dhkeys.print(log); 
}

double VesselDHEnergy::compute( const unsigned& i, const double& val, double& df ){
  plumed_dbg_assert( i==0 );
  if( mycolv->getAbsoluteIndex(0)==mycolv->getAbsoluteIndex(1) ){
      df=0.0; return 0.0;
  }

  double invdistance = 1.0 / val;
  double f=exp(-k*val)*invdistance*constant*mycolv->getCharge(0)*mycolv->getCharge(1)/epsilon;
  df=-(k+invdistance)*f;
  return f;
}

} 
}
