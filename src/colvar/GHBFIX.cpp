/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2021-2023 The plumed team
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
#include "CoordinationBase.h"
#include "tools/SwitchingFunction.h"
#include "ActionRegister.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"
#include "tools/IFile.h"

#include <iostream>

#include <string>

namespace PLMD {
namespace colvar {

//+PLUMEDOC COLVAR GHBFIX
/*
Calculate the GHBFIX interaction energy among GROUPA and GROUPB
using a potential defined in Kührová et al., Improving the performance of the AMBER RNA force field by
tuning the hydrogen-bonding interactions, JCTC, 2019. Essentially it is a switching function being -1 for small distances and 0 for large distances with a smooth interpolation in the middle. This can be scaled as desired by specifying interaction scaling parameters and energy units.

This collective variable can be used to analyze hydrogen bond interactions, or to generate bias potentials.
Notice that the value of the GHBFIX is returned in plumed units (see \ref UNITS), if not specified differently via ENERGY_UNITS.

\par Examples
This example prints the GHBFIX interaction in kcal/mol between two groups of atoms using D_0, D_MAX and C
It is applied in the functional form introduced in the pioneering paper.
The types of atoms 1-6 should be defined in typesTable_examples.dat while their interaction parameters should be defined in scalingParameters_examples.dat in kBT units.

\plumedfile
#SETTINGS AUXFOLDER=regtest/basic/rt-ghbfix
gh: GHBFIX PAIR GROUPA=1,2,3 GROUP=4,5,6 D_0=0.2 D_MAX=0.3 C=0.8 TYPES=typesTable_examples.dat PARAMS=scalingParameters_examples.dat ENERGY_UNITS=kcal/mol
PRINT FILE=output ARG=gh
\endplumedfile

*/
//+ENDPLUMEDOC

class GHBFIX : public CoordinationBase {

  double dmax;
  double dmax_squared;
  double d0;
  double c;

  std::vector<unsigned> typesTable;

  std::vector<double> etas;

  unsigned n;

  double dmax2;
  double A;
  double B;
  double C;
  double D;

public:
  explicit GHBFIX(const ActionOptions&);
// active methods:
  static void registerKeywords( Keywords& keys );
  double pairing(double distance,double&dfunc,unsigned i,unsigned j)const override;
};

PLUMED_REGISTER_ACTION(GHBFIX,"GHBFIX")

void GHBFIX::registerKeywords( Keywords& keys ) {
  CoordinationBase::registerKeywords(keys);

  keys.add("optional","ENERGY_UNITS","the value of ENERGY_UNITS in the switching function");
  keys.add("compulsory","TYPES","the value of TYPES in the switching function");
  keys.add("compulsory","PARAMS","the value of PARAMS in the switching function");
  keys.add("compulsory","D_MAX","the value of D_MAX in the switching function");
  keys.add("compulsory","D_0","the value of D_0 in the switching function");
  keys.add("compulsory","C","the value of C in the switching function");
}

GHBFIX::GHBFIX(const ActionOptions&ao):
  Action(ao),
  CoordinationBase(ao)
{
  std::string types;
  std::string params;
  std::string energy_units ="plumed" ;

  parse("D_MAX",dmax);
  dmax_squared = dmax*dmax;
  parse("D_0",d0);
  parse("C",c);
  parse("TYPES",types);
  parse("PARAMS",params);
  parse("ENERGY_UNITS",energy_units);

  dmax2 = dmax-d0;

  A = (-c*dmax2*dmax2)/((1-c)*dmax2*dmax2);
  B = (2*dmax2)/((1-c)*dmax2*dmax2);
  C = -1/((1-c)*dmax2*dmax2);
  D = 1/(c*dmax2*dmax2);

  std::map<std::string,unsigned> MapTypesTable;

  //typesTable
  IFile typesfile;
  typesfile.link(*this);
  typesfile.open(types);
  std::string itype;
  while(typesfile.scanField("itype",itype).scanField()) {
    plumed_assert(itype.empty()==false)<<"itype is empty";

    if (MapTypesTable.empty()) {
      MapTypesTable.insert({itype, 0});
    }
    else if (MapTypesTable.count(itype) == 0) {
      unsigned currentMax = 0;
      for(auto it = MapTypesTable.cbegin(); it != MapTypesTable.cend(); ++it ) {
        if (it ->second > currentMax) {
          currentMax = it->second;
        }
      }
      MapTypesTable.insert({itype, currentMax+1});
    }

    typesTable.push_back(MapTypesTable[itype]);
  }

  n = (int)*std::max_element(std::begin(typesTable), std::end(typesTable));
  n+=1;

  //scalingParameters
  etas.resize(n*n,0.0);
  IFile etafile;
  etafile.open(params);
  std::string it,jt;
  double eta;
  while(etafile.scanField("itype",it).scanField("jtype",jt).scanField("eta",eta).scanField()) {
    plumed_assert(it.empty()==false)<<"itype is empty";
    plumed_assert(jt.empty()==false)<<"jtype is empty";
    etas[n*MapTypesTable[it]+MapTypesTable[jt]]=eta;
  }

  if(energy_units!="plumed") {
    Units units;
    units.setEnergy(energy_units);
    for(auto i=0; i<etas.size(); i++) etas[i]*=units.getEnergy()/atoms.getUnits().getEnergy();
  }

}


double GHBFIX::pairing(double distance2,double&dfunc,unsigned i,unsigned j)const {

  const auto i1=getAbsoluteIndex(i).index();
  plumed_assert(i1<typesTable.size())<<"your types table only covers "<<typesTable.size()<<" atoms, but you are trying to access atom number "<<(i1+1);
  const auto t1=typesTable[i1];

  const auto i2=getAbsoluteIndex(j).index();
  plumed_assert(i2<typesTable.size())<<"your types table only covers "<<typesTable.size()<<" atoms, but you are trying to access atom number "<<(i2+1);
  const auto t2=typesTable[i2];

  const double scale=etas[n*t1+t2];

  double result;
  if(distance2>dmax_squared) {
    result=0.;
    dfunc=0.0;
    return result;
  }
  double distance=std::sqrt(distance2);
  const double rdist = (distance-d0);

  if(rdist<=0.) {
    result=-1.;
    dfunc=0.0;
  } else {
    result=-1.;
    dfunc=0.0;

    if (rdist > c*dmax2) {
      result+=(A + B*rdist + C*rdist*rdist);
      dfunc+=B+2*C*rdist;
    } else if (rdist > 0.0) {
      result+=D*(rdist*rdist);
      dfunc+=2*D*rdist;
    }

    dfunc/=distance;
  }

  result*=scale;
  dfunc*=scale;

  return result;
}

}

}
