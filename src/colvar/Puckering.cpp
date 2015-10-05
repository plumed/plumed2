/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

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
#include "Colvar.h"
#include "core/PlumedMain.h"
#include "ActionRegister.h"
#include "tools/Torsion.h"

#include <string>
#include <cmath>
#include <iostream>

using namespace std;

namespace PLMD{
namespace colvar{

//+PLUMEDOC COLVAR PUCKERING
/*
 Calculate the Nucleic Acid sugar pseudorotation coordinates.
 
 This command can be used to represent the ring's pseudorotation with polar coordinates (Phase, Amplitud) or cartesian coordinates (Zx, Zy) using the position of the five ring atoms ( see \cite huang2014improvement ).
 
 Components of this action are:
 
 \par Examples
 
 This input tells plumed to print the puckering phase angle of the 3rd nucleotide of a RNA molecule on file COLVAR.
 \verbatim
 MOLINFO STRUCTURE=rna.pdb MOLTYPE=rna
 PUCKERING ATOMS=@sugar-3 LABEL=puck
 PRINT ARG=puck.phs FILE=COLVAR
 \endverbatim

*/
//+ENDPLUMEDOC
   
class Puckering : public Colvar {

public:
  explicit Puckering(const ActionOptions&);
  virtual void calculate();
  static void registerKeywords(Keywords& keys);
};

PLUMED_REGISTER_ACTION(Puckering,"PUCKERING")

void Puckering::registerKeywords(Keywords& keys){
   Colvar::registerKeywords( keys );
   keys.remove("NOPBC");
   keys.add("atoms","ATOMS","the five atoms of the sugar ring in the order C4',O4',C1',C2',C3'");
   componentsAreNotOptional(keys);
   keys.addOutputComponent("phs","default","Pseudorotaion phase");
   keys.addOutputComponent("amp","default","Pseudorotation amplitude");
   keys.addOutputComponent("Zx","default","Pseudorotation x cartesian component");
   keys.addOutputComponent("Zy","default","Pseudorotation y cartesian component");
}

Puckering::Puckering(const ActionOptions&ao):
PLUMED_COLVAR_INIT(ao)
{
  vector<AtomNumber> atoms;
  parseAtomList("ATOMS",atoms);
  if(atoms.size()!=5) error("only for 5-membered rings");
  checkRead();
  
  plumed.cite("Huang, Giese, Lee, York, J. Chem. Theory Comput. 10, 1538 (2014)");
  
  if(atoms.size()==5){
    log.printf("  between atoms %d %d %d %d %d\n",atoms[0].serial(),atoms[1].serial(),atoms[2].serial(),atoms[3].serial(),atoms[4].serial());
  }else error("ATOMS should specify 5 atoms");
    
  addComponentWithDerivatives("phs"); componentIsPeriodic("phs","-pi","pi");
  addComponentWithDerivatives("amp"); componentIsNotPeriodic("amp");
  addComponentWithDerivatives("Zx"); componentIsNotPeriodic("Zx");
  addComponentWithDerivatives("Zy"); componentIsNotPeriodic("Zy");

  requestAtoms(atoms);
}

// calculator
void Puckering::calculate(){

  Vector d0,d1,d2,d3,d4,d5;
  makeWhole();
  
  d0=delta(getPosition(2),getPosition(1));
  d1=delta(getPosition(3),getPosition(2));
  d2=delta(getPosition(4),getPosition(3));
  d3=delta(getPosition(4),getPosition(3));
  d4=delta(getPosition(0),getPosition(4));
  d5=delta(getPosition(1),getPosition(0));
    
  Vector r[5];
  r[0]=getPosition(0);
  for(unsigned i=1;i<5;i++){
    r[i]=r[i-1]+pbcDistance(getPosition(i-1),getPosition(i));
  }
    
  Vector dd0,dd1,dd2,dd3,dd4,dd5;
    
  PLMD::Torsion t;
    
  double v1=t.compute(d0,d1,d2,dd0,dd1,dd2);
  double v3=t.compute(d3,d4,d5,dd3,dd4,dd5);

  double Zx=(v1+v3)/(2.0*cos(4.0*pi/5.0));
  double Zy=(v1-v3)/(2.0*sin(4.0*pi/5.0));
  double phase=atan2(Zy,Zx);
  double amplitude=sqrt(Zx*Zx+Zy*Zy);
    
  Vector dZx_dR[5];
  Vector dZy_dR[5];
    
  dZx_dR[0]=(dd5-dd4);
  dZx_dR[1]=(dd0-dd5);
  dZx_dR[2]=(dd1-dd0);
  dZx_dR[3]=(dd2+dd3-dd1);
  dZx_dR[4]=(dd4-dd3-dd2);
 
  dZy_dR[0]=(dd4-dd5);
  dZy_dR[1]=(dd0+dd5);
  dZy_dR[2]=(dd1-dd0);
  dZy_dR[3]=(dd2-dd3-dd1);
  dZy_dR[4]=(dd3-dd4-dd2);
        
  for(unsigned j=0;j<5;j++) dZx_dR[j]*=(1.0/(2.0*cos(4.0*pi/5.0)));
  for(unsigned j=0;j<5;j++) dZy_dR[j]*=(1.0/(2.0*sin(4.0*pi/5.0)));
    
  Vector dphase_dR[5];
    for(unsigned j=0;j<5;j++) dphase_dR[j]=(1.0/(Zx*Zx+Zy*Zy))*(-Zy*dZx_dR[j] + Zx*dZy_dR[j]);  
  
  Vector damplitude_dR[5];
    for(unsigned j=0;j<5;j++) damplitude_dR[j]=(1.0/amplitude)*(Zx*dZx_dR[j] + Zy*dZy_dR[j]);   
    
  Value* vzx=getPntrToComponent("Zx");
  vzx->set(Zx);
  setAtomsDerivatives (vzx,0, dZx_dR[0]);
  setAtomsDerivatives (vzx,1, dZx_dR[1]);
  setAtomsDerivatives (vzx,2, dZx_dR[2]);
  setAtomsDerivatives (vzx,3, dZx_dR[3]);
  setAtomsDerivatives (vzx,4, dZx_dR[4]);
  setBoxDerivativesNoPbc(vzx);
    
  Value* vzy=getPntrToComponent("Zy");
  vzy->set(Zy);
  setAtomsDerivatives (vzy,0, dZy_dR[0]);
  setAtomsDerivatives (vzy,1, dZy_dR[1]);
  setAtomsDerivatives (vzy,2, dZy_dR[2]);
  setAtomsDerivatives (vzy,3, dZy_dR[3]);
  setAtomsDerivatives (vzy,4, dZy_dR[4]);
  setBoxDerivativesNoPbc(vzy);


  Value* vph=getPntrToComponent("phs");
  vph->set(phase);
  setAtomsDerivatives (vph,0, dphase_dR[0]);
  setAtomsDerivatives (vph,1, dphase_dR[1]);
  setAtomsDerivatives (vph,2, dphase_dR[2]);
  setAtomsDerivatives (vph,3, dphase_dR[3]);
  setAtomsDerivatives (vph,4, dphase_dR[4]);
  setBoxDerivativesNoPbc(vph);
    
  Value* vam=getPntrToComponent("amp");
  vam->set(amplitude);
  setAtomsDerivatives (vam,0, damplitude_dR[0]);
  setAtomsDerivatives (vam,1, damplitude_dR[1]);
  setAtomsDerivatives (vam,2, damplitude_dR[2]);
  setAtomsDerivatives (vam,3, damplitude_dR[3]);
  setAtomsDerivatives (vam,4, damplitude_dR[4]);
  setBoxDerivativesNoPbc(vam);

    
}

}
}



