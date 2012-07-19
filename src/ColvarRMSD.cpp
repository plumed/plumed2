/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012 The PLUMED team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

   This file is part of PLUMED, version 2.0.

   PLUMED is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   PLUMED is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with PLUMED.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "Colvar.h"
#include "PlumedMain.h"
#include "ActionRegister.h"
#include "PDB.h"
#include "RMSD.h"
#include "Atoms.h"

using namespace std;

namespace PLMD{

//+PLUMEDOC COLVAR RMSD
/*
Calculate the RMSD with respect to a reference structure.  To perform
an ??optimal?? (what does this mean algorithmical speed wise?) alignment
using the Kearsley algorithm then use TYPE=OPTIMAL.  Otherwise
use TYPE=SIMPLE, which will not perform optimal alignment and will only 
remove the translation of the center of mass.

\attention
The documentation here needs some work as it is not very clear to me 
sorry GAT. 

\par Examples

The following tells plumed to calculate the RMSD distance between
the positions of the atoms in the reference file and their instantaneous
position.  The Kearseley algorithm is used so this is done optimally.

\verbatim
RMSD REFERENCE=file.pdb TYPE=OPTIMAL
\endverbatim

...

*/
//+ENDPLUMEDOC

   
class ColvarRMSD : public Colvar {
	
  RMSD rmsd;
	
  vector<Vector> derivs;

public:
  ColvarRMSD(const ActionOptions&);
  virtual void calculate();
  static void registerKeywords(Keywords& keys);
};

PLUMED_REGISTER_ACTION(ColvarRMSD,"RMSD")

void ColvarRMSD::registerKeywords(Keywords& keys){
  Colvar::registerKeywords(keys);
  keys.add("compulsory","REFERENCE","a file in pdb format containing the reference structure and the atoms involved in the CV. " + PDB::documentation() );
  keys.add("compulsory","TYPE","SIMPLE","the manner in which RMSD alignment is performed.  Should be OPTIMAL or SIMPLE.");
}

ColvarRMSD::ColvarRMSD(const ActionOptions&ao):
PLUMED_COLVAR_INIT(ao),rmsd(log)
{
  string reference;
  parse("REFERENCE",reference);
  string type;	
  type.assign("SIMPLE");
  parse("TYPE",type);

  checkRead();


  addValueWithDerivatives(); setNotPeriodic();
  PDB pdb;

  // read everything in ang and transform to nm if we are not in natural units
  pdb.read(reference,plumed.getAtoms().usingNaturalUnits(),0.1/atoms.getUnits().length);

  rmsd.set(pdb,type);

  requestAtoms(pdb.getAtomNumbers());

  derivs.resize(getNumberOfAtoms());
  log.printf("  reference from file %s\n",reference.c_str());
  log.printf("  which contains %d atoms\n",getNumberOfAtoms());
  log.printf("  method for alignment : %s \n",rmsd.getMethod().c_str() );

}


// calculator
void ColvarRMSD::calculate(){
  double r=rmsd.calculate(getPositions(),derivs);
  setValue(r);
  for(unsigned i=0;i<derivs.size();i++) setAtomsDerivatives(i,derivs[i]);
  Tensor virial;
  for(unsigned i=0;i<derivs.size();i++) virial=virial+(-1.0*Tensor(getPosition(i),derivs[i]));
  setBoxDerivatives(virial);
}

}



