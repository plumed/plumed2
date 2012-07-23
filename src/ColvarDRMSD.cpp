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
#include "PlumedMain.h"
#include "ActionRegister.h"
#include "PDB.h"
#include "DRMSD.h"
#include "Atoms.h"

using namespace std;

namespace PLMD{

//+PLUMEDOC COLVAR DRMSD
/*
Calculate the distance RMSD with respect to a reference structure. 
Only pairs of atoms whose distance in the reference structure is within 
LOWER_CUTOFF and UPPER_CUTOFF are considered.

\par Examples

The following tells plumed to calculate the distance RMSD between
the positions of the atoms in the reference file and their instantaneous
position. Only pairs of atoms whose distance in the reference structure is within 
0.1 and 0.8 nm are considered. 

\verbatim
DRMSD REFERENCE=file.pdb LOWER_CUTOFF=0.1 UPPER_CUTOFF=0.8
\endverbatim

...

*/
//+ENDPLUMEDOC

   
class ColvarDRMSD : public Colvar {
	
  vector<Vector> derivs_;
  DRMSD drmsd_;
  bool pbc_;

public:
  ColvarDRMSD(const ActionOptions&);
  virtual void calculate();
  static void registerKeywords(Keywords& keys);
};

PLUMED_REGISTER_ACTION(ColvarDRMSD,"DRMSD")

void ColvarDRMSD::registerKeywords(Keywords& keys){
  Colvar::registerKeywords(keys);
  keys.add("compulsory","REFERENCE","a file in pdb format containing the reference structure and the atoms involved in the CV. " + PDB::documentation() );
  keys.add("compulsory","LOWER_CUTOFF","only pairs of atoms further than LOWER_CUTOFF are considered in the calculation.");
  keys.add("compulsory","UPPER_CUTOFF","only pairs of atoms closer than UPPER_CUTOFF are considered in the calculation.");
}

ColvarDRMSD::ColvarDRMSD(const ActionOptions&ao):
PLUMED_COLVAR_INIT(ao), pbc_(true)
{
  string reference;
  parse("REFERENCE",reference);
  double lcutoff;	
  parse("LOWER_CUTOFF",lcutoff);
  double ucutoff;       
  parse("UPPER_CUTOFF",ucutoff);
  bool nopbc(false);
  parseFlag("NOPBC",nopbc);
  pbc_=!nopbc;

  checkRead();

  addValueWithDerivatives(); setNotPeriodic(); 

  // read everything in ang and transform to nm if we are not in natural units
  PDB pdb;
  if( !pdb.read(reference,plumed.getAtoms().usingNaturalUnits(),0.1/atoms.getUnits().length) )
      error("missing input file " + reference );

  // store target_ distance
  drmsd_.setFromPDB(pdb, lcutoff, ucutoff);

  requestAtoms(pdb.getAtomNumbers());
  derivs_.resize(getNumberOfAtoms());

  log.printf("  reference from file %s\n",reference.c_str());
  log.printf("  which contains %d atoms\n",getNumberOfAtoms());

}

// calculator
void ColvarDRMSD::calculate(){

// set derivatives to zero
 for(unsigned i=0;i<derivs_.size();++i) {derivs_[i].zero();}

 double drmsd;
 Tensor virial;

 if(pbc_){drmsd=drmsd_.calculate(getPositions(),getPbc(),derivs_,virial);}
 else{    drmsd=drmsd_.calculate(getPositions(),         derivs_,virial);}

 setValue(drmsd);

 for(unsigned i=0;i<derivs_.size();++i) {setAtomsDerivatives(i,derivs_[i]);}
 
 setBoxDerivatives(virial);

 }

}
