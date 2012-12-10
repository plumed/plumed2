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
#include "core/PlumedMain.h"
#include "ActionRegister.h"
#include "tools/PDB.h"
#include "tools/DRMSD.h"
#include "core/Atoms.h"

using namespace std;

namespace PLMD{

//+PLUMEDOC COLVAR DRMSD
/*
Calculate the distance RMSD with respect to a reference structure. 

To calculate the root-mean-square deviation between the atoms in two configurations
you must first superimpose the two structures in some ways.  Obviously, it is the internal vibrational 
motions of the structure - i.e. not the translations and rotations - that are interesting. However, 
aligning two structures by removing the translational and rotational motions is not easy.  Furthermore,
in some cases there can be alignment issues caused by so-called frame-fitting problems. It is thus
often cheaper and easier to calculate the distances between all the pairs of atoms.  The distance 
between the two structures, \f$\mathbf{X}^a\f$ and \f$\mathbf{X}^b\f$ can then be measured as:

\f[
d(\mathbf{X}^A, \mathbf{X}^B) = \frac{1}{N(N-1)} \sum_{i \ne j} [ d(\mathbf{x}_i^a,\mathbf{x}_j^a) - d(\mathbf{x}_i^b,\mathbf{x}_j^b) ]^2 
\f] 

where \f$N\f$ is the number of atoms and \f$d(\mathbf{x}_i,\mathbf{x}_j)\f$ represents the distance between
atoms \f$i\f$ and \f$j\f$.  Clearly, this representation of the configuration is invariant to translation and rotation.  
However, it can become expensive to calculate when the number of atoms is large.  This can be resolved
within the DRMSD colvar by setting LOWER_CUTOFF and UPPER_CUTOFF.  These keywords ensure that only
pairs of atoms that are within a certain range are incorporated into the above sum. 

In PDB files the atomic coordinates and box lengths should be in Angstroms unless 
you are working with natural units.  If you are working with natural units then the coordinates 
should be in your natural length unit.  For more details on the PDB file format visit http://www.wwpdb.org/docs.html

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
  keys.add("compulsory","REFERENCE","a file in pdb format containing the reference structure and the atoms involved in the CV.");
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
  if( !pdb.read(reference,plumed.getAtoms().usingNaturalUnits(),0.1/atoms.getUnits().getLength()) )
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
