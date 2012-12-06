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
#include "tools/PDB.h"
#include "tools/RMSD.h"
#include "Atoms.h"


using namespace std;

namespace PLMD{
   
class ColvarRMSD : public Colvar {
	
  RMSD rmsd;
	
  bool squared; 

  vector<Vector> derivs;

public:
  ColvarRMSD(const ActionOptions&);
  virtual void calculate();
  static void registerKeywords(Keywords& keys);
};

}


using namespace std;

namespace PLMD{

//+PLUMEDOC COLVAR RMSD
/*
Calculate the RMSD with respect to a reference structure.  

To calculate the root-mean-square deviation between the atoms in two configurations
you must first superimpose the two structures in some ways.  Obviously, it is the internal vibrational 
motions of the structure - i.e. not the translations and rotations - that are interesting.  It is 
possible to align two structures (i.e. remove the translational and rotational motions of the 
system) using the Kearsley \cite kearsley algorithm.  This algorithm calculates a 
roto-translation matrix from a set of \e alignment atoms \f${\cal A}\f$.  The amount by which 
a \e displacement set of atoms, \f${\cal B}\f$, has been moved can then be calculated using:   

\f[
d(X_j,X_i) = \sum_{a=1}^{N_{\cal B}} (X^{(j)}_a - M_{ij} X^{(i)}_a)^2
\f]

where \f$M_{ij}\f$ is the roto-translation matrix translated by the Kearsley algorithm.

When a pdb input file is required to specify the structure to align to in an RMSD calculation
the occupancy and beta columns are used to specify which atoms form part of the alignment set
and which atoms form part of the displacement set.  Values of 1.0 and 0.0 indicate that the 
atom is to be used in the \e alignment set only, while values of 0.0 and 1.0 indicates that 
the atom is to be used in the \e displacement set only.  Values of 1.0 and 1.0 indicate that 
the atom is to be used in both the \e alignment and \e displacement sets. Users can also 
use fractional values for beta and the occupancy values.  We recommend you only do this when
you really know what you are doing however as the results can be rather strange.  
When this form of RMSD is used to calculate the secondary structure variables (\ref ALPHARMSD,
\ref ANTIBETARMSD and \ref PARABETARMSD) all the atoms in the segment are assumed to be part of
both the \e alignment and \e displacement sets.

The \e OPTIMAL type of rmsd calculation removes both translational and rotational motions.  
There may, however, be ocasions where you do not need to remove the rotational motions.  That is
to say there may be times when you only need to align the centers of mass of the two structures.
This can be done by performing the \e SIMPLE form of RMSD calculation.

In PDB files the atomic coordinates and box lengths should be in Angstroms unless 
you are working with natural units.  If you are working with natural units then the coordinates 
should be in your natural length unit.  For more details on the PDB file format visit http://www.wwpdb.org/docs.html

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

PLUMED_REGISTER_ACTION(ColvarRMSD,"RMSD")

void ColvarRMSD::registerKeywords(Keywords& keys){
  Colvar::registerKeywords(keys);
  keys.add("compulsory","REFERENCE","a file in pdb format containing the reference structure and the atoms involved in the CV.");
  keys.add("compulsory","TYPE","SIMPLE","the manner in which RMSD alignment is performed.  Should be OPTIMAL or SIMPLE.");
  keys.addFlag("SQUARED",false," This should be setted if you want MSD instead of RMSD ");
}

ColvarRMSD::ColvarRMSD(const ActionOptions&ao):
PLUMED_COLVAR_INIT(ao),rmsd(log),squared(false)
{
  string reference;
  parse("REFERENCE",reference);
  string type;	
  type.assign("SIMPLE");
  parse("TYPE",type);
  parseFlag("SQUARED",squared);

  checkRead();


  addValueWithDerivatives(); setNotPeriodic();
  PDB pdb;

  // read everything in ang and transform to nm if we are not in natural units
  if( !pdb.read(reference,plumed.getAtoms().usingNaturalUnits(),0.1/atoms.getUnits().getLength()) )
      error("missing input file " + reference );

  rmsd.set(pdb,type);

  requestAtoms(pdb.getAtomNumbers());

  derivs.resize(getNumberOfAtoms());
  log.printf("  reference from file %s\n",reference.c_str());
  log.printf("  which contains %d atoms\n",getNumberOfAtoms());
  log.printf("  method for alignment : %s \n",rmsd.getMethod().c_str() );
  if(squared)log.printf("  chosen to use SQUARED option for MSD instead of RMSD\n");

}


// calculator
void ColvarRMSD::calculate(){
  double r=rmsd.calculate(getPositions(),derivs,squared);
  setValue(r);
  for(unsigned i=0;i<derivs.size();i++) setAtomsDerivatives(i,derivs[i]);
  Tensor virial;
  for(unsigned i=0;i<derivs.size();i++) virial=virial+(-1.0*Tensor(getPosition(i),derivs[i]));
  setBoxDerivatives(virial);
}

}



