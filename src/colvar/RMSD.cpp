/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013 The plumed team
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
#include "tools/PDB.h"
#include "reference/RMSDBase.h"
#include "reference/MetricRegister.h"
#include "core/Atoms.h"


using namespace std;

namespace PLMD{
namespace colvar{
   
class RMSD : public Colvar {
	
  PLMD::RMSDBase* rmsd;
  bool squared; 

public:
  RMSD(const ActionOptions&);
  ~RMSD();
  virtual void calculate();
  static void registerKeywords(Keywords& keys);
};


using namespace std;

//+PLUMEDOC DCOLVAR RMSD
/*
Calculate the RMSD with respect to a reference structure.  

To calculate the root-mean-square deviation between the atoms in two configurations
you must first superimpose the two structures in some ways.  Obviously, it is the internal vibrational 
motions of the structure - i.e. not the translations and rotations - that are interesting.  It is 
possible to align two structures (i.e. remove the translational and rotational motions of the 
system).
In general the aim of this colvar is to calculate something like:

\f[
d(X,X_r) = \vert X-X' \vert 
\f]

where \f$ X \f$ is the molecular dynamics snapshot and   
\f$ X' \f$ is a reference structure you provide as input.
Here the symbols \f$ \vert \dots \vert \f$ represent some sort of norm of choice that is selected through the OPTIMAL switch (see below for examples).

Here we discuss the switch we implemented so far.

       - If you use TYPE=SIMPLE you just calculate the root mean square deviation after reset of the geometric center  
         that reads:
         \f[
         d(X,X_r) = \sqrt{ \sum_i \sum_\alpha^{x,y,z}  \frac{w_i}{\sum_j w_j}( X_{i,\alpha}-com_\alpha(X)-{X'}_{i,\alpha}+com_\alpha(X') )^2 } 
         \f]
         with 
         \f[
         com_\alpha(X)= \sum_i  \frac{w'_{i}}{\sum_j w'_j}X_{i,\alpha}
         \f]
         and
         \f[
         com_\alpha(X')= \sum_i  \frac{w'_{i}}{\sum_j w'_j}X'_{i,\alpha}
         \f]
         where  \f$ com_\alpha(X) \f$ and  \f$ com_\alpha(X') \f$  are the center of mass
         whenever the weights $w'$ are assigned proportional to the atomic masses, otherwise are the geometric centers of the 
         running MD snapshot and the reference frame respectively whenever they are all set to one.
         Note that two set of weights exist:  \f$ w' \f$ and  \f$ w \f$. The first is used to calculate the center of mass (so it determines the \e alignment) 
         and the second is used to calculate the actual \e displacement.
         These are assigned in the reference (that you set with e.g.: REFERENCE=whatever.pdb) 
         frame which is, for this case, a simple pdb file containing the atoms
         on which you want to calculate the distance. 
         Note that in this pdb you need to set correctly the index of the atom STARTING FROM 1 (i.e. +1 shift respect 
         to VMD numbering) so that the program knows which atom needs to be compared with.
         The OCCUPANCY column (the first after the coordinates) set the \f$ w'  \f$ used for the center of mass calculation and the BETA column (the second
         after the Cartesian coordinates) set the \f$ w \f$ which is responsible of the calculation of the displacement. 
         Users can also use fractional values for beta and the occupancy values. We recommend you only do this when you really know what you are doing however as the results can be rather strange. 
         In PDB files the atomic coordinates and box lengths should be in Angstroms unless 
         you are working with natural units.  If you are working with natural units then the coordinates 
         should be in your natural length unit.  For more details on the PDB file format visit http://www.wwpdb.org/docs.html
         .
       - If you use TYPE=OPTIMAL you just the root mean square deviation after reset of the geometric center  
 and performing an optimal alignment that reads:
        \f[
         d(X,X_r) = \sqrt{ \sum_i \sum_\alpha^{x,y,z}  \frac{w_i}{\sum_j w_j}[ X_{i,\alpha}-com_\alpha(X)- \sum_\beta M(X,X',w')_{\alpha,\beta}({X'}_{i,\beta}-com_\beta(X')) ]^2 } 
        \f]
        where \f$ M(X,X',w') \f$ is the optimal alignment matrix which is calculated through the Kearsley \cite kearsley algorithm and the set of 
	weights used for the alignment of the center of mass is also used to perform the optimal alignment.
 	Note that the flexibility of setting only certain atoms to calculate the center of mass and optimal alignment which may be different from the one used in the \e displacemnt  calculation is useful whenever for example you want to calculate the motion of a ligand in a protein cavity and you want to clearly use the protein as reference system.  
        When this form of RMSD is used to calculate the secondary structure variables (\ref ALPHARMSD, \ref ANTIBETARMSD and \ref PARABETARMSD) all the atoms in the segment are assumed to be part of both the alignment and displacement sets.

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

PLUMED_REGISTER_ACTION(RMSD,"RMSD")

void RMSD::registerKeywords(Keywords& keys){
  Colvar::registerKeywords(keys);
  keys.add("compulsory","REFERENCE","a file in pdb format containing the reference structure and the atoms involved in the CV.");
  keys.add("compulsory","TYPE","SIMPLE","the manner in which RMSD alignment is performed.  Should be OPTIMAL or SIMPLE.");
  keys.addFlag("SQUARED",false," This should be setted if you want MSD instead of RMSD ");
}

RMSD::RMSD(const ActionOptions&ao):
PLUMED_COLVAR_INIT(ao),squared(false)
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

  rmsd = metricRegister().create<RMSDBase>(type,pdb);
  
  std::vector<AtomNumber> atoms;
  rmsd->getAtomRequests( atoms );
  rmsd->setNumberOfAtoms( atoms.size() );
  requestAtoms( atoms );

  log.printf("  reference from file %s\n",reference.c_str());
  log.printf("  which contains %d atoms\n",getNumberOfAtoms());
  log.printf("  method for alignment : %s \n",type.c_str() );
  if(squared)log.printf("  chosen to use SQUARED option for MSD instead of RMSD\n");
}

RMSD::~RMSD(){
  delete rmsd;
}


// calculator
void RMSD::calculate(){
  double r=rmsd->calculate( getPositions(), squared );

  setValue(r); 
  for(unsigned i=0;i<getNumberOfAtoms();i++) setAtomsDerivatives( i, rmsd->getAtomDerivative(i) );

  Tensor virial; plumed_dbg_assert( !rmsd->getVirial(virial) );
  setBoxDerivativesNoPbc();
}

}
}



