/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2020 The plumed team
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
#include "Colvar.h"
#include "core/PlumedMain.h"
#include "ActionRegister.h"
#include "tools/PDB.h"
#include "reference/RMSDBase.h"
#include "reference/MetricRegister.h"
#include "core/Atoms.h"

namespace PLMD {
namespace colvar {

class RMSD : public Colvar {

  MultiValue myvals;
  ReferenceValuePack mypack;
  std::unique_ptr<PLMD::RMSDBase> rmsd;
  bool squared;
  bool nopbc;

public:
  explicit RMSD(const ActionOptions&);
  void calculate() override;
  static void registerKeywords(Keywords& keys);
};

//+PLUMEDOC DCOLVAR RMSD
/*
Calculate the RMSD with respect to a reference structure.

The aim with this colvar it to calculate something like:

\f[
d(X,X') = \vert X-X' \vert
\f]

where \f$ X \f$ is the instantaneous position of all the atoms in the system and
\f$ X' \f$ is the positions of the atoms in some reference structure provided as input.
\f$ d(X,X') \f$ thus measures the distance all the atoms have moved away from this reference configuration.
Oftentimes, it is only the internal motions of the structure - i.e. not the translations of the center of
mass or the rotations of the reference frame - that are interesting.  Hence, when calculating the
the root-mean-square deviation between the atoms in two configurations
you must first superimpose the two structures in some way. At present PLUMED provides two distinct ways
of performing this superposition.  The first method is applied when you use TYPE=SIMPLE in the input
line.  This instruction tells PLUMED that the root mean square deviation is to be calculated after the
positions of the geometric centers in the reference and instantaneous configurations are aligned.  In
other words \f$d(X,x')\f$ is to be calculated using:

\f[
 d(X,X') = \sqrt{ \sum_i \sum_\alpha^{x,y,z}  \frac{w_i}{\sum_j w_j}( X_{i,\alpha}-com_\alpha(X)-{X'}_{i,\alpha}+com_\alpha(X') )^2 }
\f]
with
\f[
com_\alpha(X)= \sum_i  \frac{w'_{i}}{\sum_j w'_j}X_{i,\alpha}
\f]
and
\f[
com_\alpha(X')= \sum_i  \frac{w'_{i}}{\sum_j w'_j}X'_{i,\alpha}
\f]
Obviously, \f$ com_\alpha(X) \f$ and  \f$ com_\alpha(X') \f$  represent the positions of the center of mass in the reference
and instantaneous configurations if the weights $w'$ are set equal to the atomic masses.  If the weights are all set equal to
one, however, \f$com_\alpha(X) \f$ and  \f$ com_\alpha(X') \f$ are the positions of the geometric centers.
Notice that there are sets of weights:  \f$ w' \f$ and  \f$ w \f$. The first is used to calculate the position of the center of mass
(so it determines how the atoms are \e aligned).  Meanwhile, the second is used when calculating how far the atoms have actually been
\e displaced.  These weights are assigned in the reference configuration that you provide as input (i.e. the appear in the input file
to this action that you set using REFERENCE=whatever.pdb).  This input reference configuration consists of a simple pdb file
containing the set of atoms for which you want to calculate the RMSD displacement and their positions in the reference configuration.
It is important to note that the indices in this pdb need to be set correctly.  The indices in this file determine the indices of the
instantaneous atomic positions that are used by PLUMED when calculating this colvar.  As such if you want to calculate the RMSD distance
moved by the first, fourth, sixth and twenty eighth atoms in the MD codes input file then the indices of the corresponding reference positions in this pdb
file should be set equal to 1, 4, 6 and 28.

The pdb input file should also contain the values of \f$w\f$ and \f$w'\f$. In particular, the OCCUPANCY column (the first column after the coordinates)
is used provides the values of \f$ w'\f$ that are used to calculate the position of the center of mass.  The BETA column (the second column
after the Cartesian coordinates) is used to provide the \f$ w \f$ values which are used in the the calculation of the displacement.
Please note that it is possible to use fractional values for beta and for the occupancy. However, we recommend you only do this when
you really know what you are doing however as the results can be rather strange.

In PDB files the atomic coordinates and box lengths should be in Angstroms unless
you are working with natural units.  If you are working with natural units then the coordinates
should be in your natural length unit.  For more details on the PDB file format visit http://www.wwpdb.org/docs.html.
Make sure your PDB file is correctly formatted as explained \ref pdbreader "in this page".

A different method is used to calculate the RMSD distance when you use TYPE=OPTIMAL on the input line.  In this case  the root mean square
deviation is calculated after the positions of geometric centers in the reference and instantaneous configurations are aligned AND after
an optimal alignment of the two frames is performed so that motion due to rotation of the reference frame between the two structures is
removed.  The equation for \f$d(X,X')\f$ in this case reads:

\f[
d(X,X') = \sqrt{ \sum_i \sum_\alpha^{x,y,z}  \frac{w_i}{\sum_j w_j}[ X_{i,\alpha}-com_\alpha(X)- \sum_\beta M(X,X',w')_{\alpha,\beta}({X'}_{i,\beta}-com_\beta(X')) ]^2 }
\f]

where \f$ M(X,X',w') \f$ is the optimal alignment matrix which is calculated using the Kearsley \cite kearsley algorithm.  Again different sets of
weights are used for the alignment (\f$w'\f$) and for the displacement calculations (\f$w\f$).
This gives a great deal of flexibility as it allows you to use a different sets of atoms (which may or may not overlap) for the alignment and displacement
parts of the calculation. This may be very useful when you want to calculate how a ligand moves about in a protein cavity as you can use the protein as a reference
system and do no alignment of the ligand.

(Note: when this form of RMSD is used to calculate the secondary structure variables (\ref ALPHARMSD, \ref ANTIBETARMSD and \ref PARABETARMSD
all the atoms in the segment are assumed to be part of both the alignment and displacement sets and all weights are set equal to one)

Please note that there are a number of other methods for calculating the distance between the instantaneous configuration and a reference configuration
that are available in plumed.  More information on these various methods can be found in the section of the manual on \ref dists.

When running with periodic boundary conditions, the atoms should be
in the proper periodic image. This is done automatically since PLUMED 2.5,
by considering the ordered list of atoms and rebuilding molecules using a procedure
that is equivalent to that done in \ref WHOLEMOLECULES . Notice that
rebuilding is local to this action. This is different from \ref WHOLEMOLECULES
which actually modifies the coordinates stored in PLUMED.

In case you want to recover the old behavior you should use the NOPBC flag.
In that case you need to take care that atoms are in the correct
periodic image.

\par Examples

The following tells plumed to calculate the RMSD distance between
the positions of the atoms in the reference file and their instantaneous
position.  The Kearsley algorithm is used so this is done optimally.

\plumedfile
RMSD REFERENCE=file.pdb TYPE=OPTIMAL
\endplumedfile

The reference configuration is specified in a pdb file that will have a format similar to the one shown below:

\auxfile{file.pdb}
ATOM      1  CL  ALA     1      -3.171   0.295   2.045  1.00  1.00
ATOM      5  CLP ALA     1      -1.819  -0.143   1.679  1.00  1.00
ATOM      6  OL  ALA     1      -1.177  -0.889   2.401  1.00  1.00
ATOM      7  NL  ALA     1      -1.313   0.341   0.529  1.00  1.00
ATOM      8  HL  ALA     1      -1.845   0.961  -0.011  1.00  1.00
END
\endauxfile

...

*/
//+ENDPLUMEDOC

PLUMED_REGISTER_ACTION(RMSD,"RMSD")

void RMSD::registerKeywords(Keywords& keys) {
  Colvar::registerKeywords(keys);
  keys.add("compulsory","REFERENCE","a file in pdb format containing the reference structure and the atoms involved in the CV.");
  keys.add("compulsory","TYPE","SIMPLE","the manner in which RMSD alignment is performed.  Should be OPTIMAL or SIMPLE.");
  keys.addFlag("SQUARED",false," This should be set if you want mean squared displacement instead of RMSD ");
}

RMSD::RMSD(const ActionOptions&ao):
  PLUMED_COLVAR_INIT(ao),
  myvals(1,0),
  mypack(0,0,myvals),
  squared(false),
  nopbc(false)
{
  std::string reference;
  parse("REFERENCE",reference);
  std::string type;
  type.assign("SIMPLE");
  parse("TYPE",type);
  parseFlag("SQUARED",squared);
  parseFlag("NOPBC",nopbc);

  checkRead();


  addValueWithDerivatives(); setNotPeriodic();
  PDB pdb;

  // read everything in ang and transform to nm if we are not in natural units
  if( !pdb.read(reference,plumed.getAtoms().usingNaturalUnits(),0.1/atoms.getUnits().getLength()) )
    error("missing input file " + reference );

  rmsd=metricRegister().create<RMSDBase>(type,pdb);

  std::vector<AtomNumber> atoms;
  rmsd->getAtomRequests( atoms );
//  rmsd->setNumberOfAtoms( atoms.size() );
  requestAtoms( atoms );

  // Setup the derivative pack
  myvals.resize( 1, 3*atoms.size()+9 ); mypack.resize( 0, atoms.size() );
  for(unsigned i=0; i<atoms.size(); ++i) mypack.setAtomIndex( i, i );

  log.printf("  reference from file %s\n",reference.c_str());
  log.printf("  which contains %d atoms\n",getNumberOfAtoms());
  log.printf("  with indices : ");
  for(unsigned i=0; i<atoms.size(); ++i) {
    if(i%25==0) log<<"\n";
    log.printf("%d ",atoms[i].serial());
  }
  log.printf("\n");
  log.printf("  method for alignment : %s \n",type.c_str() );
  if(squared)log.printf("  chosen to use SQUARED option for MSD instead of RMSD\n");
  if(nopbc) log.printf("  without periodic boundary conditions\n");
  else      log.printf("  using periodic boundary conditions\n");
}


// calculator
void RMSD::calculate() {
  if(!nopbc) makeWhole();
  double r=rmsd->calculate( getPositions(), mypack, squared );

  setValue(r);
  for(unsigned i=0; i<getNumberOfAtoms(); i++) setAtomsDerivatives( i, mypack.getAtomDerivative(i) );

  Tensor virial; plumed_dbg_assert( !mypack.virialWasSet() );
  setBoxDerivativesNoPbc();
}

}
}



