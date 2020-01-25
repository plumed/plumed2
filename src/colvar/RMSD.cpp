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
#include "core/Colvar.h"
#include "core/PlumedMain.h"
#include "core/ActionRegister.h"
#include "tools/RMSD.h"
#include "core/Atoms.h"
#include "core/ActionSetup.h"
#include "tools/PDB.h"


using namespace std;

namespace PLMD {
namespace colvar {

class RMSD : public Colvar {
private:
  bool fixed_reference;
  Tensor rot;
  Matrix<std::vector<Vector> > DRotDPos;
  std::vector<Vector> pos, der, direction, centeredpos, centeredreference;
  bool squared;
  bool nopbc;
  bool displacement;
  bool norm_weights;
  std::string type;
  std::vector<double> align,displace,sqrtdisplace;
  PLMD::RMSD myrmsd;
  std::vector<Vector> forcesToApply;
  void setReferenceConfiguration();
public:
  explicit RMSD(const ActionOptions&);
  virtual void calculate() override;
  void apply() override;
  static void registerKeywords(Keywords& keys);
};


using namespace std;

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
  keys.add("optional","REFERENCE","a file in pdb format containing the reference structure and the atoms involved in the CV");
  keys.add("atoms","REFERENCE_ATOMS","the atom numbers for the reference configuration");
  keys.add("atoms","ATOMS","the atom numbers that you would like to consider");
  keys.add("compulsory","ALIGN","1.0","the weights to use when aligning to the reference structure");
  keys.add("compulsory","DISPLACE","1.0","the weights to use when calculating the displacement from the reference structure");
  keys.add("compulsory","TYPE","SIMPLE","the manner in which RMSD alignment is performed.  Should be OPTIMAL or SIMPLE.");
  keys.addFlag("UNORMALIZED",false,"by default the mean sequare deviation or root mean square deviation is calculated.  If this option is given no averaging is done");
  keys.addFlag("SQUARED",false," This should be setted if you want MSD instead of RMSD ");
  keys.addFlag("DISPLACEMENT",false,"Calculate the vector of displacements instead of the length of this vector");
  keys.addOutputComponent("disp","DISPLACEMENT","the vector of displacements for the atoms");
  keys.addOutputComponent("dist","DISPLACEMENT","the RMSD distance the atoms have moved");
}

RMSD::RMSD(const ActionOptions&ao):
  PLUMED_COLVAR_INIT(ao),
  fixed_reference(true),
  DRotDPos(3,3),
  squared(false),
  nopbc(false),
  displacement(false)
{
  // Check for shorcut 
  std::vector<AtomNumber> atoms_ref, atoms_conf;
  std::string reference; parse("REFERENCE",reference);
  bool unorm=false; parseFlag("UNORMALIZED",unorm); norm_weights=!unorm;
  if( reference!="") {
      // Create the input reference position
      readInputLine( getLabel() + "_ref: READ_CONFIG REFERENCE=" + reference );
      // Now create the input for the real RMSD object
      std::string rmsd_line = getLabel() + ": RMSD REFERENCE_ATOMS=" + getLabel() + "_ref";
      // Read the reference pdb file
      PDB pdb; 
      if( !pdb.read(reference,plumed.getAtoms().usingNaturalUnits(),0.1/atoms.getUnits().getLength()) ) plumed_merror("missing file " + reference );
      // Get the reference atoms
      std::vector<std::string> refname(1); refname[0]=getLabel() + "_ref"; interpretAtomList( refname, atoms_ref ); 
      // Get the atom numbers
      atoms_conf = pdb.getAtomNumbers(); align = pdb.getOccupancy(); displace = pdb.getBeta();
  } else {
      parseAtomList("REFERENCE_ATOMS",atoms_ref); parseAtomList("ATOMS",atoms_conf);
      align.resize( atoms_ref.size() ); parseVector("ALIGN",align);
      displace.resize( atoms_ref.size() ); parseVector("DISPLACE",displace);
  }

  type.assign("SIMPLE"); parse("TYPE",type);
  parseFlag("SQUARED",squared); parseFlag("NOPBC",nopbc); parseFlag("DISPLACEMENT",displacement);

  if( atoms_ref.size()!=atoms_conf.size() ) error("size mismatch between reference atoms and atoms involved");
  double wa=0, wd=0; sqrtdisplace.resize( displace.size() );
  for(unsigned i=0; i<align.size(); ++i) { wa+=align[i]; wd+=displace[i]; }
  if( unorm ) { wd = 1; }
  for(unsigned i=0; i<align.size(); ++i){ align[i] /= wa; displace[i] /= wd; sqrtdisplace[i] = sqrt(displace[i]); }

  if( displacement ) {
     std::vector<unsigned> shape(1); shape[0] = 3*atoms_conf.size();
     addComponent( "disp", shape ); componentIsNotPeriodic("disp");
     addComponentWithDerivatives( "dist" ); componentIsNotPeriodic("dist");
     getPntrToComponent(0)->alwaysStoreValues();
     forcesToApply.resize( atoms_conf.size() );
  } else {
     addValueWithDerivatives(); setNotPeriodic();
  }

  std::vector<AtomNumber> myatoms( atoms_ref );
  for(unsigned i=0;i<atoms_conf.size();++i) myatoms.push_back( atoms_conf[i] );
  requestAtoms( myatoms ); pos.resize( atoms_ref.size() ); der.resize( atoms_ref.size() );
  direction.resize( atoms_ref.size() ); ; centeredpos.resize( atoms_ref.size() ); 
  centeredreference.resize( atoms_ref.size() );

  // Determine if the reference configuration is fixed
  for(unsigned i=0;i<atoms_conf.size();++i) {
      if( atoms.isVirtualAtom(atoms_ref[i]) ) {
          ActionSetup* as = dynamic_cast<ActionSetup*>( atoms.getVirtualAtomsAction(atoms_ref[i]) );
          if( !as ) fixed_reference=false;
      } else fixed_reference=false;
  }

  // Print information to screen
  log.printf("  calculating RMSD distance between two sets of %d atoms\n", getNumberOfAtoms() / 2);
  if( fixed_reference ) {
      // Check atoms in configuration we are measuring distance from are not in reference group
      for(unsigned i=0;i<atoms_conf.size();++i) {
          for(unsigned j=0;j<atoms_ref.size();++j) {
              if( atoms_conf[i]==atoms_ref[j] ) error("atom out of range");
          }
      }
      ActionSetup* as = dynamic_cast<ActionSetup*>( atoms.getVirtualAtomsAction(atoms_ref[0]) );
      log.printf("  reference configuration is fixed and was read in by action with label %s \n", as->getLabel().c_str() );
      retrieveAtoms(); setReferenceConfiguration();
  }

  log.printf("  method for alignment : %s \n",type.c_str() );
  if(squared)log.printf("  chosen to use SQUARED option for MSD instead of RMSD\n");
  if(nopbc) log.printf("  without periodic boundary conditions\n");
  else      log.printf("  using periodic boundary conditions\n");
}

void RMSD::setReferenceConfiguration() {
  if( !fixed_reference && !nopbc ) makeWhole( 0, pos.size() );
  for(unsigned i=0;i<pos.size();++i) pos[i] = getPosition(i);

  Vector center;
  for(unsigned i=0; i<pos.size(); ++i) center+=pos[i]*align[i];
  for(unsigned i=0; i<pos.size(); ++i) pos[i] -= center;
  myrmsd.clear(); myrmsd.set(align,displace,pos,type,true,norm_weights);
}

// calculator
void RMSD::calculate() {
  // Align reference configuration and set rmsd data
  if( !fixed_reference ) setReferenceConfiguration();
  // Make the molecule whole
  if( !nopbc ) makeWhole( pos.size(), getNumberOfAtoms() );
  // Retrieve instantaneous configuration
  for(unsigned i=0;i<pos.size();++i) pos[i] = getPosition(pos.size()+i);

  if( displacement ) {
      // Calculate RMSD displacement 
      double d; 
      if(type=="SIMPLE") { 
         d = myrmsd.simpleAlignment( align, displace, pos, myrmsd.getReference(), der, direction, squared );
      } else {
         d = myrmsd.calc_PCAelements( pos, der, rot, DRotDPos, direction, centeredpos, centeredreference, squared );
         for(unsigned i=0;i<direction.size();++i) direction[i] = sqrtdisplace[i]*( direction[i] - myrmsd.getReference()[i] );
      }
      Value* val=getPntrToComponent(0);
      for(unsigned i=0;i<pos.size();++i) {
          val->set( 3*i+0, direction[i][0] ); val->set( 3*i+1, direction[i][1] ); val->set( 3*i+2, direction[i][2] );
      }
      // Set the value of the derivatives
      Value* vv = getPntrToComponent(1); vv->set(d); unsigned n=0; Tensor virial; virial.zero();
      for(unsigned i=pos.size(); i<getNumberOfAtoms(); i++){ 
          setAtomsDerivatives( vv, i, der[n] ); virial -= Tensor( pos[n], der[n] ); n++; 
      } 
      setBoxDerivatives( vv, virial );
  } else {
      // Calculate RMSD distance
      double r=myrmsd.calculate( pos, der, squared );
      // Set value and derivatives
      setValue(r); unsigned n=0; Tensor virial; virial.zero();
      for(unsigned i=pos.size(); i<getNumberOfAtoms(); i++){ 
          setAtomsDerivatives( i, der[n] ); virial -= Tensor( pos[n], der[n] ); n++; 
      }
      // And finish by working out virial
      setBoxDerivatives( virial );
  }
}

void RMSD::apply() {
  // This ensures forces are applied to dist component
  Colvar::apply();
  // If we have not calculated a displacement we are done
  if( !displacement ) return ;
  // Check if forces have been applied to displacement
  if( !getPntrToOutput(0)->forcesWereAdded() ) return ;

  Value* dval=getPntrToComponent(0); 
  if( type=="SIMPLE" ) {
      Vector comforce; comforce[0];
      for(unsigned i=0; i<pos.size(); i++) {
          for(unsigned k=0; k<3; ++k) comforce[k] += align[i]*dval->getForce(3*i+k);
      } 
      for(unsigned i=0; i<pos.size(); i++) {
          forcesToApply[i][0] = dval->getForce( 3*i+0 ) - comforce[0];
          forcesToApply[i][1] = dval->getForce( 3*i+1 ) - comforce[1];
          forcesToApply[i][2] = dval->getForce( 3*i+2 ) - comforce[2];
      }
  } else {
      Tensor trot=rot.transpose(); double prefactor = 1 / static_cast<double>( pos.size() ); Vector v1; v1.zero();
      for(unsigned n=0; n<pos.size(); n++) { 
          Vector ff; ff[0] = dval->getForce( 3*n+0 ); ff[1] = dval->getForce( 3*n+1 ); ff[2] = dval->getForce( 3*n+2 );
          v1+=prefactor*matmul(trot,ff);
      }
      for(unsigned i=0; i<pos.size(); i++) { 
          Vector ff; ff[0] = dval->getForce( 3*i+0 ); ff[1] = dval->getForce( 3*i+1 ); ff[2] = dval->getForce( 3*i+2 );
          forcesToApply[i] = sqrtdisplace[i]*( matmul(trot,ff) - v1 );
      }
      for(unsigned a=0; a<3; a++) {
        for(unsigned b=0; b<3; b++) {
          for(unsigned i=0; i<pos.size(); i++) {
            double tmp1=0.; for(unsigned m=0; m<pos.size(); m++) tmp1+=centeredpos[m][b]*dval->getForce( 3*m+a );
            forcesToApply[i] += sqrtdisplace[i]*tmp1*DRotDPos[a][b][i];
          }
        }
      }
  }
  // Make the structure whole
  if( !nopbc ) makeWhole( pos.size(), getNumberOfAtoms() );
  // Retrieve instantaneous configuration
  for(unsigned i=0;i<pos.size();++i) pos[i] = getPosition(pos.size()+i);

  Tensor& v(modifyVirial()); std::vector<Vector>& f(modifyForces()); unsigned n=pos.size();
  for(unsigned i=0; i<pos.size(); i++) {
      f[n][0] += forcesToApply[i][0]; f[n][1] += forcesToApply[i][1]; f[n][2] += forcesToApply[i][2]; n++; 
      v -= Tensor( pos[i], forcesToApply[i] );
  }
}

}
}



