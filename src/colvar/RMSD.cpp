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
#include "RMSD.h"
#include "tools/PDB.h"
#include "core/ActionRegister.h"
#include "core/ActionSetup.h"
#include "core/ActionShortcut.h"
#include "core/PlumedMain.h"

namespace PLMD {
namespace colvar {

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

PLUMED_REGISTER_ACTION(RMSD,"RMSD_CALC")

void RMSD::registerRMSD(Keywords& keys ) {
  keys.add("compulsory","TYPE","SIMPLE","the manner in which RMSD alignment is performed.  Should be OPTIMAL or SIMPLE.");
  keys.addFlag("UNORMALIZED",false,"by default the mean sequare deviation or root mean square deviation is calculated.  If this option is given no averaging is done");
  keys.addFlag("SQUARED",false," This should be setted if you want MSD instead of RMSD ");
  keys.addFlag("DISPLACEMENT",false,"Calculate the vector of displacements instead of the length of this vector");
}

void RMSD::registerKeywords(Keywords& keys) {
  Action::registerKeywords(keys); ActionWithArguments::registerKeywords(keys); ActionWithValue::registerKeywords(keys); keys.use("ARG");
  keys.add("compulsory","ALIGN","1.0","the weights to use when aligning to the reference structure");
  keys.add("compulsory","DISPLACE","1.0","the weights to use when calculating the displacement from the reference structure");
  keys.addFlag("UNFIX",false,"this is used by adaptive path to make sure that the reference structures for the RMSD are kept up to date");
  RMSD::registerRMSD( keys );
  keys.addOutputComponent("disp","DISPLACEMENT","the vector of displacements for the atoms");
  keys.addOutputComponent("dist","DISPLACEMENT","the RMSD distance the atoms have moved");
}

void RMSD::createReferenceConfiguration( const std::string& lab, const std::string& input, PlumedMain& plumed, const unsigned number ) {
  FILE* fp=std::fopen(input.c_str(),"r"); bool do_read=true; std::vector<double> vals;
  if(!fp) plumed_merror("could not open reference file " + input); unsigned natoms=0, nframes=0;
   
  while ( do_read ) {
     PDB mypdb; do_read=mypdb.readFromFilepointer(fp,plumed.usingNaturalUnits(),0.1/plumed.getUnits().getLength());
     if( !do_read && nframes>0 ) break ;

     if( natoms==0 ) natoms = mypdb.getPositions().size();
     else if( mypdb.getPositions().size()!=natoms ) plumed_merror("mismatch between sizes of reference configurations");

     if( nframes+1==number || number==0 ) {
         std::vector<double> align( mypdb.getOccupancy() ); 
         double asum=0; for(unsigned i=0;i<align.size();++i) asum += align[i];
         if( asum>epsilon ) {
             double iasum = 1 / asum; for(unsigned i=0;i<align.size();++i) align[i] *= iasum;
         } else {
             double iasum = 1 / mypdb.size(); for(unsigned i=0;i<align.size();++i) align[i] = iasum;
         }
         Vector center; center.zero(); for(unsigned i=0;i<mypdb.getPositions().size();++i) center += align[i]*mypdb.getPositions()[i];       

         for(unsigned j=0; j<3; ++j) { 
             for(unsigned i=0; i<mypdb.getPositions().size(); ++i) vals.push_back( mypdb.getPositions()[i][j] - center[j] );
         }
     }
     nframes++;
  }
  std::fclose(fp); std::string rnum; plumed_assert( vals.size()>0 );
  Tools::convert( vals[0], rnum ); std::string valstr = " VALUES=" + rnum; 
  for(unsigned i=1; i<vals.size();++i) { Tools::convert( vals[i], rnum ); valstr += "," + rnum; }    
  if( number==0 && nframes>1 ) {
      std::string nc, nr; Tools::convert( nframes, nr ); Tools::convert( 3*natoms, nc ); 
      plumed.readInputLine( lab + ": CONSTANT_VALUE NROWS=" + nr + " NCOLS=" + nc + valstr );
  } else plumed.readInputLine( lab + ": CONSTANT_VALUE" + valstr );
}

void RMSD::createPosVector( const std::string& lab, const PDB& pdb, ActionShortcut* action ) { 
  bool nopbc; action->parseFlag("NOPBC",nopbc); std::vector<AtomNumber> anum( pdb.getAtomNumbers() ); 
  if( !nopbc ) {
      std::string num; Tools::convert( anum[0].serial(), num ); std::string wm_line = "WHOLEMOLECULES ENTITY0=" + num;
      for(unsigned i=1; i<anum.size(); ++i) { Tools::convert( anum[i].serial(), num ); wm_line += "," + num; }
      action->readInputLine( wm_line ); 
  }
  std::string num; Tools::convert( anum[0].serial(), num ); std::string pos_line = lab + "_pos: POSITION NOPBC ATOMS=" + num; 
  for(unsigned i=1; i<anum.size(); ++i) { Tools::convert( anum[i].serial(), num ); pos_line += "," + num; }
  action->readInputLine( pos_line );
  // Concatenate the three positions together
  action->readInputLine( lab + ": CONCATENATE ARG1=" + lab + "_pos.x ARG2=" + lab + "_pos.y ARG3=" + lab + "_pos.z");
}

RMSD::RMSD(const ActionOptions&ao):
  Action(ao),
  ActionWithArguments(ao),
  ActionWithValue(ao),
  firsttime(true),
  squared(false),
  displacement(false),
  multiple(false)
{
  if( getNumberOfArguments()!=2 ) error("there should be exactly two arguments for this action");
  // Check for shorcut 
  unsigned natoms, ntasks; bool unorm=false; parseFlag("UNORMALIZED",unorm); norm_weights=!unorm;
  if( getPntrToArgument(0)->getRank()==1 && getPntrToArgument(1)->getRank()==1 ) {
      natoms = getPntrToArgument(1)->getShape()[0] / 3; ntasks=1; myrmsd.resize(1);
  } else if( getPntrToArgument(0)->getRank()==2 ) {
      if( getPntrToArgument(1)->getRank()!=1 ) error("if first argument is a matrix second argument should be a vector");
      natoms = getPntrToArgument(1)->getShape()[0] / 3; multiple=true; ntasks = getPntrToArgument(0)->getShape()[0]; myrmsd.resize(1);
      if( getPntrToArgument(0)->getShape()[1]!=3*natoms ) error("mismatch between numbers of pos and reference");
  } else {
      if( getPntrToArgument(1)->getRank()!=2 ) error("if first argument is a vector second argument should be a matrix");
      if( getPntrToArgument(0)->getRank()!=1 ) error("if first argument is a matrix second argument should be vector");
      natoms = getPntrToArgument(0)->getShape()[0] / 3; multiple=true; ntasks = getPntrToArgument(1)->getShape()[0]; myrmsd.resize(ntasks);
      if( getPntrToArgument(1)->getShape()[1]!=3*natoms ) error("mismatch between numbers of in pos and reference");
  }
  // Request the arguments
  requestArguments( getArguments(), false ); 
  align.resize( natoms ); parseVector("ALIGN",align);
  displace.resize( natoms ); parseVector("DISPLACE",displace);

  type.assign("SIMPLE"); parse("TYPE",type);
  parseFlag("SQUARED",squared); parseFlag("DISPLACEMENT",displacement);

  double wa=0, wd=0; sqrtdisplace.resize( displace.size() );
  for(unsigned i=0; i<align.size(); ++i) { wa+=align[i]; wd+=displace[i]; }

  if( wa>epsilon ) { 
      double iwa = 1. / wa; 
      for(unsigned i=0; i<align.size(); ++i) align[i] *= iwa; 
  } else {
      double iwa = 1. / natoms;
      for(unsigned i=0; i<align.size(); ++i) align[i] = iwa; 
  } 
  if( wd>epsilon ) {
      if( unorm ) { wd = 1; } double iwd = 1. / wd; 
      for(unsigned i=0; i<align.size(); ++i) displace[i] *= iwd; 
  } else {
      double iwd = 1. / natoms;
      for(unsigned i=0; i<align.size(); ++i) displace[i] = iwd;
  }
  for(unsigned i=0; i<align.size(); ++i) sqrtdisplace[i] = sqrt(displace[i]);
  forcesToApply.resize( 3*natoms ); 

  if( displacement ) {
     std::vector<unsigned> shape0(2), shape1(1); shape1[0] = shape0[0] = ntasks; shape0[1] = getPntrToArgument(0)->getShape()[0];
     addComponent( "disp", shape0 ); componentIsNotPeriodic("disp");
     if( ntasks==1 ) addComponentWithDerivatives( "dist" ); else addComponent( "dist", shape1 ); 
     componentIsNotPeriodic("dist");
  } else {
     if( ntasks==1 ) addValueWithDerivatives(); 
     else { std::vector<unsigned> shape(1); shape[0]=ntasks; addValue( shape ); }
     setNotPeriodic(); 
  }

  // Print information to screen
  if( ntasks==1 ) log.printf("  calculating RMSD distance between two sets of %d atoms in vectors %s and %s\n", natoms, getPntrToArgument(1)->getName().c_str(), getPntrToArgument(0)->getName().c_str() );
  else if( getPntrToArgument(1)->getRank()==2 ) log.printf("  calculating RMSD distance of %d sets of atom positions in matrix with label %s from the %d atoms positions in vector with label %s \n", ntasks, getPntrToArgument(1)->getName().c_str(), natoms, getPntrToArgument(0)->getName().c_str()  );
  else log.printf("  calculating RMSD distance of %d atom positions in vector with label %s from %d sets of atom positions in matrix with label %s \n", natoms, getPntrToArgument(1)->getName().c_str(), ntasks, getPntrToArgument(0)->getName().c_str() );

  log.printf("  method for alignment : %s \n",type.c_str() );
  if(squared)log.printf("  chosen to use SQUARED option for MSD instead of RMSD\n");
  else      log.printf("  using periodic boundary conditions\n");
}

void RMSD::setReferenceConfigurations() {
  unsigned natoms = getPntrToArgument(1)->getShape()[0] / 3; 
  if( getPntrToArgument(1)->getRank()==2 ) natoms = getPntrToArgument(1)->getShape()[1] / 3;
  Vector center; std::vector<Vector> pos( natoms );
  for(unsigned jconf=0; jconf<myrmsd.size(); ++jconf) {
      center.zero();
      for(unsigned i=0; i<pos.size(); ++i) { 
          for(unsigned j=0; j<3; ++j) pos[i][j] = getPntrToArgument(1)->get( (3*jconf+j)*pos.size() + i ); 
          center+=pos[i]*align[i];
      }
      for(unsigned i=0; i<pos.size(); ++i) pos[i] -= center;
      myrmsd[jconf].clear(); myrmsd[jconf].set(align,displace,pos,type,true,norm_weights); 
  }
}

// calculator
void RMSD::calculate() {
  // Align reference configuration and set rmsd data
  if( firsttime || !getPntrToArgument(1)->isConstant() ) { setReferenceConfigurations(); firsttime=false; }
  // Now calculate all the RMSD values
  runAllTasks();
}

bool RMSD::performTask( const std::string& controller, const unsigned& index1, const unsigned& index2, MultiValue& myvals ) const {
  // Do not perform the loop here with a loop over other matrix elements
  if( controller!=getLabel() ) return false;

  unsigned jarg = index2 - getPntrToOutput(0)->getShape()[0], natoms = getPntrToArgument(0)->getShape()[0] / 3; 
  unsigned icomp = std::floor( jarg / natoms ); unsigned iatom = jarg - icomp*natoms;
  unsigned ostrn = getPntrToOutput(0)->getPositionInStream();   
  std::vector<Vector>& pos( myvals.getFirstAtomVector() ); myvals.addValue( ostrn, pos[iatom][icomp] );
  if( !doNotCalculateDerivatives() ) { myvals.addDerivative( ostrn, jarg, 1.0 ); myvals.updateIndex( ostrn, jarg ); }

  return true;
}

void RMSD::performTask( const unsigned& task_index, MultiValue& myvals ) const {
  // Get the index of the rmsd we are calculating
  unsigned rmsdno=0, structno=0, natoms = getPntrToArgument(0)->getShape()[0] / 3;
  if( getPntrToArgument(0)->getRank()==2 ) natoms = getPntrToArgument(0)->getShape()[1] / 3;
 
  if( multiple && myrmsd.size()>1 ) rmsdno=task_index;
  else if( multiple ) structno=task_index;  

  // Retrieve instantaneous configuration
  std::vector<Vector>& pos( myvals.getFirstAtomVector() ); std::vector<Vector>& der( myvals.getSecondAtomVector() );
  if( pos.size()!=natoms ) pos.resize( natoms ); if( der.size()!=natoms ) der.resize( natoms ); 
  for(unsigned i=0;i<pos.size();++i) {
      for(unsigned j=0; j<3; ++j) pos[i][j] = getPntrToArgument(0)->get( (3*structno+j)*natoms + i );
  }

  // Calculate RMSD distance
  double r; unsigned ostrn;
  if( displacement ) {
      // Calculate RMSD displacement 
      std::vector<Vector> direction( natoms ); Value* dval = getPntrToOutput(0);
      if(type=="SIMPLE") {
         r = myrmsd[rmsdno].simpleAlignment( align, displace, pos, myrmsd[rmsdno].getReference(), der, direction, squared );
         // Notice that we can adjust the forces here because we are parallelilising the apply loop over the RMSD values we are calculating
         if( dval->forcesWereAdded() ) {
             Vector comforce; comforce.zero();
             for(unsigned i=0; i<natoms; i++) {
                 for(unsigned k=0; k<3; ++k) comforce[k] += align[i]*dval->getForce( (task_index*3+k)*natoms + i);
             } 
             for(unsigned i=0; i<natoms; i++) {
                 for(unsigned k=0; k<3; ++k) dval->addForce( (task_index*3+k)*natoms + i, -comforce[k] );
             }
         }
      } else {
         Tensor rot; Matrix<std::vector<Vector> > DRotDPos(3,3); std::vector<Vector> centeredpos( natoms ), centeredreference( natoms );
         r = myrmsd[rmsdno].calc_PCAelements( pos, der, rot, DRotDPos, direction, centeredpos, centeredreference, squared );
         for(unsigned i=0;i<direction.size();++i) direction[i] = sqrtdisplace[i]*( direction[i] - myrmsd[rmsdno].getReference()[i] );
         // Notice that we can adjust the forces here because we are parallelilising the apply loop over the RMSD values we are calculating
         if( dval->forcesWereAdded() ) {
             Tensor trot=rot.transpose(); double prefactor = 1 / static_cast<double>( natoms ); Vector v1; v1.zero();
             for(unsigned n=0; n<natoms; n++) { 
                  Vector ff; for(unsigned k=0; k<3; ++k ) ff[k] = dval->getForce( (task_index*3+k)*natoms + n ); 
                  v1+=prefactor*matmul(trot,ff);
             }
             // Notice that we use centreredreference here to accumulate the true forces
             for(unsigned n=0; n<natoms; n++) { 
                  Vector ff; for(unsigned k=0; k<3; ++k ) ff[k] = dval->getForce( (task_index*3+k)*natoms + n ); 
                  centeredreference[n] = sqrtdisplace[n]*( matmul(trot,ff) - v1 );
             }
             for(unsigned a=0; a<3; a++) {
                 for(unsigned b=0; b<3; b++) {
                     for(unsigned i=0; i<natoms; i++) {
                         double tmp1=0.; for(unsigned m=0; m<natoms; m++) tmp1+=centeredpos[m][b]*dval->getForce( (task_index*3+a)*natoms + m );
                         centeredreference[i] += sqrtdisplace[i]*tmp1*DRotDPos[a][b][i];
                     }
                 }
             }
             // Now subtract the current force and add on the true force
             for(unsigned n=0; n<natoms; n++) {
                 for(unsigned k=0; k<3; ++k) dval->addForce( (task_index*3+k)*natoms + n, centeredreference[n][k]-dval->getForce( (task_index*3+k)*natoms + n ) );
             }
         }
      }
      unsigned base = getPntrToOutput(0)->getShape()[0];
      for(unsigned j=0; j<3; ++j) {
          for(unsigned i=0; i<pos.size(); ++i) {
              pos[i][j] = direction[i][j]; 
              // This ensures that the matrix element is gathered
              runTask( getLabel(), task_index, base, myvals ); 
              // Now clear only elements that are not accumulated over whole row 
              clearMatrixElements( myvals ); base++;
          } 
      } 
      // Set the value that we are outputting on
      ostrn = getPntrToOutput(1)->getPositionInStream();
  } else {
      r=myrmsd[rmsdno].calculate( pos, der, squared ); ostrn = getPntrToOutput(0)->getPositionInStream();
  }
  myvals.setValue( ostrn, r ); if( doNotCalculateDerivatives() ) return; 

  for(unsigned i=0; i<natoms; i++){ 
      for(unsigned j=0; j<3; ++j ) { myvals.addDerivative( ostrn, j*natoms+i, der[i][j] ); myvals.updateIndex( ostrn, j*natoms+i ); } 
  }
}

void RMSD::apply() {
  if( doNotCalculateDerivatives() ) return;
  std::fill(forcesToApply.begin(),forcesToApply.end(),0); 
  //  We need a clever trick to get the forces in this special case because of a curiosity in ActionWithValue 
  if( displacement ) {
      if( getPntrToComponent(1)->getRank()==0 && getPntrToComponent(1)->forcesWereAdded() ) getPntrToComponent(1)->applyForce( forcesToApply ); 
  }
  unsigned mm=0; if( getForcesFromValues( forcesToApply ) ) setForcesOnArguments( 0, forcesToApply, mm );
}

}
}



