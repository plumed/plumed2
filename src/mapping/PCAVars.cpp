/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2018 The plumed team
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
#include "core/ActionShortcut.h"
#include "core/ActionRegister.h"
#include "core/PlumedMain.h"
#include "core/ActionWithArguments.h"
#include "tools/PDB.h"

using namespace std;

namespace PLMD {
namespace mapping {

//+PLUMEDOC COLVAR PCAVARS
/*
Projection on principal component eigenvectors or other high dimensional linear subspace

The collective variables described in \ref dists allow one to calculate the distance between the
instantaneous structure adopted by the system and some high-dimensional, reference configuration.  The
problem with doing this is that, as one gets further and further from the reference configuration, the
distance from it becomes a progressively poorer and poorer collective variable.  This happens because
the ``number" of structures at a distance \f$d\f$ from a reference configuration is proportional to \f$d^N\f$ in
an \f$N\f$ dimensional space.  Consequently, when \f$d\f$ is small the distance from the reference configuration
may well be a good collective variable.  However, when \f$d\f$ is large it is unlikely that the distance from the reference
structure is a good CV.  When the distance is large there will almost certainly be markedly different
configuration that have the same CV value and hence barriers in transverse degrees of
freedom. 
  
For these reasons dimensionality reduction is often employed so a projection \f$\mathbf{s}\f$ of a high-dimensional configuration
\f$\mathbf{X}\f$ in a lower dimensionality space using a function:

\f[
\mathbf{s} = F(\mathbf{X}-\mathbf{X}^{ref})
\f]

where here we have introduced some high-dimensional reference configuration \f$\mathbf{X}^{ref}\f$.  By far the simplest way to
do this is to use some linear operator for \f$F\f$.  That is to say we find a low-dimensional projection
by rotating the basis vectors using some linear algebra:

\f[
\mathbf{s}_i = \sum_k A_{ik} ( X_{k} - X_{k}^{ref} )
\f]

Here \f$A\f$ is a \f$d\f$ by \f$D\f$ matrix where \f$D\f$ is the dimensionality of the high dimensional space and \f$d\f$ is
the dimensionality of the lower dimensional subspace.  In plumed when this kind of projection you can use the majority
of the metrics detailed on \ref dists to calculate the displacement, \f$\mathbf{X}-\mathbf{X}^{ref}\f$, from the reference configuration.
The matrix \f$A\f$ can be found by various means including principal component analysis and normal mode analysis.  In both these methods the
rows of \f$A\f$ would be the principle eigenvectors of a square matrix.  For PCA the covariance while for normal modes the Hessian.

\bug It is not possible to use the \ref DRMSD metric with this variable.  You can get around this by listing the set of distances you wish to calculate for your DRMSD in the plumed file explicitly and using the EUCLIDEAN metric.  MAHALONOBIS and NORM-EUCLIDEAN also do not work with this variable but using these options makes little sense when projecting on a linear subspace.
  
\par Examples

The following input calculates a projection on a linear subspace where the displacements
from the reference configuration are calculated using the OPTIMAL metric.  Consequently,
both translation of the center of mass of the atoms and rotation of the reference
frame are removed from these displacements.  The matrix \f$A\f$ and the reference
configuration \f$R^{ref}\f$ are specified in the pdb input file reference.pdb and the
value of all projections (and the residual) are output to a file called colvar2.

\plumedfile
PCAVARS REFERENCE=reference.pdb TYPE=OPTIMAL LABEL=pca2
PRINT ARG=pca2.* FILE=colvar2
\endplumedfile

The reference configurations can be specified using a pdb file.  The first configuration that you provide is the reference configuration,
which is referred to in the above as \f$X^{ref}\f$ subsequent configurations give the directions of row vectors that are contained in
the matrix \f$A\f$ above.  These directions are specified by giving a second configuration that describes the components of \f$A\f$ explicitly.  

\verbatim
ATOM      2  CH3 ACE     1      12.932 -14.718  -6.016  1.00  1.00
ATOM      5  C   ACE     1      21.312  -9.928  -5.946  1.00  1.00
ATOM      9  CA  ALA     2      19.462 -11.088  -8.986  1.00  1.00
ATOM     13  HB2 ALA     2      21.112 -10.688 -12.476  1.00  1.00
ATOM     15  C   ALA     2      19.422   7.978 -14.536  1.00  1.00
ATOM     20 HH31 NME     3      20.122  -9.928 -17.746  1.00  1.00
ATOM     21 HH32 NME     3      18.572 -13.148 -16.346  1.00  1.00
END
REMARK
ATOM      2  CH3 ACE     1      0.1414  0.3334 -0.0302  1.00  0.00
ATOM      5  C   ACE     1      0.0893 -0.1095 -0.1434  1.00  0.00
ATOM      9  CA  ALA     2      0.0207 -0.321   0.0321  1.00  0.00
ATOM     13  HB2 ALA     2      0.0317 -0.6085  0.0783  1.00  0.00
ATOM     15  C   ALA     2      0.1282 -0.4792  0.0797  1.00  0.00
ATOM     20 HH31 NME     3      0.0053 -0.465   0.0309  1.00  0.00
ATOM     21 HH32 NME     3     -0.1019 -0.4261 -0.0082  1.00  0.00
END
\endverbatim

Notice that the PCAVARS command is a shortcut.  If you look at how the shortcut in the above input is expanded you should be able to see how the command works 
by calculating the RMSD displacement between the instantaneous and reference configuration and how those displacements are then projected on the eigenvector that 
was specified in the second frame of the pdb input above. Understanding the expanded version of this shortcut command allows you to recognise that you can project 
the displacement vector on any arbitrary vector.  For example in the input below two reference structures are provided in the pdb file. PLUMED calculates the RMSD
distance between these two reference configurations during setup and sets up a unit vector called eig that points along the director connecting the two RMSD structure.
During the calculation the vector connecting the instantaneous configuration and the first of the two reference configurations in computed.  This vector is then projected 
on the unit vector connecting the two initial structures:

\plumedfile
# Read in two reference configuratioms from PDB file
ref1: PDB2CONSTANT REFERENCE=two-frames.pdb NUMBER=1
ref2: PDB2CONSTANT REFERENCE=two-frames.pdb NUMBER=2
# Calculate the displacement vector that takes you from ref1 to ref2 
eigu: RMSD_DISPLACEMENT_VECTOR ARG=ref1,ref2 TYPE=OPTIMAL
# Normalise the reference vector
eigu2: CUSTOM ARG=eigu.disp FUNC=x*x PERIODIC=NO
eign2: SUM ARG=eigu2 PERIODIC=NO
eig: CUSTOM ARG=eigu.disp,eign2 FUNC=x/sqrt(y) PERIODIC=NO
# Everything prior to this point is only done in setup.  From here onwards we have the commands that are done during the main calculate loop.

# Calculate the RMSD displacement between the instaneous structure and the first reference structure
rmsd: RMSD REFERENCE=two-frames.pdb NUMBER=1 TYPE=OPTIMAL DISPLACEMENT SQUARED
# Project the displacement computed above on the director of the vector that connects reference structure ref1 to refeference structure ref2
prod: CUSTOM ARG=rmsd.disp,eig FUNC=x*y PERIODIC=NO
pca: SUM ARG=prod PERIODIC=NO

# Print the final CV to a file
PRINT ARG=pca FILE=colvar
\endplumedfile

You also project vectors of differences of arguments on reference vectors.  For example, the input below can be used to look at the projection 
of the vector connecting the instantanous configuraiton to a reference point in CV on a reference vector that has been specified in the PDB file.

\plumedfile
d1: DISTANCE ATOMS=1,2
d2: DISTANCE ATOMS=3,4
d3: DISTANCE ATOMS=5,6
pca: PCAVARS ARG=d1,d2,d3 REFERENCE=epath.pdb
PRINT ARG=pca_eig-1,pca_residual FILE=colvar
\endplumedfile

The pdb input file for this calculation might look something like this:

\verbatim
REMARK d1=0.1221 d2=0.0979 d3=0.1079
END
REMARK d1=0.078811 d2=-0.945732 d3=-0.315244
END
\endverbatim

The first set of argument values in this input file are the reference values for the arguments.  The second any subsquent sets of arguments give the 
coefficients that should be used when constructing linear combinations.

Notice, lastly, that you can also use a combination of argument values and atomic positions when specifying the reference configuration and the reference
directions.  If you are doing something this complicated, however, you are perhaps better working with the PLUMED input directly rather than this shortcut
as you will need to take special measures to ensure that all your CVs are in the same units.

*/
//+ENDPLUMEDOC

class PCAVars : public ActionShortcut {
public:
  static void registerKeywords(Keywords& keys);
  explicit PCAVars(const ActionOptions&);
};


PLUMED_REGISTER_ACTION(PCAVars,"PCAVARS")

void PCAVars::registerKeywords( Keywords& keys ) {
  ActionShortcut::registerKeywords( keys );
  keys.add("compulsory","REFERENCE","a pdb file containing the set of reference configurations");
  keys.add("compulsory","TYPE","OPTIMAL-FAST","the manner in which distances are calculated. More information on the different "
           "metrics that are available in PLUMED can be found in the section of the manual on "
           "\\ref dists");
  keys.add("optional","ARG","if there are arguments to be used specify them here");
  keys.addFlag("NOPBC",false,"do not use periodic boundary conditions when computing this quantity");
}

PCAVars::PCAVars( const ActionOptions& ao ):
  Action(ao),
  ActionShortcut(ao)
{
  std::string reference; parse("REFERENCE",reference); 
  // Create the object that holds the atomic positions by reading the first frame
  FILE* fp=std::fopen(reference.c_str(),"r"); PDB pdb; if(!fp) error("could not open reference file " + reference );
  bool do_read=pdb.readFromFilepointer(fp,false,0.1); if( !do_read ) plumed_merror("missing file " + reference );
  std::string mtype; parse("TYPE",mtype);

  if( pdb.getPositions().size()>0 ) {
      // And now create the rmsd object
      std::string rmsd_line =  getShortcutLabel() + "_at: RMSD DISPLACEMENT SQUARED NUMBER=1 REFERENCE=" + reference;
      bool nopbc; parseFlag("NOPBC",nopbc); if(nopbc) rmsd_line += " NOPBC";
      // Now create the RMSD object
      readInputLine( rmsd_line + " TYPE=" + mtype );
  }
  std::vector<std::string> argnames; parseVector("ARG",argnames); unsigned nargs=0; std::string instargs, refargs; std::vector<Value*> theargs; 
  if( argnames.size()>0 ) ActionWithArguments::interpretArgumentList( argnames, plumed.getActionSet(), this, theargs );
  for(unsigned i=0; i<theargs.size(); ++i) {
      std::vector<double> argval( theargs[i]->getNumberOfValues() ); nargs += theargs[i]->getNumberOfValues();
      if( !pdb.getArgumentValue(argnames[i], argval) ) error("argument " + argnames[i] + " was not set in pdb input");
      std::string arg_str; Tools::convert( argval[0], arg_str ); 
      for(unsigned j=1; j<argval.size(); ++j) { std::string nnn; Tools::convert( argval[j], nnn ); arg_str += "," + nnn; } 
      readInputLine( getShortcutLabel() + "_ref_" + argnames[i] + ": CONSTANT VALUES=" + arg_str );
      if( i==0 ) { instargs=" ARG1=" + argnames[i]; refargs=" ARG2=" + getShortcutLabel() + "_ref_" + argnames[i]; }
      else { instargs +="," + argnames[i]; refargs +="," + getShortcutLabel() + "_ref_" + argnames[i]; }
  }
  if( argnames.size()>0 ) readInputLine( getShortcutLabel() + "_argdist: EUCLIDEAN_DISTANCE SQUARED" + instargs + refargs );
  if( pdb.getPositions().size()>0 && argnames.size()>0 ) {
      readInputLine( getShortcutLabel() + ": CONCATENATE ARG=" + getShortcutLabel() + "_at.disp," + getShortcutLabel() + "_argdist_diffT");
      readInputLine( getShortcutLabel() + "_dist: COMBINE ARG=" + getShortcutLabel() + "_at.dist," + getShortcutLabel() + "_argdist PERIODIC=NO");
  } 

  // Get the displace stuff
  std::vector<double> displace( pdb.getBeta() ); double dtot = 0;
  for(unsigned i=0;i<displace.size();++i) dtot += displace[i];
  for(unsigned i=0;i<displace.size();++i) displace[i] = displace[i] / dtot;

  // Now read in the directions and create matheval objects to compute the pca components
  unsigned nfram=1;
  while( do_read ) {
    std::vector<double> argdir(nargs); PDB mypdb; do_read=mypdb.readFromFilepointer(fp,plumed.usingNaturalUnits(),0.1/plumed.getUnits().getLength());
    if( do_read ) {
        std::string num; Tools::convert( nfram, num ); nfram++;
        // Normalize the eigenvector in the input
        double norm=0;
        for(unsigned i=0;i<mypdb.getPositions().size();++i) {
            norm += mypdb.getPositions()[i][0]*mypdb.getPositions()[i][0];
            norm += mypdb.getPositions()[i][1]*mypdb.getPositions()[i][1];
            norm += mypdb.getPositions()[i][2]*mypdb.getPositions()[i][2];
        }
        unsigned k=0;
        for(unsigned i=0; i<theargs.size(); ++i) {
          std::vector<double> argval( theargs[i]->getNumberOfValues() );
          if( !mypdb.getArgumentValue(argnames[i], argval) ) error("argument " + argnames[i] + " was not set in pdb input");
          for(unsigned j=0; j<argval.size(); ++j) { argdir[k] = argval[j]; norm += argdir[k]*argdir[k]; k++; }
        }
        norm = sqrt( norm ); std::vector<double> normed_coeffs( 3*mypdb.getPositions().size() );
        for(unsigned i=0;i<mypdb.getPositions().size();++i) {
            if( mtype=="SIMPLE" ) {
                normed_coeffs[0*mypdb.getPositions().size()+i] = mypdb.getPositions()[i][0] / norm;
                normed_coeffs[1*mypdb.getPositions().size()+i] = mypdb.getPositions()[i][1] / norm;
                normed_coeffs[2*mypdb.getPositions().size()+i] = mypdb.getPositions()[i][2] / norm;
            } else {
                normed_coeffs[0*mypdb.getPositions().size()+i] = sqrt(displace[i])*mypdb.getPositions()[i][0] / norm;
                normed_coeffs[1*mypdb.getPositions().size()+i] = sqrt(displace[i])*mypdb.getPositions()[i][1] / norm;
                normed_coeffs[2*mypdb.getPositions().size()+i] = sqrt(displace[i])*mypdb.getPositions()[i][2] / norm;
            }
        }
        std::string coeff1, pvec; 
        if( mypdb.getPositions().size()>0 ) {
            Tools::convert( normed_coeffs[0], pvec );
            for(unsigned i=1;i<normed_coeffs.size();++i) {
                Tools::convert( normed_coeffs[i], coeff1 );
                pvec += "," + coeff1;
            }
            for(unsigned i=0; i<argdir.size(); ++i) { Tools::convert( argdir[i] / norm, coeff1 ); pvec += "," + coeff1; }
        } else if( argnames.size()>0 ) {
            Tools::convert( argdir[0] / norm, pvec ); for(unsigned i=1; i<argdir.size(); ++i) { Tools::convert( argdir[i] / norm, coeff1 ); pvec += "," + coeff1; }
        }

        // Read in eigenvector
        readInputLine( getShortcutLabel() + "_peig-" + num + ": CONSTANT VALUES=" + pvec );
        // And calculate dot product
        if( pdb.getPositions().size()>0 && argnames.size()>0 ) {
            readInputLine( getShortcutLabel() + "_prodeig-" + num + ": CUSTOM ARG=" + getShortcutLabel() + "_peig-" + num + "," + getShortcutLabel() + " FUNC=x*y PERIODIC=NO");
        } else if( pdb.getPositions().size()>0 ) {
            readInputLine( getShortcutLabel() + "_prodeig-" + num + ": CUSTOM ARG=" + getShortcutLabel() + "_peig-" + num + "," + getShortcutLabel() + "_at.disp FUNC=x*y PERIODIC=NO");
        } else if( argnames.size()>0 ) {
            readInputLine( getShortcutLabel() + "_prodeig-" + num + ": CUSTOM ARG=" + getShortcutLabel() + "_peig-" + num + "," + getShortcutLabel() + "_argdist_diffT FUNC=x*y PERIODIC=NO");
        }
        readInputLine( getShortcutLabel() + "_eig-" + num + ": SUM ARG=" + getShortcutLabel() + "_prodeig-" + num + " PERIODIC=NO");
    } else { break; }
  }
  std::fclose(fp); std::string resid_inp;
  if( pdb.getPositions().size()>0 && argnames.size()>0 ) resid_inp = getShortcutLabel() + "_residual_2: COMBINE PERIODIC=NO ARG=" + getShortcutLabel() + "_dist";
  else if( pdb.getPositions().size()>0 ) resid_inp = getShortcutLabel() + "_residual_2: COMBINE PERIODIC=NO ARG=" + getShortcutLabel() + "_at.dist";
  else if( argnames.size()>0 ) resid_inp = getShortcutLabel() + "_residual_2: COMBINE PERIODIC=NO ARG=" + getShortcutLabel() + "_argdist";
  for(unsigned i=0;i<nfram-1;++i) {
      std::string num; Tools::convert( i+1, num ); resid_inp += "," + getShortcutLabel() + "_eig-" + num;
  }
  resid_inp += " COEFFICIENTS=1"; for(unsigned i=0;i<nfram-1;++i) resid_inp +=",-1";
  resid_inp += " POWERS=1"; for(unsigned i=0;i<nfram-1;++i) resid_inp += ",2";
  readInputLine( resid_inp );
  readInputLine( getShortcutLabel() + "_residual: CUSTOM ARG=" + getShortcutLabel() + "_residual_2 FUNC=sqrt(x) PERIODIC=NO");
}

}
}


