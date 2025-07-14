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
#include "Path.h"

namespace PLMD {
namespace mapping {

//+PLUMEDOC COLVAR PCAVARS
/*
Projection on principal component eigenvectors or other high dimensional linear subspace

As discussed in the documenation for [RMSD](RMSD.md) and [the refdist module](module_refdist.md) there are various different ways
of calculating the distance between the instantaneous structure adopted by the system and some high-dimensional, reference configuration.  The
problem with all these methods, as one gets further and further from the reference configuration, the
distance from it becomes a progressively poorer and poorer collective variable.  This happens because
the ``number" of structures at a distance $d$ from a reference configuration is proportional to $d^N$ in
an $N$ dimensional space.  Consequently, when $d$ is small the distance from the reference configuration
may well be a good collective variable.  However, when $d$ is large it is unlikely that the distance from the reference
structure is a good CV.  When the distance is large there will almost certainly be markedly different
configurations that have the same CV value and hence barriers in transverse degrees of
freedom.

For these reasons dimensionality reduction is often employed so a projection $\mathbf{s}$ of a high-dimensional configuration
$\mathbf{X}$ in a lower dimensionality space is found using a function:

$$
\mathbf{s} = F(\mathbf{X}-\mathbf{X}^{ref})
$$

where here we have introduced some high-dimensional reference configuration $\mathbf{X}^{ref}$. By far the simplest way to
do this is to use some linear operator for $F$.  That is to say we find a low-dimensional projection
by rotating the basis vectors using some linear algebra:

$$
\mathbf{s}_i = \sum_k A_{ik} ( X_{k} - X_{k}^{ref} )
$$

Here $\mathbf{X}-\mathbf{X}^{ref}$ is the displacement from the refernece configuration in the high dimenisonal space, which you can
calculate using the [RMSD](RMSD.md) action or by simply calculating difference between the instaneous values for a collection of variables
and the values of the variables for a particular reference configuration.  $A$ is then a $d$ by $D$ matrix where $D$ is the dimensionality
of the high dimensional space and $d$ is the dimensionality of the lower dimensional subspace.   This matrix, $A$,
can be found by various means including principal component analysis ([PCA](PCA.md)) and normal mode analysis.  In both these methods the
rows of $A$ would be the principle eigenvectors of a square matrix.  For PCA the covariance while for normal modes the Hessian.

## Examples

The following input calculates a projection on a linear subspace where the displacements
from the reference configuration are calculated using [RMSD](RMSD.md) and TYPE=OPTIMAL.  Consequently,
both translation of the center of mass of the atoms and rotation of the reference
frame are removed from the displacements that appear in the equatiosn above.  The matrix $A$ and the reference
configuration $R^{ref}$ are specified in the pdb input file reference.pdb and the
value of all projections (and the residual) are output to a file called colvar2.

```plumed
#SETTINGS INPUTFILES=regtest/mapping/rt-pca/reference.pdb
pca2: PCAVARS REFERENCE=regtest/mapping/rt-pca/reference.pdb TYPE=OPTIMAL
PRINT ARG=pca2.* FILE=colvar2
```

As you can see from the example above, the reference configurations are specified using a pdb file.  The first configuration that you provide is the reference configuration,
which is referred to in the above expressions as $X^{ref}$.  Subsequent configurations give the directions of row vectors that are contained in
the matrix $A$ above.  These directions are specified by giving a second configuration that describes the components of $A$ explicitly.

When you use an input like the one in the example above PLUMED assumes that the atoms for which RMSD displacements are being computed together form a molecule.
A procedure akin to that in [WHOLEMOLECULES](WHOLEMOLECULES.md) is thus performed to ensure that any bonds that are broken by the periodic boundary conditions are
reformed.  If you would like to turn this feature off for any reason you add the NOPBC flag to the input line as shown below:

```plumed
#SETTINGS INPUTFILES=regtest/mapping/rt-pca/reference.pdb
pca2: PCAVARS REFERENCE=regtest/mapping/rt-pca/reference.pdb NOPBC TYPE=OPTIMAL
PRINT ARG=pca2.* FILE=colvar2
```

Notice that the PCAVARS command is a shortcut.  If you look at how the shortcut in the above input is expanded you should be able to see how the command works
by calculating the RMSD displacement between the instantaneous and reference configuration and how those displacements are then projected on the eigenvector that
was specified in the second frame of the pdb input above. Understanding the expanded version of this shortcut command allows you to recognise that you can project
the displacement vector on any arbitrary vector.  For example in the input below two reference structures are provided in the pdb file. PLUMED calculates the RMSD
distance between these two reference configurations during setup and sets up a unit vector called eig that points along the director connecting the two RMSD structure.
During the calculation the vector connecting the instantaneous configuration and the first of the two reference configurations in computed.  This vector is then projected
on the unit vector connecting the two initial structures:

```plumed
#SETTINGS INPUTFILES=regtest/mapping/rt-pca-two-frames/two-frames.pdb

# Read in two reference configuratioms from PDB file
ref1: PDB2CONSTANT REFERENCE=regtest/mapping/rt-pca-two-frames/two-frames.pdb NUMBER=1
ref1T: TRANSPOSE ARG=ref1
ref2: PDB2CONSTANT REFERENCE=regtest/mapping/rt-pca-two-frames/two-frames.pdb NUMBER=2
# Calculate the displacement vector that takes you from ref1 to ref2
eigu: RMSD ARG=ref1T,ref2 DISPLACEMENT TYPE=OPTIMAL
# Normalise the reference vector
eigu2: CUSTOM ARG=eigu.disp FUNC=x*x PERIODIC=NO
eign2: SUM ARG=eigu2 PERIODIC=NO
eig: CUSTOM ARG=eigu.disp,eign2 FUNC=x/sqrt(y) PERIODIC=NO
eigT: TRANSPOSE ARG=eig
# Everything prior to this point is only done in setup.  From here onwards we have the commands that are done during the main calculate loop.

# Calculate the RMSD displacement between the instaneous structure and the first reference structure
rmsd: RMSD REFERENCE=regtest/mapping/rt-pca-two-frames/two-frames.pdb NUMBER=1 TYPE=OPTIMAL DISPLACEMENT SQUARED
# Project the displacement computed above on the director of the vector that connects reference structure ref1 to refeference structure ref2
pca: MATRIX_VECTOR_PRODUCT ARG=eigT,rmsd.disp

# Print the final CV to a file
PRINT ARG=pca FILE=colvar
```

You can also project vectors of differences of arguments on reference vectors.  For example, the input below can be used to look at the projection
of the vector connecting the instantanous configuraiton to a reference point in CV on a reference vector that has been specified in the PDB file.

```plumed
#SETTINGS INPUTFILES=regtest/mapping/rt-pca-args/epath.pdb

d1: DISTANCE ATOMS=1,2
d2: DISTANCE ATOMS=3,4
d3: DISTANCE ATOMS=5,6
pca: PCAVARS ARG=d1,d2,d3 REFERENCE=regtest/mapping/rt-pca-args/epath.pdb
PRINT ARG=pca_eig-1,pca_residual FILE=colvar
```

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
  keys.addInputKeyword("optional","ARG","scalar/vector","if there are arguments to be used specify them here");
  keys.addFlag("NOPBC",false,"do not use periodic boundary conditions when computing this quantity");
  keys.addOutputComponent("eig","default","vector","the projections on the eigenvalues");
  keys.addOutputComponent("residual","default","scalar","the residual distance that is not projected on any of the eigenvalues");
  keys.needsAction("RMSD");
  keys.needsAction("PDB2CONSTANT");
  keys.needsAction("TRANSPOSE");
  keys.needsAction("EUCLIDEAN_DISTANCE");
  keys.needsAction("CONCATENATE");
  keys.needsAction("COMBINE");
  keys.needsAction("CONSTANT");
  keys.needsAction("COMBINE");
  keys.needsAction("MATRIX_VECTOR_PRODUCT");
  keys.needsAction("CUSTOM");
  keys.needsAction("SUM");
  keys.needsAction("SELECT_COMPONENTS");
}

PCAVars::PCAVars( const ActionOptions& ao ):
  Action(ao),
  ActionShortcut(ao) {
  std::string reference;
  parse("REFERENCE",reference);
  // Create the object that holds the atomic positions by reading the first frame
  FILE* fp=std::fopen(reference.c_str(),"r");
  PDB pdb;
  if(!fp) {
    error("could not open reference file " + reference );
  }
  bool do_read=pdb.readFromFilepointer(fp,false,0.1);
  if( !do_read ) {
    plumed_merror("missing file " + reference );
  }
  std::string mtype;
  parse("TYPE",mtype);

  if( pdb.getPositions().size()>0 ) {
    // And now create the rmsd object
    std::string rmsd_line =  getShortcutLabel() + "_at: RMSD DISPLACEMENT SQUARED NUMBER=1 REFERENCE=" + reference;
    bool nopbc;
    parseFlag("NOPBC",nopbc);
    if(nopbc) {
      rmsd_line += " NOPBC";
    }
    // Now create the RMSD object
    readInputLine( rmsd_line + " TYPE=" + mtype );
  }
  std::vector<std::string> argnames;
  parseVector("ARG",argnames);
  unsigned nargs=0;
  std::string instargs, refargs;
  std::vector<Value*> theargs;
  if( argnames.size()>0 ) {
    ActionWithArguments::interpretArgumentList( argnames, plumed.getActionSet(), this, theargs );
  }
  for(unsigned i=0; i<theargs.size(); ++i) {
    std::string iargn = Path::fixArgumentName( theargs[i]->getName() );
    nargs += theargs[i]->getNumberOfValues();
    if( theargs[i]->getNumberOfValues()>1 ) {
      readInputLine( getShortcutLabel() + "_ref_" + iargn + "T: PDB2CONSTANT NUMBER=1 REFERENCE=" + reference + " ARG=" + theargs[i]->getName() );
      readInputLine( getShortcutLabel() + "_ref_" + iargn + ": TRANSPOSE ARG=" + getShortcutLabel() + "_ref_" + iargn + "T");
    } else {
      readInputLine( getShortcutLabel() + "_ref_" + iargn + ": PDB2CONSTANT NUMBER=1 REFERENCE=" + reference + " ARG=" + theargs[i]->getName() );
    }
    if( i==0 ) {
      instargs=" ARG1=" + theargs[i]->getName();
      refargs=" ARG2=" + getShortcutLabel() + "_ref_" + iargn;
    } else {
      instargs +="," + theargs[i]->getName();
      refargs +="," + getShortcutLabel() + "_ref_" + iargn;
    }
  }
  if( theargs.size()>0 ) {
    readInputLine( getShortcutLabel() + "_argdist: EUCLIDEAN_DISTANCE SQUARED" + instargs + refargs );
  }
  if( pdb.getPositions().size()>0 && theargs.size()>0 ) {
    readInputLine( getShortcutLabel() + ": CONCATENATE ARG=" + getShortcutLabel() + "_at.disp," + getShortcutLabel() + "_argdist_diffT");
    readInputLine( getShortcutLabel() + "_dist: COMBINE ARG=" + getShortcutLabel() + "_at.dist," + getShortcutLabel() + "_argdist PERIODIC=NO");
  }

  // Get the displace stuff
  std::vector<double> displace( pdb.getBeta() );
  double dtot = 0;
  for(unsigned i=0; i<displace.size(); ++i) {
    dtot += displace[i];
  }
  for(unsigned i=0; i<displace.size(); ++i) {
    displace[i] = displace[i] / dtot;
  }

  // Now read in the directions and create matheval objects to compute the pca components
  unsigned nfram=0, ncomp=0;
  std::string pvec;
  while( do_read ) {
    std::vector<double> argdir(nargs);
    PDB mypdb;
    do_read=mypdb.readFromFilepointer(fp,plumed.usingNaturalUnits(),0.1/plumed.getUnits().getLength());
    if( do_read ) {
      nfram++;
      // Normalize the eigenvector in the input
      double norm=0;
      for(unsigned i=0; i<mypdb.getPositions().size(); ++i) {
        norm += mypdb.getPositions()[i][0]*mypdb.getPositions()[i][0];
        norm += mypdb.getPositions()[i][1]*mypdb.getPositions()[i][1];
        norm += mypdb.getPositions()[i][2]*mypdb.getPositions()[i][2];
      }
      unsigned k=0;
      for(unsigned i=0; i<theargs.size(); ++i) {
        std::vector<double> argval( theargs[i]->getNumberOfValues() );
        if( !mypdb.getArgumentValue(theargs[i]->getName(), argval) ) {
          error("argument " + theargs[i]->getName() + " was not set in pdb input");
        }
        for(unsigned j=0; j<argval.size(); ++j) {
          argdir[k] = argval[j];
          norm += argdir[k]*argdir[k];
          k++;
        }
      }
      norm = sqrt( norm );
      std::vector<double> normed_coeffs( 3*mypdb.getPositions().size() );
      for(unsigned i=0; i<mypdb.getPositions().size(); ++i) {
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
      std::string coeff1;
      if( mypdb.getPositions().size()>0 ) {
        if( nfram==1 ) {
          Tools::convert( normed_coeffs[0], pvec );
        } else {
          Tools::convert( normed_coeffs[0], coeff1 );
          pvec += "," + coeff1;
        }
        for(unsigned i=1; i<normed_coeffs.size(); ++i) {
          Tools::convert( normed_coeffs[i], coeff1 );
          pvec += "," + coeff1;
        }
        for(unsigned i=0; i<argdir.size(); ++i) {
          Tools::convert( argdir[i] / norm, coeff1 );
          pvec += "," + coeff1;
        }
      } else if( theargs.size()>0 ) {
        if( nfram==1 ) {
          Tools::convert( argdir[0] / norm, pvec );
        } else {
          Tools::convert( argdir[0] / norm, coeff1 );
          pvec += "," + coeff1;
        }
        for(unsigned i=1; i<argdir.size(); ++i) {
          Tools::convert( argdir[i] / norm, coeff1 );
          pvec += "," + coeff1;
        }
      }
      ncomp = 3*mypdb.getPositions().size() + nargs;
    } else {
      break;
    }
  }
  std::fclose(fp);
  std::string neig, ncols;
  Tools::convert( nfram, neig );
  Tools::convert( ncomp, ncols );
  readInputLine( getShortcutLabel() + "_peig: CONSTANT VALUES=" + pvec + " NROWS=" + neig + " NCOLS=" + ncols );
  if( pdb.getPositions().size()>0 && theargs.size()>0 ) {
    readInputLine( getShortcutLabel() + "_eig: MATRIX_VECTOR_PRODUCT ARG=" + getShortcutLabel() + "_peig," + getShortcutLabel() );
  } else if( pdb.getPositions().size()>0 ) {
    readInputLine( getShortcutLabel() + "_eig: MATRIX_VECTOR_PRODUCT ARG=" + getShortcutLabel() + "_peig," + getShortcutLabel() + "_at.disp");
  } else if( theargs.size()>0 ) {
    readInputLine( getShortcutLabel() + "_eig: MATRIX_VECTOR_PRODUCT ARG=" + getShortcutLabel() + "_peig," + getShortcutLabel() + "_argdist_diffT");
  }
  for(unsigned i=0; i<nfram; ++i) {
    std::string num;
    Tools::convert( i+1, num );
    readInputLine( getShortcutLabel() + "_eig-" + num + ": SELECT_COMPONENTS ARG=" + getShortcutLabel() + "_eig COMPONENTS=" + num );
  }
  readInputLine( getShortcutLabel() + "_eig2: CUSTOM ARG=" + getShortcutLabel() + "_eig FUNC=x*x PERIODIC=NO");
  readInputLine( getShortcutLabel() + "_eigsum2: SUM ARG=" +  getShortcutLabel() + "_eig2 PERIODIC=NO");
  if( pdb.getPositions().size()>0 && theargs.size()>0 ) {
    readInputLine( getShortcutLabel() + "_residual: CUSTOM ARG=" + getShortcutLabel() + "_dist," + getShortcutLabel() + "_eigsum2 FUNC=sqrt(x-y) PERIODIC=NO");
  } else if( pdb.getPositions().size()>0 ) {
    readInputLine( getShortcutLabel() + "_residual: CUSTOM ARG=" + getShortcutLabel() + "_at.dist," + getShortcutLabel() + "_eigsum2 FUNC=sqrt(x-y) PERIODIC=NO");
  } else if( theargs.size()>0 ) {
    readInputLine( getShortcutLabel() + "_residual: CUSTOM ARG=" + getShortcutLabel() + "_argdist," + getShortcutLabel() + "_eigsum2 FUNC=sqrt(x-y) PERIODIC=NO");
  }
}

}
}


