/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2019 The plumed team
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
#include "PCA.h"
#include "tools/Matrix.h"
#include "reference/MetricRegister.h"
#include "reference/ReferenceValuePack.h"
#include "analysis/ReadAnalysisFrames.h"
#include "core/ActionRegister.h"

//+PLUMEDOC DIMRED PCA
/*
Perform principal component analysis (PCA) using either the positions of the atoms a large number of collective variables as input.

Principal component analysis is a statistical technique that uses an orthogonal transformation to convert a set of observations of
poorly correlated variables into a set of linearly uncorrelated variables.  You can read more about the specifics of this technique
here: https://en.wikipedia.org/wiki/Principal_component_analysis

When used with molecular dynamics simulations a set of frames taken from the trajectory, \f$\{X_i\}\f$, or the values of
a number of collective variables which are calculated from the trajectory frames are used as input.  In this second instance your
input to the PCA analysis algorithm is thus a set of high-dimensional vectors of collective variables.  However, if
collective variables are calculated from the positions of the atoms or if the positions are used directly the assumption is that
this input trajectory is a set of poorly correlated (high-dimensional) vectors.  After principal component analysis has been
performed the output is a set of orthogonal vectors that describe the directions in which the largest motions have been seen.
In other words, principal component analysis provides a method for lowering the dimensionality of the data contained in a trajectory.
These output directions are some linear combination of the \f$x\f$, \f$y\f$ and \f$z\f$ positions if the positions were used as input
or some linear combination of the input collective variables if a high-dimensional vector of collective variables was used as input.

As explained on the Wikipedia page you must calculate the average and covariance for each of the input coordinates.  In other words, you must
calculate the average structure and the amount the system fluctuates around this average structure.  The problem in doing so when the
\f$x\f$, \f$y\f$ and \f$z\f$ coordinates of a molecule are used as input is that the majority of the changes in the positions of the
atoms comes from the translational and rotational degrees of freedom of the molecule.  The first six principal components will thus, most likely,
be uninteresting.  Consequently, to remedy this problem PLUMED provides the functionality to perform an RMSD alignment of the all the structures
to be analyzed to the first frame in the trajectory.  This can be used to effectively remove translational and/or rotational motions from
consideration.  The resulting principal components thus describe vibrational motions of the molecule.

If you wish to calculate the projection of a trajectory on a set of principal components calculated from this PCA action then the output can be
used as input for the \ref PCAVARS action.

\par Examples

The following input instructs PLUMED to perform a principal component analysis in which the covariance matrix is calculated from changes in the positions
of the first 22 atoms.  The TYPE=OPTIMAL instruction ensures that translational and rotational degrees of freedom are removed from consideration.
The first two principal components will be output to a file called PCA-comp.pdb.  Trajectory frames will be collected on every step and the PCA calculation
will be performed at the end of the simulation.

\plumedfile
PCA METRIC=OPTIMAL ATOMS=1-22 STRIDE=1 NLOW_DIM=2 OFILE=PCA-comp.pdb
\endplumedfile

The following input instructs PLUMED to perform a principal component analysis in which the covariance matrix is calculated from changes in the six distances
seen in the previous lines.  Notice that here the TYPE=EUCLIDEAN keyword is used to indicate that no alignment has to be done when calculating the various
elements of the covariance matrix from the input vectors.  In this calculation the first two principal components will be output to a file called PCA-comp.pdb.
Trajectory frames will be collected every five steps and the PCA calculation is performed every 1000 steps.  Consequently, if you run a 2000 step simulation the
PCA analysis will be performed twice.  The REWEIGHT_BIAS keyword in this input tells PLUMED that rather that ascribing a weight of one to each of the frames
when calculating averages and covariance matrices a reweighting should be performed based and each frames' weight in these calculations should be determined based on
the current value of the instantaneous bias (see \ref REWEIGHT_BIAS).

\plumedfile
d1: DISTANCE ATOMS=1,2
d2: DISTANCE ATOMS=1,3
d3: DISTANCE ATOMS=1,4
d4: DISTANCE ATOMS=2,3
d5: DISTANCE ATOMS=2,4
d6: DISTANCE ATOMS=3,4

PCA ARG=d1,d2,d3,d4,d5,d6 METRIC=EUCLIDEAN STRIDE=5 RUN=1000 NLOW_DIM=2 REWEIGHT_BIAS OFILE=PCA-comp.pdb
\endplumedfile

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace dimred {

PLUMED_REGISTER_ACTION(PCA,"PCA")

void PCA::registerKeywords( Keywords& keys ) {
  DimensionalityReductionBase::registerKeywords( keys ); keys.use("ARG"); keys.reset_style("ARG","optional");
  keys.add("compulsory","METRIC","EUCLIDEAN","the method that you are going to use to measure the distances between points");
  keys.add("atoms","ATOMS","the list of atoms that you are going to use in the measure of distance that you are using");
}

PCA::PCA(const ActionOptions&ao):
  Action(ao),
  DimensionalityReductionBase(ao)
{
  // Get the input PDB file from the underlying ReadAnalysisFrames object
  analysis::ReadAnalysisFrames* myframes = dynamic_cast<analysis::ReadAnalysisFrames*>( my_input_data );
  if( !myframes ) error("input to PCA should be ReadAnalysisFrames object");
  parse("METRIC",mtype); std::vector<AtomNumber> atoms;
  log.printf("  performing PCA analysis using %s metric \n", mtype.c_str() );
  if( my_input_data->getNumberOfAtoms()>0 ) {
    parseAtomList("ATOMS",atoms);
    if( atoms.size()!=0 ) {
      mypdb.setAtomNumbers( atoms );
      for(unsigned i=0; i<atoms.size(); ++i) {
        bool found=false;
        for(unsigned j=0; j<my_input_data->getAtomIndexes().size(); ++j) {
          if( my_input_data->getAtomIndexes()[j]==atoms[i] ) { found=true; break; }
        }
        if( !found ) {
          std::string num; Tools::convert( atoms[i].serial(), num );
          error("atom number " + num + " is not stored in any action that has been input");
        }
      }
      mypdb.addBlockEnd( atoms.size() );
    } else if( getNumberOfArguments()==0 ) {
      mypdb.setAtomNumbers( my_input_data->getAtomIndexes() );
      mypdb.addBlockEnd( my_input_data->getAtomIndexes().size() );
      if( mtype=="EUCLIDEAN" ) mtype="OPTIMAL";
    }
  }
  if( my_input_data->getArgumentNames().size()>0 ) {
    if( getNumberOfArguments()==0 && atoms.size()==0 ) {
      std::vector<std::string> argnames( my_input_data->getArgumentNames() );
      mypdb.setArgumentNames( argnames ); requestArguments( my_input_data->getArgumentList() );
    } else {
      std::vector<Value*> myargs( getArguments() );
      std::vector<std::string> inargnames( my_input_data->getArgumentNames() );
      std::vector<std::string> argnames( myargs.size() );
      for(unsigned i=0; i<myargs.size(); ++i) {
        argnames[i]=myargs[i]->getName();
        bool found=false;
        for(unsigned j=0; j<inargnames.size(); ++j) {
          if( argnames[i]==inargnames[j] ) { found=true; break; }
        }
        if( !found ) error("input named " + my_input_data->getLabel() + " does not store/calculate quantity named " + argnames[i] );
      }
      mypdb.setArgumentNames( argnames ); requestArguments( myargs );
    }
  }
  if( nlow==0 ) error("dimensionality of output not set");
  checkRead();
}

void PCA::performAnalysis() {
  // Align everything to the first frame
  my_input_data->getStoredData( 0, false ).transferDataToPDB( mypdb );
  auto myconf0=metricRegister().create<ReferenceConfiguration>(mtype, mypdb);
  MultiValue myval( 1, myconf0->getNumberOfReferenceArguments() + 3*myconf0->getNumberOfReferencePositions() + 9 );
  ReferenceValuePack mypack( myconf0->getNumberOfReferenceArguments(), myconf0->getNumberOfReferencePositions(), myval );
  for(unsigned i=0; i<myconf0->getNumberOfReferencePositions(); ++i) mypack.setAtomIndex( i, i );
  // Setup some PCA storage
  myconf0->setupPCAStorage ( mypack );
  std::vector<double> displace( myconf0->getNumberOfReferencePositions() );
  if( myconf0->getNumberOfReferencePositions()>0 ) {
    ReferenceAtoms* at = dynamic_cast<ReferenceAtoms*>( myconf0.get() );
    displace = at->getDisplace();
  }

  // Create some arrays to store the average position
  std::vector<double> sarg( myconf0->getNumberOfReferenceArguments(), 0 );
  std::vector<Vector> spos( myconf0->getNumberOfReferencePositions() );
  for(unsigned i=0; i<myconf0->getNumberOfReferencePositions(); ++i) spos[i].zero();

  // Calculate the average displacement from the first frame
  double norm=getWeight(0); std::vector<double> args( getNumberOfArguments() );
  for(unsigned i=1; i<getNumberOfDataPoints(); ++i) {
    my_input_data->getStoredData( i, false ).transferDataToPDB( mypdb );
    for(unsigned j=0; j<getArguments().size(); ++j) mypdb.getArgumentValue( getArguments()[j]->getName(), args[j] );
    double d = myconf0->calc( mypdb.getPositions(), getPbc(), getArguments(), args, mypack, true );
    // Accumulate average displacement of arguments (Here PBC could do fucked up things - really needs Berry Phase ) GAT
    for(unsigned j=0; j<myconf0->getNumberOfReferenceArguments(); ++j) sarg[j] += 0.5*getWeight(i)*mypack.getArgumentDerivative(j);
    // Accumulate average displacement of position
    for(unsigned j=0; j<myconf0->getNumberOfReferencePositions(); ++j) spos[j] += getWeight(i)*mypack.getAtomsDisplacementVector()[j] / displace[j];
    norm += getWeight(i);
  }
  // Now normalise the displacements to get the average and add these to the first frame
  double inorm = 1.0 / norm ;
  for(unsigned j=0; j<myconf0->getNumberOfReferenceArguments(); ++j) sarg[j] = inorm*sarg[j] + myconf0->getReferenceArguments()[j];
  for(unsigned j=0; j<myconf0->getNumberOfReferencePositions(); ++j) spos[j] = inorm*spos[j] + myconf0->getReferencePositions()[j];
  // Now accumulate the covariance
  unsigned narg=myconf0->getNumberOfReferenceArguments(), natoms=myconf0->getNumberOfReferencePositions();
  Matrix<double> covar( narg+3*natoms, narg+3*natoms ); covar=0;
  for(unsigned i=0; i<getNumberOfDataPoints(); ++i) {
    my_input_data->getStoredData( i, false ).transferDataToPDB( mypdb );
    for(unsigned j=0; j<getArguments().size(); ++j) mypdb.getArgumentValue( getArguments()[j]->getName(), args[j] );
    double d = myconf0->calc( mypdb.getPositions(), getPbc(), getArguments(), args, mypack, true );
    for(unsigned jarg=0; jarg<narg; ++jarg) {
      // Need sorting for PBC with GAT
      double jarg_d = 0.5*mypack.getArgumentDerivative(jarg) + myconf0->getReferenceArguments()[jarg] - sarg[jarg];
      for(unsigned karg=0; karg<narg; ++karg) {
        // Need sorting for PBC with GAT
        double karg_d = 0.5*mypack.getArgumentDerivative(karg) + myconf0->getReferenceArguments()[karg] - sarg[karg];
        covar( jarg, karg ) += 0.25*getWeight(i)*jarg_d*karg_d; // mypack.getArgumentDerivative(jarg)*mypack.getArgumentDerivative(karg);
      }
    }
    for(unsigned jat=0; jat<natoms; ++jat) {
      for(unsigned jc=0; jc<3; ++jc) {
        double jdisplace = mypack.getAtomsDisplacementVector()[jat][jc] / displace[jat] + myconf0->getReferencePositions()[jat][jc] - spos[jat][jc];
        for(unsigned kat=0; kat<natoms; ++kat) {
          for(unsigned kc=0; kc<3; ++kc) {
            double kdisplace = mypack.getAtomsDisplacementVector()[kat][kc] /displace[kat] + myconf0->getReferencePositions()[kat][kc] - spos[kat][kc];
            covar( narg+3*jat + jc, narg+3*kat + kc ) += getWeight(i)*jdisplace*kdisplace;
          }
        }
      }
    }
  }
  // Normalise
  for(unsigned i=0; i<covar.nrows(); ++i) {
    for(unsigned j=0; j<covar.ncols(); ++j) covar(i,j) *= inorm;
  }

  // Diagonalise the covariance
  std::vector<double> eigval( narg+3*natoms );
  Matrix<double> eigvec( narg+3*natoms, narg+3*natoms );
  diagMat( covar, eigval, eigvec );

  // Output the reference configuration
  mypdb.setAtomPositions( spos );
  for(unsigned j=0; j<sarg.size(); ++j) mypdb.setArgumentValue( getArguments()[j]->getName(), sarg[j] );
  // Reset the reference configuration
  myref = metricRegister().create<ReferenceConfiguration>( mtype, mypdb );

  // Store and print the eigenvectors
  std::vector<Vector> tmp_atoms( natoms );
  for(unsigned dim=0; dim<nlow; ++dim) {
    unsigned idim = covar.ncols() - 1 - dim;
    for(unsigned i=0; i<narg; ++i) mypdb.setArgumentValue( getArguments()[i]->getName(), eigvec(idim,i) );
    for(unsigned i=0; i<natoms; ++i) {
      for(unsigned k=0; k<3; ++k) tmp_atoms[i][k]=eigvec(idim,narg+3*i+k);
    }
    mypdb.setAtomPositions( tmp_atoms );
    // Create a direction object so that we can calculate other PCA components
    directions.push_back( Direction(ReferenceConfigurationOptions("DIRECTION")));
    directions[dim].read( mypdb );
  }
}

void PCA::getProjection( const unsigned& idata, std::vector<double>& point, double& weight ) {
  if( point.size()!=nlow ) point.resize( nlow );
  // Retrieve the weight
  weight = getWeight(idata);
  // Retrieve the data point
  getProjection( my_input_data->getStoredData( idata, false ), point );
}

void PCA::getProjection( analysis::DataCollectionObject& myidata, std::vector<double>& point ) {
  myidata.transferDataToPDB( mypdb ); std::vector<double> args( getArguments().size() );
  for(unsigned j=0; j<getArguments().size(); ++j) mypdb.getArgumentValue( getArguments()[j]->getName(), args[j] );
  // Create some storage space
  MultiValue myval( 1, 3*mypdb.getPositions().size() + args.size() + 9);
  ReferenceValuePack mypack( args.size(), mypdb.getPositions().size(), myval );
  for(unsigned i=0; i<mypdb.getPositions().size(); ++i) mypack.setAtomIndex( i, i );
  myref->setupPCAStorage( mypack );
  // And calculate
  myref->calculate( mypdb.getPositions(), getPbc(), getArguments(), mypack, true );
  for(unsigned i=0; i<nlow; ++i) point[i]=myref->projectDisplacementOnVector( directions[i], getArguments(), args, mypack );
}

}
}
