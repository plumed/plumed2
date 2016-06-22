/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016 The plumed team
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
#include "Analysis.h"
#include "tools/Matrix.h"
#include "reference/Direction.h"
#include "reference/MetricRegister.h"
#include "reference/ReferenceConfiguration.h"
#include "reference/ReferenceValuePack.h"
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
to be analysed to the first frame in the trajectory.  This can be used to effectively remove translational and/or rotational motions from 
consideration.  The resulting principal components thus describe vibrational motions of the molecule. 

If you wish to calculate the projection of a trajectory on a set of principal components calculated from this PCA action then the output can be 
used as input for the \ref PCAVARS action.

\par Examples

The following input instructs PLUMED to perform a principal component analysis in which the covariance matrix is calculated from changes in the positions
of the first 22 atoms.  The TYPE=OPTIMAL instruction ensures that translational and rotational degrees of freedom are removed from consideration.
The first two principal components will be output to a file called pca-comp.pdb.  Trajectory frames will be collected on every step and the PCA calculation
will be performed at the end of the simulation.

\verbatim
PCA METRIC=OPTIMAL ATOMS=1-22 STRIDE=1 USE_ALL_DATA NLOW_DIM=2 OFILE=pca-comp.pdb
\endverbatim

The following input instructs PLUMED to perform a principal component analysis in which the covariance matrix is calculated from chnages in the six distances
seen in the previous lines.  Notice that here the TYPE=EUCLIDEAN keyword is used to indicate that no alighment has to be done when calculating the various
elements of the covariance matrix from the input vectors.  In this calculation the first two principal components will be output to a file called pca-comp.pdb.
Trajectory frames will be collected every five steps and the PCA calculation is performed every 1000 steps.  Consequently, if you run a 2000 step simulation the 
PCA analysis will be performed twice.  The REWEIGHT_BIAS keyword in this input tells PLUMED that rather that ascribing a weight of one to each of the frames
when calculating averages and covariances a reweighting should be performed based and each frames' weight in these calculations should be determined based on 
the current value of the instantaneous bias (see \ref reweighting).  

\verbatim
d1: DISTANCE ATOMS=1,2
d2: DISTANCE ATOMS=1,3
d3: DISTANCE ATOMS=1,4
d4: DISTNACE ATOMS=2,3
d5: DISTANCE ATOMS=2,4
d6: DISTANCE ATOMS=3,4

PCA ARG=d1,d2,d3,d4,d5,d6 METRIC=EUCLIDEAN STRIDE=5 RUN=1000 NLOW_DIM=2 REWEIGHT_BIAS OFILE=pca-comp.pdb
\endverbatim

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace analysis {

class PCA : public Analysis {
private:
  unsigned ndim;
/// The position of the reference configuration (the one we align to)
  ReferenceConfiguration* myref;
/// The eigenvectors for the atomic displacements
  Matrix<Vector> atom_eigv;
/// The eigenvectors for the displacements in argument space
  Matrix<double> arg_eigv;
  std::string ofilename;
public:
  static void registerKeywords( Keywords& keys );
  explicit PCA(const ActionOptions&ao);
  ~PCA();
  void performAnalysis();
  void performTask( const unsigned& , const unsigned& , MultiValue& ) const { plumed_error(); }
};

PLUMED_REGISTER_ACTION(PCA,"PCA")

void PCA::registerKeywords( Keywords& keys ){
  Analysis::registerKeywords( keys );
  keys.add("compulsory","NLOW_DIM","number of PCA coordinates required");
  keys.add("compulsory","OFILE","the file on which to output the eigenvectors");
}

PCA::PCA(const ActionOptions&ao):
PLUMED_ANALYSIS_INIT(ao)
{
  // Setup reference configuration
  log.printf("  performing PCA analysis using %s metric \n", getMetricName().c_str() );
  myref = metricRegister().create<ReferenceConfiguration>( getMetricName() );
  std::vector<std::string> argnames( getNumberOfArguments() );
  for(unsigned i=0;i<argnames.size();++i){
     if( getArguments()[i]->isPeriodic() ) error("cannot run PCA with periodic variables");
     argnames[i] = getArguments()[i]->getName();
  }
  myref->setNamesAndAtomNumbers( getAbsoluteIndexes(), argnames );

  parse("NLOW_DIM",ndim);
  if( getNumberOfAtoms()>0 ) atom_eigv.resize( ndim, getNumberOfAtoms() );
  if( getNumberOfArguments()>0 ) arg_eigv.resize( ndim, getNumberOfArguments() );

  // Read stuff for output file
  parseOutputFile("OFILE",ofilename); 
  checkRead();
}

PCA::~PCA(){
  delete myref;
}

void PCA::performAnalysis(){
  // Align everything to the first frame
  MultiValue myval( 1, getNumberOfArguments() + 3*getNumberOfAtoms() + 9 );
  ReferenceValuePack mypack( getNumberOfArguments(), getNumberOfAtoms(), myval );
  for(unsigned i=0;i<getNumberOfAtoms();++i) mypack.setAtomIndex( i, i );
  // Setup some PCA storage 
  data[0]->setupPCAStorage ( mypack );

  // Create some arrays to store the average position
  std::vector<double> sarg( getNumberOfArguments(), 0 );
  std::vector<Vector> spos( getNumberOfAtoms() );
  for(unsigned i=0;i<getNumberOfAtoms();++i) spos[i].zero();  
 
  // Calculate the average displacement from the first frame
  double norm=getWeight(0); 
  for(unsigned i=1;i<getNumberOfDataPoints();++i){
      double d = data[0]->calc( data[i]->getReferencePositions(), getPbc(), getArguments(), data[i]->getReferenceArguments(), mypack, true );
      // Accumulate average displacement of arguments (Here PBC could do fucked up things - really needs Berry Phase ) GAT
      for(unsigned j=0;j<getNumberOfArguments();++j) sarg[j] += 0.5*getWeight(i)*mypack.getArgumentDerivative(j); 
      // Accumulate average displacement of position
      for(unsigned j=0;j<getNumberOfAtoms();++j) spos[j] += getWeight(i)*mypack.getAtomsDisplacementVector()[j];
      norm += getWeight(i);
  }
  // Now normalise the displacements to get the average and add these to the first frame
  double inorm = 1.0 / norm ; 
  for(unsigned j=0;j<getNumberOfArguments();++j) sarg[j] = inorm*sarg[j] + data[0]->getReferenceArguments()[j];
  for(unsigned j=0;j<getNumberOfAtoms();++j) spos[j] = inorm*spos[j] + data[0]->getReferencePositions()[j]; 
  // And set the reference configuration
  std::vector<double> empty( getNumberOfArguments(), 1.0 ); myref->setReferenceConfig( spos, sarg, empty ); 

  // Now accumulate the covariance
  unsigned narg=getNumberOfArguments();
  Matrix<double> covar( getNumberOfArguments()+3*getNumberOfAtoms(), getNumberOfArguments()+3*getNumberOfAtoms() ); covar=0;
  for(unsigned i=0;i<getNumberOfDataPoints();++i){
      // double d = data[i]->calc( spos, getPbc(), getArguments(), sarg, mypack, true );
      double d = data[0]->calc( data[i]->getReferencePositions(), getPbc(), getArguments(), data[i]->getReferenceArguments(), mypack, true );
      for(unsigned jarg=0;jarg<getNumberOfArguments();++jarg){
         // Need sorting for PBC with GAT 
         double jarg_d = 0.5*mypack.getArgumentDerivative(jarg) + data[0]->getReferenceArguments()[jarg] - sarg[jarg];
         for(unsigned karg=0;karg<getNumberOfArguments();++karg){
            // Need sorting for PBC with GAT 
            double karg_d = 0.5*mypack.getArgumentDerivative(karg) + data[0]->getReferenceArguments()[karg] - sarg[karg];
            covar( jarg, karg ) += 0.25*getWeight(i)*jarg_d*karg_d; // mypack.getArgumentDerivative(jarg)*mypack.getArgumentDerivative(karg); 
         }
      }
      for(unsigned jat=0;jat<getNumberOfAtoms();++jat){ 
        for(unsigned jc=0;jc<3;++jc){
             double jdisplace = mypack.getAtomsDisplacementVector()[jat][jc] + data[0]->getReferencePositions()[jat][jc] - spos[jat][jc];
             for(unsigned kat=0;kat<getNumberOfAtoms();++kat){ 
                 for(unsigned kc=0;kc<3;++kc){
                    double kdisplace = mypack.getAtomsDisplacementVector()[kat][kc] + data[0]->getReferencePositions()[kat][kc] - spos[kat][kc];
                    covar( narg+3*jat + jc, narg+3*kat + kc ) += getWeight(i)*jdisplace*kdisplace; 
                 }
             }
         }
      }
  }
  // Normalise
  for(unsigned i=0;i<covar.nrows();++i){
      for(unsigned j=0;j<covar.ncols();++j) covar(i,j) *= inorm; 
  }

  // Diagonalise the covariance
  std::vector<double> eigval( getNumberOfArguments()+3*getNumberOfAtoms() );
  Matrix<double> eigvec( getNumberOfArguments()+3*getNumberOfAtoms(), getNumberOfArguments()+3*getNumberOfAtoms() );
  diagMat( covar, eigval, eigvec );

  // Open an output file
  OFile ofile; ofile.link(*this); ofile.setBackupString("analysis");
  ofile.open( ofilename ); 
  // Output the reference configuration
  myref->print( ofile, getOutputFormat(), atoms.getUnits().getLength()/0.1 );   

  // Store and print the eigenvectors
  std::vector<Vector> tmp_atoms( getNumberOfAtoms() );
  std::vector<double> tmp_args( getNumberOfArguments() );
  Direction* tref = metricRegister().create<Direction>( "DIRECTION" );
  tref->setNamesAndAtomNumbers( getAbsoluteIndexes(), argument_names );
  for(unsigned dim=0;dim<ndim;++dim){
     unsigned idim = covar.ncols() - 1 - dim;
     for(unsigned i=0;i<getNumberOfArguments();++i) tmp_args[i]=arg_eigv(dim,i)=eigvec(idim,i);
     for(unsigned i=0;i<getNumberOfAtoms();++i){
         for(unsigned k=0;k<3;++k) tmp_atoms[i][k]=atom_eigv(dim,i)[k]=eigvec(idim,narg+3*i+k);
     }  
     tref->setDirection( tmp_atoms, tmp_args );
     tref->print( ofile, getOutputFormat(), atoms.getUnits().getLength()/0.1 ); 
  } 
  // Close the output file   
  delete tref; ofile.close();
}

}
}
