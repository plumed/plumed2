/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2014-2020 The plumed team
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
#include "core/ActionWithValue.h"
#include "core/ActionAtomistic.h"
#include "core/ActionWithArguments.h"
#include "reference/MetricRegister.h"
#include "core/ActionRegister.h"
#include "core/PlumedMain.h"
#include "reference/Direction.h"
#include "tools/Pbc.h"

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
the matrix \f$A\f$ above.  These directions can be specified by specifying a second configuration - in this case a vector will
be constructed by calculating the displacement of this second configuration from the reference configuration.  A pdb input prepared
in this way would look as follows:

\auxfile{reference.pdb}
REMARK TYPE=OPTIMAL
ATOM      2  CH3 ACE     1      12.932 -14.718  -6.016  1.00  1.00
ATOM      5  C   ACE     1      21.312  -9.928  -5.946  1.00  1.00
ATOM      9  CA  ALA     2      19.462 -11.088  -8.986  1.00  1.00
ATOM     13  HB2 ALA     2      21.112 -10.688 -12.476  1.00  1.00
ATOM     15  C   ALA     2      19.422   7.978 -14.536  1.00  1.00
ATOM     20 HH31 NME     3      20.122  -9.928 -17.746  1.00  1.00
ATOM     21 HH32 NME     3      18.572 -13.148 -16.346  1.00  1.00
END
REMARK TYPE=OPTIMAL
ATOM      2  CH3 ACE     1      13.932 -14.718  -6.016  1.00  1.00
ATOM      5  C   ACE     1      20.312  -9.928  -5.946  1.00  1.00
ATOM      9  CA  ALA     2      18.462 -11.088  -8.986  1.00  1.00
ATOM     13  HB2 ALA     2      20.112 -11.688 -12.476  1.00  1.00
ATOM     15  C   ALA     2      19.422   7.978 -12.536  1.00  1.00
ATOM     20 HH31 NME     3      20.122  -9.928 -17.746  1.00  1.00
ATOM     21 HH32 NME     3      18.572 -13.148 -16.346  1.00  1.00
END
\endauxfile

Alternatively, the second configuration can specify the components of \f$A\f$ explicitly.  In this case you need to include the
keyword TYPE=DIRECTION in the remarks to the pdb as shown below.

\verbatim
ATOM      2  CH3 ACE     1      12.932 -14.718  -6.016  1.00  1.00
ATOM      5  C   ACE     1      21.312  -9.928  -5.946  1.00  1.00
ATOM      9  CA  ALA     2      19.462 -11.088  -8.986  1.00  1.00
ATOM     13  HB2 ALA     2      21.112 -10.688 -12.476  1.00  1.00
ATOM     15  C   ALA     2      19.422   7.978 -14.536  1.00  1.00
ATOM     20 HH31 NME     3      20.122  -9.928 -17.746  1.00  1.00
ATOM     21 HH32 NME     3      18.572 -13.148 -16.346  1.00  1.00
END
REMARK TYPE=DIRECTION
ATOM      2  CH3 ACE     1      0.1414  0.3334 -0.0302  1.00  0.00
ATOM      5  C   ACE     1      0.0893 -0.1095 -0.1434  1.00  0.00
ATOM      9  CA  ALA     2      0.0207 -0.321   0.0321  1.00  0.00
ATOM     13  HB2 ALA     2      0.0317 -0.6085  0.0783  1.00  0.00
ATOM     15  C   ALA     2      0.1282 -0.4792  0.0797  1.00  0.00
ATOM     20 HH31 NME     3      0.0053 -0.465   0.0309  1.00  0.00
ATOM     21 HH32 NME     3     -0.1019 -0.4261 -0.0082  1.00  0.00
END
\endverbatim

If your metric involves arguments the labels of these arguments in your plumed input file should be specified in the REMARKS
for each of the frames of your path.  An input file in this case might look like this:

\verbatim
DESCRIPTION: a pca eigenvector specified using the start point and direction in the HD space.
REMARK WEIGHT=1.0
REMARK ARG=d1,d2
REMARK d1=1.0 d2=1.0
END
REMARK TYPE=DIRECTION
REMARK ARG=d1,d2
REMARK d1=0.1 d2=0.25
END
\endverbatim

Here we are working with the EUCLIDEAN metric and notice that we have specified the components of \f$A\f$ using DIRECTION.
Consequently, the values of d1 and d2 in the second frame above do not specify a particular coordinate in the high-dimensional
space as in they do in the first frame.  Instead these values are the coefficients that can be used to construct a linear combination of d1 and d2.
If we wanted to specify the direction in this metric using the start and end point of the vector we would write:

\verbatim
DESCRIPTION: a pca eigenvector specified using the start and end point of a vector in the HD space.
REMARK WEIGHT=1.0
REMARK ARG=d1,d2
REMARK d1=1.0 d2=1.0
END
REMARK ARG=d1,d2
REMARK d1=1.1 d2=1.25
END
\endverbatim

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace mapping {

class PCAVars :
  public ActionWithValue,
  public ActionAtomistic,
  public ActionWithArguments
{
private:
/// The holders for the derivatives
  MultiValue myvals;
  ReferenceValuePack mypack;
/// The position of the reference configuration (the one we align to)
  std::unique_ptr<ReferenceConfiguration> myref;
/// The eigenvectors we are interested in
  std::vector<Direction> directions;
/// Stuff for applying forces
  std::vector<double> forces, forcesToApply;
  bool nopbc;
public:
  static void registerKeywords( Keywords& keys );
  explicit PCAVars(const ActionOptions&);
  unsigned getNumberOfDerivatives() override;
  void lockRequests() override;
  void unlockRequests() override;
  void calculateNumericalDerivatives( ActionWithValue* a ) override;
  void calculate() override;
  void apply() override;
};

PLUMED_REGISTER_ACTION(PCAVars,"PCAVARS")

void PCAVars::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys );
  ActionWithValue::registerKeywords( keys );
  ActionAtomistic::registerKeywords( keys );
  ActionWithArguments::registerKeywords( keys );
  componentsAreNotOptional(keys); keys.use("ARG");
  keys.addOutputComponent("eig","default","the projections on each eigenvalue are stored on values labeled eig-1, eig-2, ...");
  keys.addOutputComponent("residual","default","the distance of the configuration from the linear subspace defined "
                          "by the vectors, \\f$e_i\\f$, that are contained in the rows of \\f$A\\f$.  In other words this is "
                          "\\f$\\sqrt( r^2 - \\sum_i [\\mathbf{r}.\\mathbf{e_i}]^2)\\f$ where "
                          "\\f$r\\f$ is the distance between the instantaneous position and the "
                          "reference point.");
  keys.add("compulsory","REFERENCE","a pdb file containing the reference configuration and configurations that define the directions for each eigenvector");
  keys.add("compulsory","TYPE","OPTIMAL","The method we are using for alignment to the reference structure");
  keys.addFlag("NOPBC",false,"ignore the periodic boundary conditions when calculating distances");
}

PCAVars::PCAVars(const ActionOptions& ao):
  Action(ao),
  ActionWithValue(ao),
  ActionAtomistic(ao),
  ActionWithArguments(ao),
  myvals(1,0),
  mypack(0,0,myvals),
  nopbc(false)
{

  // What type of distance are we calculating
  std::string mtype; parse("TYPE",mtype);

  parseFlag("NOPBC",nopbc);

  // Open reference file
  std::string reference; parse("REFERENCE",reference);
  FILE* fp=fopen(reference.c_str(),"r");
  if(!fp) error("could not open reference file " + reference );

  // Read all reference configurations
  // MultiReferenceBase myframes( "", false );
  std::vector<std::unique_ptr<ReferenceConfiguration> > myframes;
  bool do_read=true; unsigned nfram=0;
  while (do_read) {
    PDB mypdb;
    // Read the pdb file
    do_read=mypdb.readFromFilepointer(fp,plumed.getAtoms().usingNaturalUnits(),0.1/atoms.getUnits().getLength());
    // Fix argument names
    if(do_read) {
      if( nfram==0 ) {
        myref=metricRegister().create<ReferenceConfiguration>( mtype, mypdb );
        Direction* tdir = dynamic_cast<Direction*>( myref.get() );
        if( tdir ) error("first frame should be reference configuration - not direction of vector");
        if( !myref->pcaIsEnabledForThisReference() ) error("can't do PCA with reference type " + mtype );
        // std::vector<std::string> remarks( mypdb.getRemark() ); std::string rtype;
        // bool found=Tools::parse( remarks, "TYPE", rtype );
        // if(!found){ std::vector<std::string> newrem(1); newrem[0]="TYPE="+mtype; mypdb.addRemark(newrem); }
        // myframes.push_back( metricRegister().create<ReferenceConfiguration>( "", mypdb ) );
      } else {
        auto mymsd = metricRegister().create<ReferenceConfiguration>( "", mypdb );
        myframes.emplace_back( std::move(mymsd) );
      }
      nfram++;
    } else {
      break;
    }
  }
  fclose(fp);

  if( nfram<=1 ) error("no eigenvectors were specified");
  log.printf("  found %u eigenvectors in file %s \n",nfram-1,reference.c_str() );

  // Finish the setup of the mapping object
  // Get the arguments and atoms that are required
  std::vector<AtomNumber> atoms; myref->getAtomRequests( atoms, false );
  std::vector<std::string> args; myref->getArgumentRequests( args, false );
  if( atoms.size()>0 ) {
    log.printf("  found %z atoms in input \n",atoms.size());
    log.printf("  with indices : ");
    for(unsigned i=0; i<atoms.size(); ++i) {
      if(i%25==0) log<<"\n";
      log.printf("%d ",atoms[i].serial());
    }
    log.printf("\n");
  }
  requestAtoms( atoms ); std::vector<Value*> req_args;
  interpretArgumentList( args, req_args ); requestArguments( req_args );

  // And now check that the atoms/arguments are the same in all the eigenvalues
  for(unsigned i=0; i<myframes.size(); ++i) { myframes[i]->getAtomRequests( atoms, false ); myframes[i]->getArgumentRequests( args, false ); }

  // Setup the derivative pack
  if( atoms.size()>0 ) myvals.resize( 1, args.size() + 3*atoms.size() + 9 );
  else myvals.resize( 1, args.size() );
  mypack.resize( args.size(), atoms.size() );
  for(unsigned i=0; i<atoms.size(); ++i) mypack.setAtomIndex( i, i );
  /// This sets up all the storage data required by PCA in the pack
  myref->setupPCAStorage( mypack );

  // Check there are no periodic arguments
  for(unsigned i=0; i<getNumberOfArguments(); ++i) {
    if( getPntrToArgument(i)->isPeriodic() ) error("cannot use periodic variables in pca projections");
  }
  // Work out if the user wants to normalise the input vector
  checkRead();

  if(nopbc) log.printf("  without periodic boundary conditions\n");
  else      log.printf("  using periodic boundary conditions\n");

  // Resize the matrices that will hold our eivenvectors
  PDB mypdb; mypdb.setAtomNumbers( atoms ); mypdb.addBlockEnd( atoms.size() );
  if( args.size()>0 ) mypdb.setArgumentNames( args );
  // Resize the matrices that will hold our eivenvectors
  for(unsigned i=0; i<myframes.size(); ++i) {
    directions.push_back( Direction(ReferenceConfigurationOptions("DIRECTION"))); directions[i].read( mypdb );
  }

  // Create fake periodic boundary condition (these would only be used for DRMSD which is not allowed)
  // Now calculate the eigenvectors
  for(unsigned i=0; i<myframes.size(); ++i) {
    // Calculate distance from reference configuration
    myframes[i]->extractDisplacementVector( myref->getReferencePositions(), getArguments(), myref->getReferenceArguments(), true, directions[i] );
    // Create a component to store the output
    std::string num; Tools::convert( i+1, num );
    addComponentWithDerivatives("eig-"+num); componentIsNotPeriodic("eig-"+num);
  }
  addComponentWithDerivatives("residual"); componentIsNotPeriodic("residual");

  // Get appropriate number of derivatives
  unsigned nder;
  if( getNumberOfAtoms()>0 ) {
    nder = 3*getNumberOfAtoms() + 9 + getNumberOfArguments();
  } else {
    nder = getNumberOfArguments();
  }

  // Resize all derivative arrays
  forces.resize( nder ); forcesToApply.resize( nder );
  for(int i=0; i<getNumberOfComponents(); ++i) getPntrToComponent(i)->resizeDerivatives(nder);
}

unsigned PCAVars::getNumberOfDerivatives() {
  if( getNumberOfAtoms()>0 ) {
    return 3*getNumberOfAtoms() + 9 + getNumberOfArguments();
  }
  return getNumberOfArguments();
}

void PCAVars::lockRequests() {
  ActionWithArguments::lockRequests();
  ActionAtomistic::lockRequests();
}

void PCAVars::unlockRequests() {
  ActionWithArguments::unlockRequests();
  ActionAtomistic::unlockRequests();
}

void PCAVars::calculate() {

  if(!nopbc && getNumberOfAtoms()>0) makeWhole();

  // Clear the reference value pack
  mypack.clear();
  // Calculate distance between instaneous configuration and reference
  double dist = myref->calculate( getPositions(), getPbc(), getArguments(), mypack, true );

  // Start accumulating residual by adding derivatives of distance
  Value* resid=getPntrToComponent( getNumberOfComponents()-1 ); unsigned nargs=getNumberOfArguments();
  for(unsigned j=0; j<getNumberOfArguments(); ++j) resid->addDerivative( j, mypack.getArgumentDerivative(j) );
  for(unsigned j=0; j<getNumberOfAtoms(); ++j) {
    Vector ader=mypack.getAtomDerivative( j );
    for(unsigned k=0; k<3; ++k) resid->addDerivative( nargs +3*j+k, ader[k] );
  }
  // Retrieve the values of all arguments
  std::vector<double> args( getNumberOfArguments() ); for(unsigned i=0; i<getNumberOfArguments(); ++i) args[i]=getArgument(i);

  // Now calculate projections on pca vectors
  Vector adif, ader; Tensor fvir, tvir;
  for(int i=0; i<getNumberOfComponents()-1; ++i) { // One less component as we also have residual
    double proj=myref->projectDisplacementOnVector( directions[i], getArguments(), args, mypack );
    // And now accumulate derivatives
    Value* eid=getPntrToComponent(i);
    for(unsigned j=0; j<getNumberOfArguments(); ++j) eid->addDerivative( j, mypack.getArgumentDerivative(j) );
    if( getNumberOfAtoms()>0 ) {
      tvir.zero();
      for(unsigned j=0; j<getNumberOfAtoms(); ++j) {
        Vector myader=mypack.getAtomDerivative(j);
        for(unsigned k=0; k<3; ++k) {
          eid->addDerivative( nargs + 3*j+k, myader[k] );
          resid->addDerivative( nargs + 3*j+k, -2*proj*myader[k] );
        }
        tvir += -1.0*Tensor( getPosition(j), myader );
      }
      for(unsigned j=0; j<3; ++j) {
        for(unsigned k=0; k<3; ++k) eid->addDerivative( nargs + 3*getNumberOfAtoms() + 3*j + k, tvir(j,k) );
      }
    }
    dist -= proj*proj; // Subtract square from total squared distance to get residual squared
    // Derivatives of residual
    for(unsigned j=0; j<getNumberOfArguments(); ++j) resid->addDerivative( j, -2*proj*eid->getDerivative(j) );
    // for(unsigned j=0;j<getNumberOfArguments();++j) resid->addDerivative( j, -2*proj*arg_eigv(i,j) );
    // And set final value
    getPntrToComponent(i)->set( proj );
  }
  dist=sqrt(dist);
  resid->set( dist );

  // Take square root of residual derivatives
  double prefactor = 0.5 / dist;
  for(unsigned j=0; j<getNumberOfArguments(); ++j) resid->setDerivative( j, prefactor*resid->getDerivative(j) );
  for(unsigned j=0; j<getNumberOfAtoms(); ++j) {
    for(unsigned k=0; k<3; ++k) resid->setDerivative( nargs + 3*j+k, prefactor*resid->getDerivative( nargs+3*j+k ) );
  }

  // And finally virial for residual
  if( getNumberOfAtoms()>0 ) {
    tvir.zero();
    for(unsigned j=0; j<getNumberOfAtoms(); ++j) {
      Vector ader; for(unsigned k=0; k<3; ++k) ader[k]=resid->getDerivative( nargs + 3*j+k );
      tvir += -1.0*Tensor( getPosition(j), ader );
    }
    for(unsigned j=0; j<3; ++j) {
      for(unsigned k=0; k<3; ++k) resid->addDerivative( nargs + 3*getNumberOfAtoms() + 3*j + k, tvir(j,k) );
    }
  }

}

void PCAVars::calculateNumericalDerivatives( ActionWithValue* a ) {
  if( getNumberOfArguments()>0 ) {
    ActionWithArguments::calculateNumericalDerivatives( a );
  }
  if( getNumberOfAtoms()>0 ) {
    Matrix<double> save_derivatives( getNumberOfComponents(), getNumberOfArguments() );
    for(int j=0; j<getNumberOfComponents(); ++j) {
      for(unsigned i=0; i<getNumberOfArguments(); ++i) save_derivatives(j,i)=getPntrToComponent(j)->getDerivative(i);
    }
    calculateAtomicNumericalDerivatives( a, getNumberOfArguments() );
    for(int j=0; j<getNumberOfComponents(); ++j) {
      for(unsigned i=0; i<getNumberOfArguments(); ++i) getPntrToComponent(j)->addDerivative( i, save_derivatives(j,i) );
    }
  }
}

void PCAVars::apply() {

  bool wasforced=false; forcesToApply.assign(forcesToApply.size(),0.0);
  for(int i=0; i<getNumberOfComponents(); ++i) {
    if( getPntrToComponent(i)->applyForce( forces ) ) {
      wasforced=true;
      for(unsigned i=0; i<forces.size(); ++i) forcesToApply[i]+=forces[i];
    }
  }
  if( wasforced ) {
    addForcesOnArguments( forcesToApply );
    if( getNumberOfAtoms()>0 ) setForcesOnAtoms( forcesToApply, getNumberOfArguments() );
  }

}

}
}
