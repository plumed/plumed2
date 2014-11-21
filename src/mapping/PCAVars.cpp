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
#include "core/ActionWithValue.h"
#include "core/ActionAtomistic.h"
#include "core/ActionWithArguments.h"
#include "reference/MultiReferenceBase.h"
#include "reference/MetricRegister.h"
#include "core/ActionRegister.h"
#include "core/PlumedMain.h"
#include "tools/Pbc.h"

namespace PLMD {
namespace mapping {

class PCAVars :
  public ActionWithValue,
  public ActionAtomistic,
  public ActionWithArguments
  {
private:
/// The position of the reference configuration (the one we align to)
  ReferenceConfiguration* myref; 
/// The eigenvectors for the atomic displacements
  Matrix<Vector> atom_eigv;
/// The eigenvectors for the displacements in argument space
  Matrix<double> arg_eigv;
/// Tempory vector used for storing derivatives
  std::vector<Vector> tmpder;
/// Stuff for applying forces
  std::vector<double> forces, forcesToApply;
public:
  static void registerKeywords( Keywords& keys );
  PCAVars(const ActionOptions&);
  ~PCAVars();
  unsigned getNumberOfDerivatives();
  void lockRequests();
  void unlockRequests();
  void calculateNumericalDerivatives( ActionWithValue* a );
  void calculate();
  void apply();
};

PLUMED_REGISTER_ACTION(PCAVars,"PCAVARS")

void PCAVars::registerKeywords( Keywords& keys ){
  Action::registerKeywords( keys );
  ActionWithValue::registerKeywords( keys );
  ActionAtomistic::registerKeywords( keys );
  ActionWithArguments::registerKeywords( keys );
  componentsAreNotOptional(keys);
  keys.addOutputComponent("eig-","default","the projections on each eigenvalue are stored on values labeled eig-1, eig-2, ...");
  keys.add("compulsory","REFERENCE","a pdb file containing the reference configuration and configurations that define the directions for each eigenvector");
  keys.add("compulsory","TYPE","OPTIMAL","The method we are using for alignment to the reference structure");
}

PCAVars::PCAVars(const ActionOptions& ao):
Action(ao),
ActionWithValue(ao),
ActionAtomistic(ao),
ActionWithArguments(ao)
{

  // What type of distance are we calculating
  std::string mtype; parse("TYPE",mtype);

  // Open reference file
  std::string reference; parse("REFERENCE",reference);
  FILE* fp=fopen(reference.c_str(),"r");
  if(!fp) error("could not open reference file " + reference );

  // Read all reference configurations 
  MultiReferenceBase myframes( "", false );
  bool do_read=true; unsigned nfram=0;
  while (do_read){
     PDB mypdb;
     // Read the pdb file
     do_read=mypdb.readFromFilepointer(fp,plumed.getAtoms().usingNaturalUnits(),0.1/atoms.getUnits().getLength());
     // Fix argument names
     expandArgKeywordInPDB( mypdb );
     if(do_read){
        if( nfram==0 ){
           myref = metricRegister().create<ReferenceConfiguration>( mtype, mypdb );
           if( myref->isDirection() ) error("first frame should be reference configuration - not direction of vector");
           if( !myref->pcaIsEnabledForThisReference() ) error("can't do PCA with reference type " + mtype );
           std::vector<std::string> remarks( mypdb.getRemark() ); std::string rtype;
           bool found=Tools::parse( remarks, "TYPE", rtype ); 
           if(!found){ std::vector<std::string> newrem(1); newrem[0]="TYPE="+mtype; mypdb.addRemark(newrem); }
           myframes.readFrame( mypdb );
        } else myframes.readFrame( mypdb ); 
        nfram++;
     } else {
        break;
     }
  }
  fclose(fp);

  if( nfram<2 ) error("no eigenvectors were specified");
  log.printf("  found %d eigenvectors in file %s \n",nfram-1,reference.c_str() );

  // Finish the setup of the mapping object
  // Get the arguments and atoms that are required
  std::vector<AtomNumber> atoms; std::vector<std::string> args;
  myframes.getAtomAndArgumentRequirements( atoms, args );
  requestAtoms( atoms ); std::vector<Value*> req_args;
  interpretArgumentList( args, req_args ); requestArguments( req_args );
  // Resize all derivative arrays
  myframes.setNumberOfAtomsAndArguments( atoms.size(), args.size() );

  // Retrieve the position of the first frame, as we use this for alignment
  myref->setNamesAndAtomNumbers( atoms, args );
  myref->setNumberOfAtoms( atoms.size() ); myref->setNumberOfArguments( args.size() );
  // Check there are no periodic arguments
  for(unsigned i=0;i<getNumberOfArguments();++i){
      if( getPntrToArgument(i)->isPeriodic() ) error("cannot use periodic variables in pca projections");
  }
  checkRead();

  // Resize the matrices that will hold our eivenvectors 
  if( getNumberOfAtoms()>0 ) atom_eigv.resize( nfram-1, getNumberOfAtoms() ); 
  if( getNumberOfArguments()>0 ) arg_eigv.resize( nfram-1, getNumberOfArguments() );

  // Create fake periodic boundary condition (these would only be used for DRMSD which is not allowed)
  Pbc fake_pbc; 
  // Now calculate the eigenvectors 
  for(unsigned i=1;i<nfram;++i){
      // Calculate distance from reference configuration
      double dist=myframes.getFrame(i)->calc( myref->getReferencePositions(), fake_pbc, getArguments(), myref->getReferenceArguments(), true );

      // Calculate the length of the vector for normalization
      double tmp, norm=0.0;
      for(unsigned j=0;j<getNumberOfAtoms();++j){ 
         for(unsigned k=0;k<3;++k){ tmp = myframes.getFrame(i)->getAtomicDisplacement(j)[k]; norm+=tmp*tmp; } 
      }
      for(unsigned j=0;j<getNumberOfArguments();++j){ tmp = 0.5*myframes.getFrame(i)->getArgumentDerivative(j); norm+=tmp*tmp; }

      // Normalize the eigevector
      norm = 1.0; //  / sqrt(norm);
      for(unsigned j=0;j<getNumberOfAtoms();++j) atom_eigv(i-1,j) = norm*myframes.getFrame(i)->getAtomicDisplacement(j); 
      for(unsigned j=0;j<getNumberOfArguments();++j) arg_eigv(i-1,j) = -0.5*norm*myframes.getFrame(i)->getArgumentDerivative(j); 

      // Create a component to store the output
      std::string num; Tools::convert( i, num );
      addComponentWithDerivatives("eig-"+num); componentIsNotPeriodic("eig-"+num);
  }

  // Get appropriate number of derivatives
  unsigned nder;
  if( getNumberOfAtoms()>0 ){
      nder = 3*getNumberOfAtoms() + 9 + getNumberOfArguments();
  } else {
      nder = getNumberOfArguments();
  }

  // Resize all derivative arrays
  forces.resize( nder ); forcesToApply.resize( nder ); tmpder.resize( getNumberOfAtoms() );
  for(unsigned i=0;i<getNumberOfComponents();++i) getPntrToComponent(i)->resizeDerivatives(nder);
}

PCAVars::~PCAVars(){
   delete myref;
}

unsigned PCAVars::getNumberOfDerivatives(){
  if( getNumberOfAtoms()>0 ){
      return 3*getNumberOfAtoms() + 9 + getNumberOfArguments();
  } 
  return getNumberOfArguments();
}    

void PCAVars::lockRequests(){
  ActionWithArguments::lockRequests();
  ActionAtomistic::lockRequests();
}

void PCAVars::unlockRequests(){          
  ActionWithArguments::unlockRequests(); 
  ActionAtomistic::unlockRequests();
} 

void PCAVars::calculate(){
  // Calculate distance between instaneous configuration and reference
  double dist = myref->calculate( getPositions(), getPbc(), getArguments(), true );

  // Now calculate projections on pca vectors
  unsigned nargs=getNumberOfArguments(); Vector adif, ader; Tensor fvir, tvir;
  for(unsigned i=0;i<getNumberOfComponents();++i){
      double proj=0; tvir.zero();
      for(unsigned j=0;j<getNumberOfArguments();++j){
          proj+=arg_eigv(i,j)*0.5*myref->getArgumentDerivative(j);
          getPntrToComponent(i)->addDerivative( j, arg_eigv(i,j) );
      }
      if( getNumberOfAtoms()>0 ){
         proj += myref->projectAtomicDisplacementOnVector( i, atom_eigv, getPositions(), tmpder );
         for(unsigned j=0;j<getNumberOfAtoms();++j){
            for(unsigned k=0;k<3;++k) getPntrToComponent(i)->addDerivative( nargs + 3*j+k, tmpder[j][k] );
            tvir += -1.0*Tensor( getPosition(j), tmpder[j] );
         }
         plumed_assert( !myref->getVirial( fvir ) );
         for(unsigned j=0;j<3;++j){
            for(unsigned k=0;k<3;++k) getPntrToComponent(i)->addDerivative( nargs + 3*getNumberOfAtoms() + 3*j + k, tvir(j,k) );
         }
      }
      getPntrToComponent(i)->set( proj );
  }
}

void PCAVars::calculateNumericalDerivatives( ActionWithValue* a ){
  if( getNumberOfArguments()>0 ){
     ActionWithArguments::calculateNumericalDerivatives( a );
  }
  if( getNumberOfAtoms()>0 ){
     Matrix<double> save_derivatives( getNumberOfComponents(), getNumberOfArguments() );
     for(unsigned j=0;j<getNumberOfComponents();++j){
        for(unsigned i=0;i<getNumberOfArguments();++i) save_derivatives(j,i)=getPntrToComponent(j)->getDerivative(i);
     }
     calculateAtomicNumericalDerivatives( a, getNumberOfArguments() );
     for(unsigned j=0;j<getNumberOfComponents();++j){
        for(unsigned i=0;i<getNumberOfArguments();++i) getPntrToComponent(j)->addDerivative( i, save_derivatives(j,i) );
     }
  }
} 

void PCAVars::apply(){

  bool wasforced=false;
  for(unsigned i=0;i<getNumberOfComponents();++i){
     if( getPntrToComponent(i)->applyForce( forces ) ){
         wasforced=true;
         for(unsigned i=0;i<forces.size();++i) forcesToApply[i]+=forces[i];
     }
  }
  if( wasforced ){
     addForcesOnArguments( forcesToApply );
     if( getNumberOfAtoms()>0 ) setForcesOnAtoms( forcesToApply, getNumberOfArguments() );
  }

}

}
}
