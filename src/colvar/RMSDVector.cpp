/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2023 The plumed team
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
#include "RMSDVector.h"
#include "core/PlumedMain.h"
#include "core/ActionRegister.h"
#include "tools/PDB.h"

//+PLUMEDOC DCOLVAR RMSD_VECTOR
/*
Calculate the RMSD distance between the instaneous configuration and multiple reference configurations


\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace colvar {

PLUMED_REGISTER_ACTION(RMSDVector,"RMSD_VECTOR")

void RMSDVector::registerKeywords(Keywords& keys) {
  ActionWithVector::registerKeywords(keys);
  keys.setDisplayName("RMSD");
  keys.addInputKeyword("compulsory","ARG","vector/matrix","the labels of two actions that you are calculating the RMSD between");
  keys.add("compulsory","TYPE","SIMPLE","the manner in which RMSD alignment is performed.  Should be OPTIMAL or SIMPLE.");
  keys.add("compulsory","ALIGN","1.0","the weights to use when aligning to the reference structure");
  keys.add("compulsory","DISPLACE","1.0","the weights to use when calculating the displacement from the reference structure");
  keys.addFlag("SQUARED",false," This should be set if you want mean squared displacement instead of RMSD ");
  keys.addFlag("UNORMALIZED",false,"by default the mean sequare deviation or root mean square deviation is calculated.  If this option is given no averaging is done");
  keys.addFlag("DISPLACEMENT",false,"Calculate the vector of displacements instead of the length of this vector");
  keys.addOutputComponent("disp","DISPLACEMENT","vector/matrix","the vector of displacements for the atoms");
  keys.addOutputComponent("dist","DISPLACEMENT","scalar/vector","the RMSD distance the atoms have moved");
  PTM::registerKeywords( keys );
  keys.setValueDescription("scalar/vector","a vector containing the RMSD between the instantaneous structure and each of the reference structures that were input");
}

RMSDVector::RMSDVector(const ActionOptions&ao):
  Action(ao),
  ActionWithVector(ao),
  firststep(true),
  input(getPbc()),
  taskmanager(this) {
  if( getPntrToArgument(0)->getRank()!=1 ) {
    error("first argument should be vector");
  }
  if( getPntrToArgument(1)->getRank()<1 ) {
    error("second argument should be matrix or a vector");
  }
  if( getPntrToArgument(1)->getRank()==1 ) {
    if( getPntrToArgument(0)->getNumberOfValues()!=getPntrToArgument(1)->getNumberOfValues() ) {
      error("mismatch between sizes of input vectors");
    }
  } else if( getPntrToArgument(1)->getRank()==2 ) {
    if( getPntrToArgument(0)->getNumberOfValues()!=getPntrToArgument(1)->getShape()[1] ) {
      error("mismatch between sizes of input vectors");
    }
  }
  if( getPntrToArgument(0)->getNumberOfValues()%3!=0 ) {
    error("number of components in input arguments should be multiple of three");
  }

  RMSDVectorData myinput;
  unsigned natoms = getPntrToArgument(0)->getNumberOfValues() / 3;
  myinput.type.assign("SIMPLE");
  parse("TYPE",myinput.type);
  parseFlag("SQUARED",myinput.squared);
  myinput.align.resize( natoms );
  parseVector("ALIGN",myinput.align);
  myinput.displace.resize( natoms );
  parseVector("DISPLACE",myinput.displace);
  bool unorm=false;
  parseFlag("UNORMALIZED",unorm);
  norm_weights=!unorm;
  double wa=0, wd=0;
  myinput.sqrtdisplace.resize( myinput.displace.size() );
  for(unsigned i=0; i<myinput.align.size(); ++i) {
    wa+=myinput.align[i];
    wd+=myinput.displace[i];
  }

  if( wa>epsilon ) {
    double iwa = 1. / wa;
    for(unsigned i=0; i<myinput.align.size(); ++i) {
      myinput.align[i] *= iwa;
    }
  } else {
    double iwa = 1. / natoms;
    for(unsigned i=0; i<myinput.align.size(); ++i) {
      myinput.align[i] = iwa;
    }
  }
  if( wd>epsilon ) {
    if( !norm_weights ) {
      wd = 1;
    }
    double iwd = 1. / wd;
    for(unsigned i=0; i<myinput.align.size(); ++i) {
      myinput.displace[i] *= iwd;
    }
  } else {
    double iwd = 1. / natoms;
    for(unsigned i=0; i<myinput.align.size(); ++i) {
      myinput.displace[i] = iwd;
    }
  }
  for(unsigned i=0; i<myinput.align.size(); ++i) {
    myinput.sqrtdisplace[i] = sqrt(myinput.displace[i]);
  }

  parseFlag("DISPLACEMENT",myinput.displacement);
  if( myinput.displacement && (getPntrToArgument(1)->getRank()==1 || getPntrToArgument(1)->getShape()[0]<=1) ) {
    addComponentWithDerivatives("dist");
    componentIsNotPeriodic("dist");
    std::vector<std::size_t> shape( 1, getPntrToArgument(0)->getNumberOfValues() );
    addComponent( "disp", shape );
    componentIsNotPeriodic("disp");
  } else if( myinput.displacement ) {
    std::vector<std::size_t> shape( 1, getPntrToArgument(1)->getShape()[0] );
    addComponent( "dist", shape );
    componentIsNotPeriodic("dist");
    shape.resize(2);
    shape[0] = getPntrToArgument(1)->getShape()[0];
    shape[1] = getPntrToArgument(0)->getNumberOfValues();
    addComponent( "disp", shape );
    getPntrToComponent(1)->reshapeMatrixStore( shape[1] );
    componentIsNotPeriodic("disp");
  } else if( (getPntrToArgument(1)->getRank()==1 || getPntrToArgument(1)->getShape()[0]==1) ) {
    addValue();
    setNotPeriodic();
  } else {
    std::vector<std::size_t> shape( 1, getPntrToArgument(1)->getShape()[0] );
    addValue( shape );
    setNotPeriodic();
  }
  if( getPntrToArgument(1)->getRank()==1 || getPntrToArgument(1)->getNumberOfValues()==0 ) {
    myinput.myrmsd.resize(1);
  } else {
    myinput.myrmsd.resize( getPntrToArgument(1)->getShape()[0] );
  }

  if( getPntrToArgument(1)->getRank()==1 )
    log.printf("  calculating RMSD distance between %d atoms. Distance between the avectors of atoms in %s and %s\n",
               natoms, getPntrToArgument(0)->getName().c_str(), getPntrToArgument(1)->getName().c_str() );
  else
    log.printf("  calculating RMSD distance between %d sets of %d atoms. Distance between vector %s of atoms and matrix of configurations in %s\n",
               getPntrToArgument(1)->getShape()[0], natoms, getPntrToArgument(0)->getName().c_str(), getPntrToArgument(1)->getName().c_str() );
  log.printf("  method for alignment : %s \n",myinput.type.c_str() );
  if(myinput.squared) {
    log.printf("  chosen to use SQUARED option for MSD instead of RMSD\n");
  } else {
    log.printf("  using periodic boundary conditions\n");
  }
  // Setup the task manager
  int stride = 1;
  if( myinput.displacement ) {
    stride = 1 + 3*natoms;
  }
  taskmanager.setupParallelTaskManager( 0, 3*natoms, -1, stride );
  taskmanager.setActionInput( myinput );
}

unsigned RMSDVector::getNumberOfDerivatives() {
  return getPntrToArgument(0)->getNumberOfValues() + getPntrToArgument(1)->getNumberOfValues();
}

int RMSDVector::checkTaskIsActive( const unsigned& itask ) const {
  return 1;
}

void RMSDVector::setReferenceConfigurations() {
  unsigned natoms = getPntrToArgument(0)->getShape()[0] / 3;
  Vector center;
  RMSDVectorData& myinput=taskmanager.getActionInput();
  std::vector<Vector> pos( natoms );
  for(unsigned jconf=0; jconf<myinput.myrmsd.size(); ++jconf) {
    center.zero();
    for(unsigned i=0; i<pos.size(); ++i) {
      for(unsigned j=0; j<3; ++j) {
        pos[i][j] = getPntrToArgument(1)->get( (3*jconf+j)*pos.size() + i );
      }
      center+=pos[i]*myinput.align[i];
    }
    for(unsigned i=0; i<pos.size(); ++i) {
      pos[i] -= center;
    }
    myinput.myrmsd[jconf].clear();
    myinput.myrmsd[jconf].set(myinput.align,myinput.displace,pos,myinput.type,true,norm_weights);
  }
}

// calculator
void RMSDVector::calculate() {
  if( firststep || !getPntrToArgument(1)->isConstant() ) {
    setReferenceConfigurations();
    firststep=false;
  }

  if( getPntrToComponent(0)->getRank()==0 ) {
    input.noderiv = false;
    input.ncomponents = getNumberOfComponents();
    input.nscalars= 0;
    for(unsigned i=0; i<input.ncomponents; ++i) {
      input.nscalars += getPntrToComponent(i)->getNumberOfStoredValues();
    }
    getInputData( input_buffer );
    input.dataSize = input_buffer.size();
    input.inputdata = input_buffer.data();
    ArgumentsBookkeeping abk;
    abk.setupArguments( this );
    input.setupArguments( abk );

    std::vector<double> buffer;
    std::vector<double> values( input.nscalars );
    std::vector<double> deriv( input_buffer.size(), 0 );
    ParallelActionsOutput output( input.nscalars, values.data(), input_buffer.size(), deriv.data(), 0, buffer.data() );
    performTask( 0, taskmanager.getActionInput(), input, output );

    getPntrToComponent(0)->set( values[0] );
    if( getNumberOfComponents()==2 ) {
      Value* mydisp = getPntrToComponent(1);
      for(unsigned i=0; i<mydisp->getNumberOfStoredValues(); i++) {
        mydisp->set( i, values[i+1] );
      }
    }
    if( doNotCalculateDerivatives() ) {
      return;
    }

    Value* myval = getPntrToComponent(0);
    for(unsigned i=0; i<deriv.size(); i++) {
      myval->setDerivative( i, deriv[i] );
    }
  } else {
    taskmanager.runAllTasks();
  }
}

void RMSDVector::apply() {
  if( doNotCalculateDerivatives() ) {
    return;
  }

  if( getPntrToComponent(0)->getRank()==0 ) {
    Value* myval = getPntrToComponent(0);
    std::vector<double> deriv( input_buffer.size(), 0 );
    for(unsigned i=0; i<deriv.size(); ++i) {
      deriv[i] = myval->getDerivative( i );
    }

    std::vector<double> f( input.nscalars );
    const Value* disp = getConstPntrToComponent(1);
    f[0] = getConstPntrToComponent(0)->getForce(0);
    for(unsigned i=0; i<disp->getNumberOfStoredValues(); ++i) {
      f[i+1] = disp->getForce( i );
    }

    std::vector<double> ignore;
    std::vector<double> forces( getNumberOfDerivatives(), 0 );
    ForceOutput fout( forces, ignore );
    gatherForces( 0, taskmanager.getActionInput(), input, ForceInput( getNumberOfComponents(), input.nscalars, f.data(), deriv.size(), deriv.data() ), fout );

    unsigned ss=0;
    addForcesOnArguments( 0, forces, ss  );
  } else {
    ActionWithVector::apply();
  }
}

void RMSDVector::getPositionsFromInputData( const ParallelActionsInput& input, std::vector<Vector>& pos ) {
  View2D<const double,helpers::dynamic_extent,helpers::dynamic_extent> argpos( input.inputdata, 3, pos.size() );
  // At some stage it would be good to get rid of this and to be able to pass the 2D view directly to the methods in RMSD
  for(unsigned i=0; i<pos.size(); ++i) {
    pos[i][0] = argpos[0][i];
    pos[i][1] = argpos[1][i];
    pos[i][2] = argpos[2][i];
  }
}

void RMSDVector::performTask( std::size_t task_index,
                              const RMSDVectorData& actiondata,
                              ParallelActionsInput& input,
                              ParallelActionsOutput& output ) {
  std::size_t natoms = actiondata.align.size();
  std::vector<Vector> der(natoms), pos(natoms);
  getPositionsFromInputData( input, pos );
  if( actiondata.displacement && actiondata.type=="SIMPLE" ) {
    std::vector<Vector> direction( natoms );
    output.values[0] = actiondata.myrmsd[task_index].simpleAlignment( actiondata.align, actiondata.displace, View(pos.data(),pos.size()), actiondata.myrmsd[task_index].getReference(), der, direction, actiondata.squared );
    for(unsigned i=0; i<direction.size(); ++i) {
      output.values[1+i] = direction[i][0];
      output.values[1+natoms+i] = direction[i][1];
      output.values[1+2*natoms+i] = direction[i][2];
    }
  } else if( actiondata.displacement ) {
    if( !input.noderiv && actiondata.myrmsd.size()>1 ) {
      // This is activated on the backwards loop to ensure that we don't calculate the RMSD twice during that backwards loop
      return;
    }
    Tensor rot;
    std::vector<Vector> direction( natoms );
    Matrix<std::vector<Vector> > DRotDPos(3,3);
    std::vector<Vector> centeredpos( natoms ), centeredreference( natoms );
    output.values[0] = actiondata.myrmsd[task_index].calc_PCAelements( pos, der, rot, DRotDPos, direction, centeredpos, centeredreference, actiondata.squared );
    std::vector<Vector> ref( actiondata.myrmsd[task_index].getReference() );
    for(unsigned i=0; i<direction.size(); ++i) {
      output.values[1+i] = actiondata.sqrtdisplace[i]*( direction[i][0] - ref[i][0] );
      output.values[1+natoms+i] = actiondata.sqrtdisplace[i]*( direction[i][1] - ref[i][1] );
      output.values[1+2*natoms+i] = actiondata.sqrtdisplace[i]*( direction[i][2] - ref[i][2] );
    }
  } else {
    output.values[0] = actiondata.myrmsd[task_index].calculate( pos, der, actiondata.squared );
  }
  if( !input.noderiv ) {
    // This also could be removed if someone was willing to do some additional work
    for(unsigned i=0; i<natoms; ++i) {
      output.derivatives[i] = der[i][0];
      output.derivatives[natoms+i] = der[i][1];
      output.derivatives[2*natoms+i] = der[i][2];
    }
  }
}

void RMSDVector::transferStashToValues( const std::vector<double>& stash ) {
  if( getNumberOfComponents()==1 ) {
    ActionWithVector::transferStashToValues( stash );
    return;
  }
  std::size_t k=0;
  Value* dist = getPntrToComponent(0);
  Value* disp = getPntrToComponent(1);
  std::size_t ss = disp->getShape()[1];
  for(unsigned i=0; i<dist->getNumberOfStoredValues(); ++i) {
    dist->set(i,stash[k]);
    k++;
    for(unsigned j=0; j<ss; ++j) {
      disp->set( ss*i + j, stash[k] );
      k++;
    }
  }
}

void RMSDVector::transferForcesToStash( std::vector<double>& stash ) const {
  if( getNumberOfComponents()==1 ) {
    ActionWithVector::transferForcesToStash( stash );
    return;
  }
  std::size_t k=0;
  const Value* dist = getConstPntrToComponent(0);
  const Value* disp = getConstPntrToComponent(1);
  std::size_t ss = disp->getShape()[1];
  for(unsigned i=0; i<dist->getNumberOfStoredValues(); ++i) {
    stash[k] = dist->getForce(i);
    k++;
    for(unsigned j=0; j<ss; ++j) {
      stash[k] = disp->getForce( ss*i + j );
      k++;
    }
  }
}

void RMSDVector::gatherForces( std::size_t task_index,
                               const RMSDVectorData& actiondata,
                               const ParallelActionsInput& input,
                               const ForceInput& fdata,
                               ForceOutput forces ) {
  std::size_t natoms = actiondata.align.size();
  if( actiondata.displacement && actiondata.type=="SIMPLE" ) {
    Vector comforce;
    comforce.zero();
    for(unsigned i=0; i<natoms; ++i) {
      comforce[0] += actiondata.align[i]*fdata.force[1+i];
      comforce[1] += actiondata.align[i]*fdata.force[1+natoms+i];
      comforce[2] += actiondata.align[i]*fdata.force[1+2*natoms+i];
    }
    for(unsigned i=0; i<natoms; ++i) {
      forces.thread_safe[i] += fdata.force[1+i] - comforce[0];
      forces.thread_safe[natoms+i] += fdata.force[1+natoms+i] - comforce[1];
      forces.thread_safe[2*natoms+i] += fdata.force[1+2*natoms+i] - comforce[2];
    }
  } else if( actiondata.displacement ) {
    Tensor rot;
    std::vector<Vector> der(natoms), pos(natoms);
    getPositionsFromInputData( input, pos );
    std::vector<Vector> direction( natoms );
    Matrix<std::vector<Vector> > DRotDPos(3,3);
    std::vector<Vector> centeredpos( natoms ), centeredreference( natoms );
    double rmsd = actiondata.myrmsd[task_index].calc_PCAelements( pos, der, rot, DRotDPos, direction, centeredpos, centeredreference, actiondata.squared );
    Tensor trot=rot.transpose();
    double prefactor = 1 / static_cast<double>( natoms );
    Vector v1;
    v1.zero();
    for(unsigned n=0; n<natoms; n++) {
      v1+=prefactor*matmul(trot, Vector(fdata.force[1+n],fdata.force[1+natoms+n],fdata.force[1+2*natoms+n]) );
    }
    for(unsigned n=0; n<natoms; n++) {
      Vector ff(fdata.force[1+n],fdata.force[1+natoms+n],fdata.force[1+2*natoms+n]);
      Vector oforce = actiondata.sqrtdisplace[n]*( matmul(trot,ff) - v1 );
      forces.thread_safe[n] += oforce[0];
      forces.thread_safe[natoms + n] += oforce[1];
      forces.thread_safe[2*natoms + n] += oforce[2];
    }
    for(unsigned a=0; a<3; a++) {
      for(unsigned b=0; b<3; b++) {
        double tmp1=0.;
        for(unsigned m=0; m<natoms; m++) {
          tmp1+=centeredpos[m][b]*fdata.force[1+a*natoms + m];
        }
        for(unsigned i=0; i<natoms; i++) {
          forces.thread_safe[i] += actiondata.sqrtdisplace[i]*tmp1*DRotDPos[a][b][i][0];
          forces.thread_safe[natoms + i] += actiondata.sqrtdisplace[i]*tmp1*DRotDPos[a][b][i][1];
          forces.thread_safe[2*natoms + i] += actiondata.sqrtdisplace[i]*tmp1*DRotDPos[a][b][i][2];
        }
      }
    }
  }
  double ff = fdata.force[0];
  for(unsigned j=0; j<fdata.deriv[0].size(); ++j ) {
    forces.thread_safe[j] += ff*fdata.deriv[0][j];
  }
}

void RMSDVector::applyNonZeroRankForces( std::vector<double>& outforces ) {
  taskmanager.applyForces( outforces );
}

}
}



