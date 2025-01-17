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
  keys.use("ARG");
  keys.setDisplayName("RMSD");
  keys.add("compulsory","TYPE","SIMPLE","the manner in which RMSD alignment is performed.  Should be OPTIMAL or SIMPLE.");
  keys.add("compulsory","ALIGN","1.0","the weights to use when aligning to the reference structure");
  keys.add("compulsory","DISPLACE","1.0","the weights to use when calculating the displacement from the reference structure");
  keys.addFlag("SQUARED",false," This should be set if you want mean squared displacement instead of RMSD ");
  keys.addFlag("UNORMALIZED",false,"by default the mean sequare deviation or root mean square deviation is calculated.  If this option is given no averaging is done");
  keys.addFlag("DISPLACEMENT",false,"Calculate the vector of displacements instead of the length of this vector");
  keys.addOutputComponent("disp","DISPLACEMENT","the vector of displacements for the atoms");
  keys.addOutputComponent("dist","DISPLACEMENT","the RMSD distance the atoms have moved");
  keys.setValueDescription("a vector containing the RMSD between the instantaneous structure and each of the reference structures that were input");
}

RMSDVector::RMSDVector(const ActionOptions&ao):
  Action(ao),
  ActionWithVector(ao),
  firststep(true) {
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

  unsigned natoms = getPntrToArgument(0)->getNumberOfValues() / 3;
  type.assign("SIMPLE");
  parse("TYPE",type);
  parseFlag("SQUARED",squared);
  align.resize( natoms );
  parseVector("ALIGN",align);
  displace.resize( natoms );
  parseVector("DISPLACE",displace);
  bool unorm=false;
  parseFlag("UNORMALIZED",unorm);
  norm_weights=!unorm;
  double wa=0, wd=0;
  sqrtdisplace.resize( displace.size() );
  for(unsigned i=0; i<align.size(); ++i) {
    wa+=align[i];
    wd+=displace[i];
  }

  if( wa>epsilon ) {
    double iwa = 1. / wa;
    for(unsigned i=0; i<align.size(); ++i) {
      align[i] *= iwa;
    }
  } else {
    double iwa = 1. / natoms;
    for(unsigned i=0; i<align.size(); ++i) {
      align[i] = iwa;
    }
  }
  if( wd>epsilon ) {
    if( !norm_weights ) {
      wd = 1;
    }
    double iwd = 1. / wd;
    for(unsigned i=0; i<align.size(); ++i) {
      displace[i] *= iwd;
    }
  } else {
    double iwd = 1. / natoms;
    for(unsigned i=0; i<align.size(); ++i) {
      displace[i] = iwd;
    }
  }
  for(unsigned i=0; i<align.size(); ++i) {
    sqrtdisplace[i] = sqrt(displace[i]);
  }

  parseFlag("DISPLACEMENT",displacement);
  if( displacement && (getPntrToArgument(1)->getRank()==1 || getPntrToArgument(1)->getShape()[0]<=1) ) {
    addComponentWithDerivatives("dist");
    componentIsNotPeriodic("dist");
    std::vector<unsigned> shape( 1, getPntrToArgument(0)->getNumberOfValues() );
    addComponent( "disp", shape );
    getPntrToComponent(1)->buildDataStore();
    componentIsNotPeriodic("disp");
  } else if( displacement ) {
    std::vector<unsigned> shape( 1, getPntrToArgument(1)->getShape()[0] );
    addComponent( "dist", shape );
    getPntrToComponent(0)->buildDataStore();
    componentIsNotPeriodic("dist");
    shape.resize(2);
    shape[0] = getPntrToArgument(1)->getShape()[0];
    shape[1] = getPntrToArgument(0)->getNumberOfValues();
    addComponent( "disp", shape );
    getPntrToComponent(1)->buildDataStore();
    getPntrToComponent(1)->reshapeMatrixStore( shape[1] );
    componentIsNotPeriodic("disp");
  } else if( (getPntrToArgument(1)->getRank()==1 || getPntrToArgument(1)->getShape()[0]==1) ) {
    addValue();
    setNotPeriodic();
  } else {
    std::vector<unsigned> shape( 1, getPntrToArgument(1)->getShape()[0] );
    addValue( shape );
    setNotPeriodic();
  }
  if( getPntrToArgument(1)->getRank()==1 || getPntrToArgument(1)->getNumberOfValues()==0 ) {
    myrmsd.resize(1);
  } else {
    myrmsd.resize( getPntrToArgument(1)->getShape()[0] );
  }

  if( getPntrToArgument(1)->getRank()==1 )
    log.printf("  calculating RMSD distance between %d atoms. Distance between the avectors of atoms in %s and %s\n",
               natoms, getPntrToArgument(0)->getName().c_str(), getPntrToArgument(1)->getName().c_str() );
  else
    log.printf("  calculating RMSD distance between %d sets of %d atoms. Distance between vector %s of atoms and matrix of configurations in %s\n",
               getPntrToArgument(1)->getShape()[0], natoms, getPntrToArgument(0)->getName().c_str(), getPntrToArgument(1)->getName().c_str() );
  log.printf("  method for alignment : %s \n",type.c_str() );
  if(squared) {
    log.printf("  chosen to use SQUARED option for MSD instead of RMSD\n");
  } else {
    log.printf("  using periodic boundary conditions\n");
  }
}

unsigned RMSDVector::getNumberOfDerivatives() {
  return getPntrToArgument(0)->getNumberOfValues() + getPntrToArgument(1)->getNumberOfValues();
}

void RMSDVector::setReferenceConfigurations() {
  unsigned natoms = getPntrToArgument(0)->getShape()[0] / 3;
  Vector center;
  std::vector<Vector> pos( natoms );
  for(unsigned jconf=0; jconf<myrmsd.size(); ++jconf) {
    center.zero();
    for(unsigned i=0; i<pos.size(); ++i) {
      for(unsigned j=0; j<3; ++j) {
        pos[i][j] = getPntrToArgument(1)->get( (3*jconf+j)*pos.size() + i );
      }
      center+=pos[i]*align[i];
    }
    for(unsigned i=0; i<pos.size(); ++i) {
      pos[i] -= center;
    }
    myrmsd[jconf].clear();
    myrmsd[jconf].set(align,displace,pos,type,true,norm_weights);
  }
}

double RMSDVector::calculateRMSD( const unsigned& current, std::vector<Vector>& pos, std::vector<Vector>& der, std::vector<Vector>& direction ) const {
  unsigned natoms = pos.size();
  for(unsigned i=0; i<natoms; ++i) {
    for(unsigned j=0; j<3; ++j) {
      pos[i][j] = getPntrToArgument(0)->get( j*natoms + i );
    }
  }

  if( displacement && type=="SIMPLE" ) {
    const Value* myval = getConstPntrToComponent(1);
    double r = myrmsd[current].simpleAlignment( align, displace, pos, myrmsd[current].getReference(), der, direction, squared );
    if( !doNotCalculateDerivatives() && myval->forcesWereAdded() ) {
      Vector comforce;
      comforce.zero();
      for(unsigned i=0; i<natoms; i++) {
        for(unsigned k=0; k<3; ++k) {
          comforce[k] += align[i]*myval->getForce( (3*current+k)*natoms + i);
        }
      }
      for(unsigned i=0; i<natoms; i++) {
        for(unsigned k=0; k<3; ++k) {
          direction[i][k] = myval->getForce( (3*current+k)*natoms + i ) - comforce[k];
        }
      }
    }
    return r;
  } else if( displacement ) {
    const Value* myval = getConstPntrToComponent(1);
    Tensor rot;
    Matrix<std::vector<Vector> > DRotDPos(3,3);
    std::vector<Vector> centeredpos( natoms ), centeredreference( natoms );
    double r = myrmsd[current].calc_PCAelements( pos, der, rot, DRotDPos, direction, centeredpos, centeredreference, squared );
    std::vector<Vector> ref( myrmsd[current].getReference() );
    if( !doNotCalculateDerivatives() && myval->forcesWereAdded() ) {
      Tensor trot=rot.transpose();
      double prefactor = 1 / static_cast<double>( natoms );
      Vector v1;
      v1.zero();
      for(unsigned n=0; n<natoms; n++) {
        Vector ff;
        for(unsigned k=0; k<3; ++k ) {
          ff[k] = myval->getForce( (3*current+k)*natoms + n );
        }
        v1+=prefactor*matmul(trot,ff);
      }
      // Notice that we use centreredreference here to accumulate the true forces
      for(unsigned n=0; n<natoms; n++) {
        Vector ff;
        for(unsigned k=0; k<3; ++k ) {
          ff[k] = myval->getForce( (3*current+k)*natoms + n );
        }
        centeredreference[n] = sqrtdisplace[n]*( matmul(trot,ff) - v1 );
      }
      for(unsigned a=0; a<3; a++) {
        for(unsigned b=0; b<3; b++) {
          double tmp1=0.;
          for(unsigned m=0; m<natoms; m++) {
            tmp1+=centeredpos[m][b]*myval->getForce( (3*current+a)*natoms + m );
          }
          for(unsigned i=0; i<natoms; i++) {
            centeredreference[i] += sqrtdisplace[i]*tmp1*DRotDPos[a][b][i];
          }
        }
      }
      // Now subtract the current force and add on the true force
      for(unsigned n=0; n<natoms; n++) {
        for(unsigned k=0; k<3; ++k) {
          direction[n][k] = centeredreference[n][k];
        }
      }
    } else {
      for(unsigned i=0; i<direction.size(); ++i) {
        direction[i] = sqrtdisplace[i]*( direction[i] - ref[i] );
      }
    }
    return r;
  }
  return myrmsd[current].calculate( pos, der, squared );
}

// calculator
void RMSDVector::calculate() {
  if( firststep || !getPntrToArgument(1)->isConstant() ) {
    setReferenceConfigurations();
    firststep=false;
  }

  if( getPntrToComponent(0)->getRank()==0 ) {
    unsigned natoms = getPntrToArgument(0)->getShape()[0] / 3;
    std::vector<Vector> pos( natoms ), der( natoms ), direction( natoms );
    double r = calculateRMSD( 0, pos, der, direction );

    getPntrToComponent(0)->set( r );
    if( getNumberOfComponents()==2 ) {
      Value* mydisp = getPntrToComponent(1);
      for(unsigned i=0; i<natoms; i++) {
        for(unsigned j=0; j<3; ++j ) {
          mydisp->set( j*natoms+i, direction[i][j] );
        }
      }
    }
    if( doNotCalculateDerivatives() ) {
      return;
    }

    Value* myval = getPntrToComponent(0);
    for(unsigned i=0; i<natoms; i++) {
      for(unsigned j=0; j<3; ++j ) {
        myval->setDerivative( j*natoms+i, der[i][j] );
      }
    }
  } else {
    runAllTasks();
  }
}

bool RMSDVector::checkForTaskForce( const unsigned& itask, const Value* myval ) const {
  if( myval->getRank()<2 ) {
    return ActionWithVector::checkForTaskForce( itask, myval );
  }
  unsigned nelements = myval->getShape()[1], startr = itask*nelements;
  for(unsigned j=0; j<nelements; ++j ) {
    if( fabs( myval->getForce( startr + j ) )>epsilon ) {
      return true;
    }
  }
  return false;
}

void RMSDVector::apply() {
  if( doNotCalculateDerivatives() ) {
    return;
  }

  if( getPntrToComponent(0)->getRank()==0 ) {
    std::vector<double> forces( getNumberOfDerivatives(), 0 );
    bool wasforced = getPntrToComponent(0)->applyForce( forces );

    if( getNumberOfComponents()==2 && getPntrToComponent(1)->forcesWereAdded() ) {
      unsigned natoms = getPntrToArgument(0)->getShape()[0] / 3;
      std::vector<Vector> pos( natoms ), der( natoms ), direction( natoms );
      double r = calculateRMSD( 0, pos, der, direction );
      for(unsigned i=0; i<natoms; ++i) {
        for(unsigned j=0; j<3; ++j ) {
          forces[j*natoms+i] += direction[i][j];
        }
      }
      wasforced=true;
    }
    if( wasforced ) {
      unsigned ss=0;
      addForcesOnArguments( 0, forces, ss, getLabel() );
    }
  } else {
    ActionWithVector::apply();
  }
}

void RMSDVector::performTask( const unsigned& current, MultiValue& myvals ) const {
  unsigned natoms = getPntrToArgument(0)->getShape()[0] / 3;
  std::vector<Vector>& pos( myvals.getFirstAtomVector() );
  std::vector<std::vector<Vector> > & allder( myvals.getFirstAtomDerivativeVector() );
  if( allder.size()!=2 ) {
    allder.resize(2);
  }
  std::vector<Vector>& der( allder[0] );
  std::vector<Vector>& direction( allder[1] );
  if( pos.size()!=natoms ) {
    pos.resize( natoms );
    der.resize( natoms );
    direction.resize( natoms );
  }
  for(unsigned i=0; i<pos.size(); ++i) {
    for(unsigned j=0; j<3; ++j) {
      pos[i][j] = getPntrToArgument(0)->get( j*natoms + i );
    }
  }
  double r = calculateRMSD( current, pos, der, direction );
  unsigned ostrn = getConstPntrToComponent(0)->getPositionInStream();
  myvals.setValue( ostrn, r );

  if( doNotCalculateDerivatives() ) {
    return;
  }

  for(unsigned i=0; i<natoms; i++) {
    for(unsigned j=0; j<3; ++j ) {
      myvals.addDerivative( ostrn, j*natoms+i, der[i][j] );
      myvals.updateIndex( ostrn, j*natoms+i );
    }
  }
}

void RMSDVector::gatherStoredValue( const unsigned& valindex, const unsigned& code, const MultiValue& myvals,
                                    const unsigned& bufstart, std::vector<double>& buffer ) const {
  if( getConstPntrToComponent(valindex)->getRank()==1 ) {
    ActionWithVector::gatherStoredValue( valindex, code, myvals, bufstart, buffer );
    return;
  }
  const std::vector<Vector>& direction( myvals.getConstFirstAtomDerivativeVector()[1] );
  unsigned natoms = direction.size();
  unsigned vindex = bufstart + 3*code*natoms;
  for(unsigned i=0; i<natoms; ++i) {
    for(unsigned j=0; j<3; ++j ) {
      buffer[vindex + j*natoms + i] += direction[i][j];
    }
  }
}

void RMSDVector::gatherForcesOnStoredValue( const Value* myval, const unsigned& itask, const MultiValue& myvals, std::vector<double>& forces ) const {
  if( myval->getRank()==1 ) {
    ActionWithVector::gatherForcesOnStoredValue( myval, itask, myvals, forces );
    return;
  }
  const std::vector<Vector>& direction( myvals.getConstFirstAtomDerivativeVector()[1] );
  unsigned natoms = direction.size();
  for(unsigned i=0; i<natoms; ++i) {
    for(unsigned j=0; j<3; ++j ) {
      forces[j*natoms+i] += direction[i][j];
    }
  }
}


}
}



