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
#include "PDB.h"
#include "Log.h"
#include "Exception.h"
#include <cmath>
#include <iostream>
#include "Tools.h"

namespace PLMD {

RMSD::RMSD() : alignmentMethod(SIMPLE),reference_center_is_calculated(false),reference_center_is_removed(false),positions_center_is_calculated(false),positions_center_is_removed(false) {}

///
/// general method to set all the rmsd property at once by using a pdb where occupancy column sets the weights for the atoms involved in the
/// alignment and beta sets the weight that are used for calculating the displacement.
///
void RMSD::set(const PDB&pdb, const std::string & mytype, bool remove_center, bool normalize_weights ) {

  set(pdb.getOccupancy(),pdb.getBeta(),pdb.getPositions(),mytype,remove_center,normalize_weights);

}
void RMSD::set(const std::vector<double> & align, const std::vector<double> & displace, const std::vector<Vector> & reference, const std::string & mytype, bool remove_center, bool normalize_weights ) {

  setReference(reference); // this by default remove the com and assumes uniform weights
  setAlign(align, normalize_weights, remove_center); // this recalculates the com with weights. If remove_center=false then it restore the center back
  setDisplace(displace, normalize_weights);  // this is does not affect any calculation of the weights
  setType(mytype);

}

void RMSD::setType(const std::string & mytype) {

  alignmentMethod=SIMPLE; // initialize with the simplest case: no rotation
  if (mytype=="SIMPLE") {
    alignmentMethod=SIMPLE;
  }
  else if (mytype=="OPTIMAL") {
    alignmentMethod=OPTIMAL;
  }
  else if (mytype=="OPTIMAL-FAST") {
    alignmentMethod=OPTIMAL_FAST;
  }
  else plumed_merror("unknown RMSD type" + mytype);

}

void RMSD::clear() {
  reference.clear();
  reference_center.zero();
  reference_center_is_calculated=false;
  reference_center_is_removed=false;
  align.clear();
  displace.clear();
  positions_center.zero();
  positions_center_is_calculated=false;
  positions_center_is_removed=false;
}

std::string RMSD::getMethod() {
  std::string mystring;
  switch(alignmentMethod) {
  case SIMPLE: mystring.assign("SIMPLE"); break;
  case OPTIMAL: mystring.assign("OPTIMAL"); break;
  case OPTIMAL_FAST: mystring.assign("OPTIMAL-FAST"); break;
  }
  return mystring;
}
///
/// this calculates the center of mass for the reference and removes it from the reference itself
/// considering uniform weights for alignment
///
void RMSD::setReference(const std::vector<Vector> & reference) {
  unsigned n=reference.size();
  this->reference=reference;
  plumed_massert(align.empty(),"you should first clear() an RMSD object, then set a new reference");
  plumed_massert(displace.empty(),"you should first clear() an RMSD object, then set a new reference");
  align.resize(n,1.0/n);
  displace.resize(n,1.0/n);
  for(unsigned i=0; i<n; i++) reference_center+=this->reference[i]*align[i];
  #pragma omp simd
  for(unsigned i=0; i<n; i++) this->reference[i]-=reference_center;
  reference_center_is_calculated=true;
  reference_center_is_removed=true;
}
std::vector<Vector> RMSD::getReference() {
  return reference;
}
///
/// the alignment weights are here normalized to 1 and  the center of the reference is removed accordingly
///
void RMSD::setAlign(const std::vector<double> & align, bool normalize_weights, bool remove_center) {
  unsigned n=reference.size();
  plumed_massert(this->align.size()==align.size(),"mismatch in dimension of align/displace arrays");
  this->align=align;
  if(normalize_weights) {
    double w=0.0;
    #pragma omp simd reduction(+:w)
    for(unsigned i=0; i<n; i++) w+=this->align[i];
    if(w>epsilon) {
      double inv=1.0/w;
      #pragma omp simd
      for(unsigned i=0; i<n; i++) this->align[i]*=inv;
    } else {
      double inv=1.0/n;
      #pragma omp simd
      for(unsigned i=0; i<n; i++) this->align[i]=inv;
    }
  }
  // recalculate the center anyway
  // just remove the center if that is asked
  // if the center was removed before, then add it and store the new one
  if(reference_center_is_removed) {
    plumed_massert(reference_center_is_calculated," seems that the reference center has been removed but not calculated and stored!");
    addCenter(reference,reference_center);
  }
  reference_center=calculateCenter(reference,this->align);
  reference_center_is_calculated=true;
  if(remove_center) {
    removeCenter(reference,reference_center);
    reference_center_is_removed=true;
  } else {
    reference_center_is_removed=false;
  }
}
std::vector<double> RMSD::getAlign() {
  return align;
}
///
/// here the weigth for normalized weighths are normalized and set
///
void RMSD::setDisplace(const std::vector<double> & displace, bool normalize_weights) {
  unsigned n=reference.size();
  plumed_massert(this->displace.size()==displace.size(),"mismatch in dimension of align/displace arrays");
  this->displace=displace;
  if(normalize_weights) {
    double w=0.0;
    #pragma omp simd reduction(+:w)
    for(unsigned i=0; i<n; i++) w+=this->displace[i];
    if(w>epsilon) {
      double inv=1.0/w;
      #pragma omp simd
      for(unsigned i=0; i<n; i++) this->displace[i]*=inv;
    } else {
      double inv=1.0/n;
      #pragma omp simd
      for(unsigned i=0; i<n; i++) this->displace[i]=inv;
    }
  }
}
std::vector<double> RMSD::getDisplace() {
  return displace;
}
///
/// This is the main workhorse for rmsd that decides to use specific optimal alignment versions
///
double RMSD::calculate(const std::vector<Vector> & positions,std::vector<Vector> &derivatives, bool squared)const {

  double ret=0.;

  switch(alignmentMethod) {
  case SIMPLE : {
    //	do a simple alignment without rotation
    std::vector<Vector> displacement( derivatives.size() );
    ret=simpleAlignment(align,displace,positions,reference,derivatives,displacement,squared);
    break;
  } case OPTIMAL_FAST : {
    // this is calling the fastest option:
    if(align==displace) ret=optimalAlignment<false,true>(align,displace,positions,reference,derivatives,squared);
    else                ret=optimalAlignment<false,false>(align,displace,positions,reference,derivatives,squared);
    break;

  } case OPTIMAL : {
    // this is the fast routine but in the "safe" mode, which gives less numerical error:
    if(align==displace) ret=optimalAlignment<true,true>(align,displace,positions,reference,derivatives,squared);
    else ret=optimalAlignment<true,false>(align,displace,positions,reference,derivatives,squared);
    break;
  }
  }

  return ret;

}


/// convenience method for calculating the standard derivatives and the derivative of the rmsd respect to the reference position
double RMSD::calc_DDistDRef( const std::vector<Vector>& positions, std::vector<Vector> &derivatives, std::vector<Vector>& DDistDRef, const bool squared  ) {
  double ret=0.;
  switch(alignmentMethod) {
  case SIMPLE:
    plumed_merror("derivative of the refreence frame not implemented for SIMPLE alignmentMethod \n");
    break;
  case OPTIMAL_FAST:
    if(align==displace) ret=optimalAlignment_DDistDRef<false,true>(align,displace,positions,reference,derivatives,DDistDRef, squared);
    else                ret=optimalAlignment_DDistDRef<false,false>(align,displace,positions,reference,derivatives,DDistDRef,squared);
    break;
  case OPTIMAL:
    if(align==displace) ret=optimalAlignment_DDistDRef<true,true>(align,displace,positions,reference,derivatives,DDistDRef,squared);
    else                ret=optimalAlignment_DDistDRef<true,false>(align,displace,positions,reference,derivatives,DDistDRef,squared);
    break;
  }
  return ret;

}

/// convenience method for calculating the standard derivatives and the derivative of the rmsd respect to the reference position without the matrix contribution
/// as required by SOMA
double RMSD::calc_SOMA( const std::vector<Vector>& positions, std::vector<Vector> &derivatives, std::vector<Vector>& DDistDRef, const bool squared  ) {
  double ret=0.;
  switch(alignmentMethod) {
  case SIMPLE:
    plumed_merror("derivative of the refreence frame not implemented for SIMPLE alignmentMethod \n");
    break;
  case OPTIMAL_FAST:
    if(align==displace) ret=optimalAlignment_SOMA<false,true>(align,displace,positions,reference,derivatives,DDistDRef, squared);
    else                ret=optimalAlignment_SOMA<false,false>(align,displace,positions,reference,derivatives,DDistDRef,squared);
    break;
  case OPTIMAL:
    if(align==displace) ret=optimalAlignment_SOMA<true,true>(align,displace,positions,reference,derivatives,DDistDRef,squared);
    else                ret=optimalAlignment_SOMA<true,false>(align,displace,positions,reference,derivatives,DDistDRef,squared);
    break;
  }
  return ret;

}

double RMSD::calc_DDistDRef_Rot_DRotDPos( const std::vector<Vector>& positions, std::vector<Vector> &derivatives, std::vector<Vector>& DDistDRef, Tensor & Rot, Matrix<std::vector<Vector> > &DRotDPos, const bool squared  ) {
  double ret=0.;
  switch(alignmentMethod) {
  case SIMPLE:
    plumed_merror("derivative of the refreence frame not implemented for SIMPLE alignmentMethod \n");
    break;
  case OPTIMAL_FAST:
    if(align==displace) ret=optimalAlignment_DDistDRef_Rot_DRotDPos<false,true>(align,displace,positions,reference,derivatives,DDistDRef, Rot, DRotDPos,  squared);
    else                ret=optimalAlignment_DDistDRef_Rot_DRotDPos<false,false>(align,displace,positions,reference,derivatives,DDistDRef, Rot, DRotDPos, squared);
    break;
  case OPTIMAL:
    if(align==displace) ret=optimalAlignment_DDistDRef_Rot_DRotDPos<true,true>(align,displace,positions,reference,derivatives,DDistDRef, Rot, DRotDPos, squared);
    else                ret=optimalAlignment_DDistDRef_Rot_DRotDPos<true,false>(align,displace,positions,reference,derivatives,DDistDRef, Rot, DRotDPos, squared);
    break;
  }
  return ret;
}

double RMSD::calc_DDistDRef_Rot_DRotDPos_DRotDRef( const std::vector<Vector>& positions, std::vector<Vector> &derivatives, std::vector<Vector>& DDistDRef, Tensor & Rot, Matrix<std::vector<Vector> > &DRotDPos,  Matrix<std::vector<Vector> > &DRotDRef, const bool squared  ) {
  double ret=0.;
  switch(alignmentMethod) {
  case SIMPLE:
    plumed_merror("derivative of the refreence frame not implemented for SIMPLE alignmentMethod \n");
    break;
  case OPTIMAL_FAST:
    if(align==displace) ret=optimalAlignment_DDistDRef_Rot_DRotDPos_DRotDRef<false,true>(align,displace,positions,reference,derivatives,DDistDRef, Rot, DRotDPos, DRotDRef,   squared);
    else                ret=optimalAlignment_DDistDRef_Rot_DRotDPos_DRotDRef<false,false>(align,displace,positions,reference,derivatives,DDistDRef, Rot, DRotDPos, DRotDRef,  squared);
    break;
  case OPTIMAL:
    if(align==displace) ret=optimalAlignment_DDistDRef_Rot_DRotDPos_DRotDRef<true,true>(align,displace,positions,reference,derivatives,DDistDRef, Rot, DRotDPos, DRotDRef, squared);
    else                ret=optimalAlignment_DDistDRef_Rot_DRotDPos_DRotDRef<true,false>(align,displace,positions,reference,derivatives,DDistDRef, Rot, DRotDPos, DRotDRef, squared);
    break;
  }
  return ret;
}

double RMSD::calc_Rot_DRotDRr01( const std::vector<Vector>& positions, Tensor & Rotation, std::array<std::array<Tensor,3>,3> & DRotDRr01, const bool squared) {
  double ret=0.;
  switch(alignmentMethod) {
  case SIMPLE:
    plumed_merror("derivative of the refreence frame not implemented for SIMPLE alignmentMethod \n");
    break;
  case OPTIMAL_FAST:
    if(align==displace) ret=optimalAlignment_Rot_DRotDRr01<false,true>(align,displace,positions,reference, Rotation, DRotDRr01,   squared);
    else                ret=optimalAlignment_Rot_DRotDRr01<false,false>(align,displace,positions,reference, Rotation, DRotDRr01,  squared);
    break;
  case OPTIMAL:
    if(align==displace) ret=optimalAlignment_Rot_DRotDRr01<true,true>(align,displace,positions,reference, Rotation, DRotDRr01, squared);
    else                ret=optimalAlignment_Rot_DRotDRr01<true,false>(align,displace,positions,reference, Rotation, DRotDRr01, squared);
    break;
  }
  return ret;
}

double RMSD::calc_Rot( const std::vector<Vector>& positions, std::vector<Vector> &derivatives, Tensor & Rotation, const bool squared) {
  double ret=0.;
  switch(alignmentMethod) {
  case SIMPLE:
    plumed_merror("derivative of the refreence frame not implemented for SIMPLE alignmentMethod \n");
    break;
  case OPTIMAL_FAST:
    if(align==displace) ret=optimalAlignment_Rot<false,true>(align,displace,positions,reference,derivatives, Rotation, squared);
    else                ret=optimalAlignment_Rot<false,false>(align,displace,positions,reference,derivatives, Rotation, squared);
    break;
  case OPTIMAL:
    if(align==displace) ret=optimalAlignment_Rot<true,true>(align,displace,positions,reference,derivatives, Rotation, squared);
    else                ret=optimalAlignment_Rot<true,false>(align,displace,positions,reference,derivatives, Rotation, squared);
    break;
  }
  return ret;
}

double RMSD::calculateWithCloseStructure( const std::vector<Vector>& positions, std::vector<Vector> &derivatives, Tensor & rotationPosClose, Tensor & rotationRefClose, std::array<std::array<Tensor,3>,3> & drotationPosCloseDrr01, const bool squared) {
  double ret=0.;
  switch(alignmentMethod) {
  case SIMPLE:
    plumed_merror("derivative of the refreence frame not implemented for SIMPLE alignmentMethod \n");
    break;
  case OPTIMAL_FAST:
    if(align==displace) ret=optimalAlignmentWithCloseStructure<false,true>(align,displace,positions,reference,derivatives, rotationPosClose, rotationRefClose, drotationPosCloseDrr01, squared);
    else                ret=optimalAlignmentWithCloseStructure<false,false>(align,displace,positions,reference,derivatives, rotationPosClose, rotationRefClose, drotationPosCloseDrr01, squared);
    break;
  case OPTIMAL:
    if(align==displace) ret=optimalAlignmentWithCloseStructure<true,true>(align,displace,positions,reference,derivatives, rotationPosClose, rotationRefClose, drotationPosCloseDrr01, squared);
    else                ret=optimalAlignmentWithCloseStructure<true,false>(align,displace,positions,reference,derivatives, rotationPosClose, rotationRefClose, drotationPosCloseDrr01, squared);
    break;
  }
  return ret;
}

double RMSD::calc_PCAelements( const std::vector<Vector>& positions, std::vector<Vector> &DDistDPos, Tensor & Rotation, Matrix<std::vector<Vector> > & DRotDPos,std::vector<Vector>  & alignedpositions, std::vector<Vector> & centeredpositions, std::vector<Vector> &centeredreference, const bool& squared  ) const {
  double ret=0.;
  switch(alignmentMethod) {
  case SIMPLE:
    plumed_merror("derivative of the refreence frame not implemented for SIMPLE alignmentMethod \n");
    break;
  case OPTIMAL_FAST:
    if(align==displace) ret=optimalAlignment_PCA<false,true>(align,displace,positions,reference, alignedpositions, centeredpositions,centeredreference,Rotation,DDistDPos,DRotDPos,squared);
    else                ret=optimalAlignment_PCA<false,false>(align,displace,positions,reference, alignedpositions, centeredpositions,centeredreference,Rotation,DDistDPos,DRotDPos,squared);
    break;
  case OPTIMAL:
    if(align==displace) ret=optimalAlignment_PCA<true,true>(align,displace,positions,reference, alignedpositions, centeredpositions,centeredreference,Rotation,DDistDPos,DRotDPos,squared);
    else                ret=optimalAlignment_PCA<true,false>(align,displace,positions,reference, alignedpositions, centeredpositions,centeredreference,Rotation,DDistDPos,DRotDPos,squared);
    break;
  }
  return ret;
}


double RMSD::calc_FitElements( const std::vector<Vector>& positions, Tensor & Rotation, Matrix<std::vector<Vector> > & DRotDPos, std::vector<Vector> & centeredpositions, Vector &center_positions, const bool& squared  ) {
  double ret=0.;
  switch(alignmentMethod) {
  case SIMPLE:
    plumed_merror("derivative of the refreence frame not implemented for SIMPLE alignmentMethod \n");
    break;
  case OPTIMAL_FAST:
    if(align==displace)ret=optimalAlignment_Fit<false,true>(align,displace,positions,reference, Rotation,DRotDPos,centeredpositions,center_positions,squared);
    else               ret=optimalAlignment_Fit<false,false>(align,displace,positions,reference, Rotation,DRotDPos,centeredpositions,center_positions,squared);
    break;
  case OPTIMAL:
    if(align==displace)ret=optimalAlignment_Fit<true,true>(align,displace,positions,reference,Rotation,DRotDPos,centeredpositions,center_positions,squared);
    else               ret=optimalAlignment_Fit<true,false>(align,displace,positions,reference,Rotation,DRotDPos,centeredpositions,center_positions,squared);
    break;
  }
  return ret;
}






double RMSD::simpleAlignment(const  std::vector<double>  & align,
                             const  std::vector<double>  & displace,
                             const std::vector<Vector> & positions,
                             const std::vector<Vector> & reference,
                             std::vector<Vector>  & derivatives,
                             std::vector<Vector>  & displacement,
                             bool squared)const {

  double dist(0);
  unsigned n=reference.size();

  Vector apositions;
  Vector areference;
  Vector dpositions;
  Vector dreference;

  for(unsigned i=0; i<n; i++) {
    double aw=align[i];
    double dw=displace[i];
    apositions+=positions[i]*aw;
    areference+=reference[i]*aw;
    dpositions+=positions[i]*dw;
    dreference+=reference[i]*dw;
  }

  Vector shift=((apositions-areference)-(dpositions-dreference));
  for(unsigned i=0; i<n; i++) {
    displacement[i]=(positions[i]-apositions)-(reference[i]-areference);
    dist+=displace[i]*displacement[i].modulo2();
    derivatives[i]=2*(displace[i]*displacement[i]+align[i]*shift);
  }

  if(!squared) {
    // sqrt
    dist=std::sqrt(dist);
    ///// sqrt on derivatives
    for(unsigned i=0; i<n; i++) {derivatives[i]*=(0.5/dist);}
  }
  return dist;
}

// this below enable the standard case for rmsd where the rmsd is calculated and the derivative of rmsd respect to positions is retrieved
// additionally this assumes that the com of the reference is already subtracted.
#define OLDRMSD
#ifdef OLDRMSD
// notice that in the current implementation the safe argument only makes sense for
// align==displace
template <bool safe,bool alEqDis>
double RMSD::optimalAlignment(const  std::vector<double>  & align,
                              const  std::vector<double>  & displace,
                              const std::vector<Vector> & positions,
                              const std::vector<Vector> & reference,
                              std::vector<Vector>  & derivatives, bool squared)const {
  const unsigned n=reference.size();
// This is the trace of positions*positions + reference*reference
  double rr00(0);
  double rr11(0);
// This is positions*reference
  Tensor rr01;

  derivatives.resize(n);

  Vector cpositions;

// first expensive loop: compute centers
  for(unsigned iat=0; iat<n; iat++) {
    double w=align[iat];
    cpositions+=positions[iat]*w;
  }

// second expensive loop: compute second moments wrt centers
  for(unsigned iat=0; iat<n; iat++) {
    double w=align[iat];
    rr00+=dotProduct(positions[iat]-cpositions,positions[iat]-cpositions)*w;
    rr11+=dotProduct(reference[iat],reference[iat])*w;
    rr01+=Tensor(positions[iat]-cpositions,reference[iat])*w;
  }

  Tensor4d m;

  m[0][0]=2.0*(-rr01[0][0]-rr01[1][1]-rr01[2][2]);
  m[1][1]=2.0*(-rr01[0][0]+rr01[1][1]+rr01[2][2]);
  m[2][2]=2.0*(+rr01[0][0]-rr01[1][1]+rr01[2][2]);
  m[3][3]=2.0*(+rr01[0][0]+rr01[1][1]-rr01[2][2]);
  m[0][1]=2.0*(-rr01[1][2]+rr01[2][1]);
  m[0][2]=2.0*(+rr01[0][2]-rr01[2][0]);
  m[0][3]=2.0*(-rr01[0][1]+rr01[1][0]);
  m[1][2]=2.0*(-rr01[0][1]-rr01[1][0]);
  m[1][3]=2.0*(-rr01[0][2]-rr01[2][0]);
  m[2][3]=2.0*(-rr01[1][2]-rr01[2][1]);
  m[1][0] = m[0][1];
  m[2][0] = m[0][2];
  m[2][1] = m[1][2];
  m[3][0] = m[0][3];
  m[3][1] = m[1][3];
  m[3][2] = m[2][3];

  Tensor dm_drr01[4][4];
  if(!alEqDis) {
    dm_drr01[0][0] = 2.0*Tensor(-1.0, 0.0, 0.0,  0.0,-1.0, 0.0,  0.0, 0.0,-1.0);
    dm_drr01[1][1] = 2.0*Tensor(-1.0, 0.0, 0.0,  0.0,+1.0, 0.0,  0.0, 0.0,+1.0);
    dm_drr01[2][2] = 2.0*Tensor(+1.0, 0.0, 0.0,  0.0,-1.0, 0.0,  0.0, 0.0,+1.0);
    dm_drr01[3][3] = 2.0*Tensor(+1.0, 0.0, 0.0,  0.0,+1.0, 0.0,  0.0, 0.0,-1.0);
    dm_drr01[0][1] = 2.0*Tensor( 0.0, 0.0, 0.0,  0.0, 0.0,-1.0,  0.0,+1.0, 0.0);
    dm_drr01[0][2] = 2.0*Tensor( 0.0, 0.0,+1.0,  0.0, 0.0, 0.0, -1.0, 0.0, 0.0);
    dm_drr01[0][3] = 2.0*Tensor( 0.0,-1.0, 0.0, +1.0, 0.0, 0.0,  0.0, 0.0, 0.0);
    dm_drr01[1][2] = 2.0*Tensor( 0.0,-1.0, 0.0, -1.0, 0.0, 0.0,  0.0, 0.0, 0.0);
    dm_drr01[1][3] = 2.0*Tensor( 0.0, 0.0,-1.0,  0.0, 0.0, 0.0, -1.0, 0.0, 0.0);
    dm_drr01[2][3] = 2.0*Tensor( 0.0, 0.0, 0.0,  0.0, 0.0,-1.0,  0.0,-1.0, 0.0);
    dm_drr01[1][0] = dm_drr01[0][1];
    dm_drr01[2][0] = dm_drr01[0][2];
    dm_drr01[2][1] = dm_drr01[1][2];
    dm_drr01[3][0] = dm_drr01[0][3];
    dm_drr01[3][1] = dm_drr01[1][3];
    dm_drr01[3][2] = dm_drr01[2][3];
  }

  double dist=0.0;
  Vector4d q;

  Tensor dq_drr01[4];
  if(!alEqDis) {
    Vector4d eigenvals;
    Tensor4d eigenvecs;
    diagMatSym(m, eigenvals, eigenvecs );
    dist=eigenvals[0]+rr00+rr11;
    q=Vector4d(eigenvecs[0][0],eigenvecs[0][1],eigenvecs[0][2],eigenvecs[0][3]);
    double dq_dm[4][4][4];
    for(unsigned i=0; i<4; i++) for(unsigned j=0; j<4; j++) for(unsigned k=0; k<4; k++) {
          double tmp=0.0;
// perturbation theory for matrix m
          for(unsigned l=1; l<4; l++) tmp+=eigenvecs[l][j]*eigenvecs[l][i]/(eigenvals[0]-eigenvals[l])*eigenvecs[0][k];
          dq_dm[i][j][k]=tmp;
        }
// propagation to _drr01
    for(unsigned i=0; i<4; i++) {
      Tensor tmp;
      for(unsigned j=0; j<4; j++) for(unsigned k=0; k<4; k++) {
          tmp+=dq_dm[i][j][k]*dm_drr01[j][k];
        }
      dq_drr01[i]=tmp;
    }
  } else {
    VectorGeneric<1> eigenvals;
    TensorGeneric<1,4> eigenvecs;
    diagMatSym(m, eigenvals, eigenvecs );
    dist=eigenvals[0]+rr00+rr11;
    q=Vector4d(eigenvecs[0][0],eigenvecs[0][1],eigenvecs[0][2],eigenvecs[0][3]);
  }


// This is the rotation matrix that brings reference to positions
// i.e. matmul(rotation,reference[iat])+shift is fitted to positions[iat]

  Tensor rotation;
  rotation[0][0]=q[0]*q[0]+q[1]*q[1]-q[2]*q[2]-q[3]*q[3];
  rotation[1][1]=q[0]*q[0]-q[1]*q[1]+q[2]*q[2]-q[3]*q[3];
  rotation[2][2]=q[0]*q[0]-q[1]*q[1]-q[2]*q[2]+q[3]*q[3];
  rotation[0][1]=2*(+q[0]*q[3]+q[1]*q[2]);
  rotation[0][2]=2*(-q[0]*q[2]+q[1]*q[3]);
  rotation[1][2]=2*(+q[0]*q[1]+q[2]*q[3]);
  rotation[1][0]=2*(-q[0]*q[3]+q[1]*q[2]);
  rotation[2][0]=2*(+q[0]*q[2]+q[1]*q[3]);
  rotation[2][1]=2*(-q[0]*q[1]+q[2]*q[3]);


  std::array<std::array<Tensor,3>,3> drotation_drr01;
  if(!alEqDis) {
    drotation_drr01[0][0]=2*q[0]*dq_drr01[0]+2*q[1]*dq_drr01[1]-2*q[2]*dq_drr01[2]-2*q[3]*dq_drr01[3];
    drotation_drr01[1][1]=2*q[0]*dq_drr01[0]-2*q[1]*dq_drr01[1]+2*q[2]*dq_drr01[2]-2*q[3]*dq_drr01[3];
    drotation_drr01[2][2]=2*q[0]*dq_drr01[0]-2*q[1]*dq_drr01[1]-2*q[2]*dq_drr01[2]+2*q[3]*dq_drr01[3];
    drotation_drr01[0][1]=2*(+(q[0]*dq_drr01[3]+dq_drr01[0]*q[3])+(q[1]*dq_drr01[2]+dq_drr01[1]*q[2]));
    drotation_drr01[0][2]=2*(-(q[0]*dq_drr01[2]+dq_drr01[0]*q[2])+(q[1]*dq_drr01[3]+dq_drr01[1]*q[3]));
    drotation_drr01[1][2]=2*(+(q[0]*dq_drr01[1]+dq_drr01[0]*q[1])+(q[2]*dq_drr01[3]+dq_drr01[2]*q[3]));
    drotation_drr01[1][0]=2*(-(q[0]*dq_drr01[3]+dq_drr01[0]*q[3])+(q[1]*dq_drr01[2]+dq_drr01[1]*q[2]));
    drotation_drr01[2][0]=2*(+(q[0]*dq_drr01[2]+dq_drr01[0]*q[2])+(q[1]*dq_drr01[3]+dq_drr01[1]*q[3]));
    drotation_drr01[2][1]=2*(-(q[0]*dq_drr01[1]+dq_drr01[0]*q[1])+(q[2]*dq_drr01[3]+dq_drr01[2]*q[3]));
  }

  double prefactor=2.0;

  if(!squared && alEqDis) prefactor*=0.5/std::sqrt(dist);

// if "safe", recompute dist here to a better accuracy
  if(safe || !alEqDis) dist=0.0;

// If safe is set to "false", MSD is taken from the eigenvalue of the M matrix
// If safe is set to "true", MSD is recomputed from the rotational matrix
// For some reason, this last approach leads to less numerical noise but adds an overhead

  Tensor ddist_drotation;
  Vector ddist_dcpositions;

// third expensive loop: derivatives
  for(unsigned iat=0; iat<n; iat++) {
    Vector d(positions[iat]-cpositions - matmul(rotation,reference[iat]));
    if(alEqDis) {
// there is no need for derivatives of rotation and shift here as it is by construction zero
// (similar to Hellman-Feynman forces)
      derivatives[iat]= prefactor*align[iat]*d;
      if(safe) dist+=align[iat]*modulo2(d);
    } else {
// the case for align != displace is different, sob:
      dist+=displace[iat]*modulo2(d);
// these are the derivatives assuming the roto-translation as frozen
      derivatives[iat]=2*displace[iat]*d;
// here I accumulate derivatives wrt rotation matrix ..
      ddist_drotation+=-2*displace[iat]*extProduct(d,reference[iat]);
// .. and cpositions
      ddist_dcpositions+=-2*displace[iat]*d;
    }
  }

  if(!alEqDis) {
    Tensor ddist_drr01;
    for(unsigned i=0; i<3; i++) for(unsigned j=0; j<3; j++) ddist_drr01+=ddist_drotation[i][j]*drotation_drr01[i][j];
    for(unsigned iat=0; iat<n; iat++) {
// this is propagating to positions.
// I am implicitly using the derivative of rr01 wrt positions here
      derivatives[iat]+=matmul(ddist_drr01,reference[iat])*align[iat];
      derivatives[iat]+=ddist_dcpositions*align[iat];
    }
  }
  if(!squared) {
    dist=std::sqrt(dist);
    if(!alEqDis) {
      double xx=0.5/dist;
      for(unsigned iat=0; iat<n; iat++) derivatives[iat]*=xx;
    }
  }

  return dist;
}
#else
/// note that this method is intended to be repeatedly invoked
/// when the reference does already have the center subtracted
/// but the position has not calculated center and not subtracted
template <bool safe,bool alEqDis>
double RMSD::optimalAlignment(const  std::vector<double>  & align,
                              const  std::vector<double>  & displace,
                              const std::vector<Vector> & positions,
                              const std::vector<Vector> & reference,
                              std::vector<Vector>  & derivatives,
                              bool squared) const {
  //std::cerr<<"setting up the core data \n";
  RMSDCoreData cd(align,displace,positions,reference);

  // transfer the settings for the center to let the CoreCalc deal with it
  cd.setPositionsCenterIsRemoved(positions_center_is_removed);
  if(positions_center_is_calculated) {cd.setPositionsCenter(positions_center);}
  else {cd.calcPositionsCenter();};

  cd.setReferenceCenterIsRemoved(reference_center_is_removed);
  if(!reference_center_is_calculated) {cd.calcReferenceCenter();}
  else {cd.setReferenceCenter(reference_center);}

  // Perform the diagonalization and all the needed stuff
  cd.doCoreCalc(safe,alEqDis);
  // make the core calc distance
  double dist=cd.getDistance(squared);
//  make the derivatives by using pieces calculated in coreCalc (probably the best is just to copy the vector...)
  derivatives=cd.getDDistanceDPositions();
  return dist;
}
#endif
template <bool safe,bool alEqDis>
double RMSD::optimalAlignment_DDistDRef(const  std::vector<double>  & align,
                                        const  std::vector<double>  & displace,
                                        const std::vector<Vector> & positions,
                                        const std::vector<Vector> & reference,
                                        std::vector<Vector>  & derivatives,
                                        std::vector<Vector> & ddistdref,
                                        bool squared) const {
  //initialize the data into the structure
  // typically the positions do not have the com neither calculated nor subtracted. This layer takes care of this business
  RMSDCoreData cd(align,displace,positions,reference);
  // transfer the settings for the center to let the CoreCalc deal with it
  // transfer the settings for the center to let the CoreCalc deal with it
  cd.setPositionsCenterIsRemoved(positions_center_is_removed);
  if(positions_center_is_calculated) {cd.setPositionsCenter(positions_center);}
  else {cd.calcPositionsCenter();};

  cd.setReferenceCenterIsRemoved(reference_center_is_removed);
  if(!reference_center_is_calculated) {cd.calcReferenceCenter();}
  else {cd.setReferenceCenter(reference_center);}

  // Perform the diagonalization and all the needed stuff
  cd.doCoreCalc(safe,alEqDis);
  // make the core calc distance
  double dist=cd.getDistance(squared);
//  make the derivatives by using pieces calculated in coreCalc (probably the best is just to copy the vector...)
  derivatives=cd.getDDistanceDPositions();
  ddistdref=cd.getDDistanceDReference();
  return dist;
}

template <bool safe,bool alEqDis>
double RMSD::optimalAlignment_SOMA(const  std::vector<double>  & align,
                                   const  std::vector<double>  & displace,
                                   const std::vector<Vector> & positions,
                                   const std::vector<Vector> & reference,
                                   std::vector<Vector>  & derivatives,
                                   std::vector<Vector> & ddistdref,
                                   bool squared) const {
  //initialize the data into the structure
  // typically the positions do not have the com neither calculated nor subtracted. This layer takes care of this business
  RMSDCoreData cd(align,displace,positions,reference);
  // transfer the settings for the center to let the CoreCalc deal with it
  // transfer the settings for the center to let the CoreCalc deal with it
  cd.setPositionsCenterIsRemoved(positions_center_is_removed);
  if(positions_center_is_calculated) {cd.setPositionsCenter(positions_center);}
  else {cd.calcPositionsCenter();};

  cd.setReferenceCenterIsRemoved(reference_center_is_removed);
  if(!reference_center_is_calculated) {cd.calcReferenceCenter();}
  else {cd.setReferenceCenter(reference_center);}

  // Perform the diagonalization and all the needed stuff
  cd.doCoreCalc(safe,alEqDis);
  // make the core calc distance
  double dist=cd.getDistance(squared);
//  make the derivatives by using pieces calculated in coreCalc (probably the best is just to copy the vector...)
  derivatives=cd.getDDistanceDPositions();
  ddistdref=cd.getDDistanceDReferenceSOMA();
  return dist;
}


template <bool safe,bool alEqDis>
double RMSD::optimalAlignment_DDistDRef_Rot_DRotDPos(const  std::vector<double>  & align,
    const  std::vector<double>  & displace,
    const std::vector<Vector> & positions,
    const std::vector<Vector> & reference,
    std::vector<Vector>  & derivatives,
    std::vector<Vector> & ddistdref,
    Tensor & Rotation,
    Matrix<std::vector<Vector> > &DRotDPos,
    bool squared) const {
  //initialize the data into the structure
  // typically the positions do not have the com neither calculated nor subtracted. This layer takes care of this business
  RMSDCoreData cd(align,displace,positions,reference);
  // transfer the settings for the center to let the CoreCalc deal with it
  // transfer the settings for the center to let the CoreCalc deal with it
  cd.setPositionsCenterIsRemoved(positions_center_is_removed);
  if(positions_center_is_calculated) {cd.setPositionsCenter(positions_center);}
  else {cd.calcPositionsCenter();};

  cd.setReferenceCenterIsRemoved(reference_center_is_removed);
  if(!reference_center_is_calculated) {cd.calcReferenceCenter();}
  else {cd.setReferenceCenter(reference_center);}

  // Perform the diagonalization and all the needed stuff
  cd.doCoreCalc(safe,alEqDis);
  // make the core calc distance
  double dist=cd.getDistance(squared);
//  make the derivatives by using pieces calculated in coreCalc (probably the best is just to copy the vector...)
  derivatives=cd.getDDistanceDPositions();
  ddistdref=cd.getDDistanceDReference();
  // get the rotation matrix
  Rotation=cd.getRotationMatrixReferenceToPositions();
  // get its derivative
  DRotDPos=cd.getDRotationDPositions();
  return dist;
}

template <bool safe,bool alEqDis>
double RMSD::optimalAlignment_DDistDRef_Rot_DRotDPos_DRotDRef(const  std::vector<double>  & align,
    const  std::vector<double>  & displace,
    const std::vector<Vector> & positions,
    const std::vector<Vector> & reference,
    std::vector<Vector>  & derivatives,
    std::vector<Vector> & ddistdref,
    Tensor & Rotation,
    Matrix<std::vector<Vector> > &DRotDPos,
    Matrix<std::vector<Vector> > &DRotDRef,
    bool squared) const {
  //initialize the data into the structure
  // typically the positions do not have the com neither calculated nor subtracted. This layer takes care of this business
  RMSDCoreData cd(align,displace,positions,reference);
  // transfer the settings for the center to let the CoreCalc deal with it
  // transfer the settings for the center to let the CoreCalc deal with it
  cd.setPositionsCenterIsRemoved(positions_center_is_removed);
  if(positions_center_is_calculated) {cd.setPositionsCenter(positions_center);}
  else {cd.calcPositionsCenter();};

  cd.setReferenceCenterIsRemoved(reference_center_is_removed);
  if(!reference_center_is_calculated) {cd.calcReferenceCenter();}
  else {cd.setReferenceCenter(reference_center);}

  // Perform the diagonalization and all the needed stuff
  cd.doCoreCalc(safe,alEqDis);
  // make the core calc distance
  double dist=cd.getDistance(squared);
//  make the derivatives by using pieces calculated in coreCalc (probably the best is just to copy the vector...)
  derivatives=cd.getDDistanceDPositions();
  ddistdref=cd.getDDistanceDReference();
  // get the rotation matrix
  Rotation=cd.getRotationMatrixReferenceToPositions();
  // get its derivative
  DRotDPos=cd.getDRotationDPositions();
  DRotDRef=cd.getDRotationDReference();
  return dist;
}

template <bool safe,bool alEqDis>
double RMSD::optimalAlignment_Rot_DRotDRr01(const  std::vector<double>  & align,
    const  std::vector<double>  & displace,
    const std::vector<Vector> & positions,
    const std::vector<Vector> & reference,
    Tensor & Rotation,
    std::array<std::array<Tensor,3>,3> & DRotDRr01,
    bool squared) const {
  //initialize the data into the structure
  // typically the positions do not have the com neither calculated nor subtracted. This layer takes care of this business
  RMSDCoreData cd(align,displace,positions,reference);
  // transfer the settings for the center to let the CoreCalc deal with it
  cd.setPositionsCenterIsRemoved(positions_center_is_removed);
  if(positions_center_is_calculated) {cd.setPositionsCenter(positions_center);}
  else {cd.calcPositionsCenter();};

  cd.setReferenceCenterIsRemoved(reference_center_is_removed);
  if(!reference_center_is_calculated) {cd.calcReferenceCenter();}
  else {cd.setReferenceCenter(reference_center);}

  // Perform the diagonalization and all the needed stuff
  cd.doCoreCalc(safe,alEqDis);
  // make the core calc distance
  double dist=cd.getDistance(squared);
  // get the rotation matrix
  Rotation=cd.getRotationMatrixReferenceToPositions();
  //get detivative w.r.t. rr01
  DRotDRr01=cd.getDRotationDRr01();
  return dist;
}

template <bool safe,bool alEqDis>
double RMSD::optimalAlignment_Rot(const  std::vector<double>  & align,
                                  const  std::vector<double>  & displace,
                                  const std::vector<Vector> & positions,
                                  const std::vector<Vector> & reference,
                                  std::vector<Vector>  & derivatives,
                                  Tensor & Rotation,
                                  bool squared) const {
  //initialize the data into the structure
  // typically the positions do not have the com neither calculated nor subtracted. This layer takes care of this business
  RMSDCoreData cd(align,displace,positions,reference);
  // transfer the settings for the center to let the CoreCalc deal with it
  cd.setPositionsCenterIsRemoved(positions_center_is_removed);
  if(positions_center_is_calculated) {cd.setPositionsCenter(positions_center);}
  else {cd.calcPositionsCenter();};

  cd.setReferenceCenterIsRemoved(reference_center_is_removed);
  if(!reference_center_is_calculated) {cd.calcReferenceCenter();}
  else {cd.setReferenceCenter(reference_center);}

  // Perform the diagonalization and all the needed stuff
  cd.doCoreCalc(safe,alEqDis);
  // make the core calc distance
  double dist=cd.getDistance(squared);
  //  make the derivatives by using pieces calculated in coreCalc (probably the best is just to copy the vector...)
  derivatives=cd.getDDistanceDPositions();
  // get the rotation matrix
  Rotation=cd.getRotationMatrixReferenceToPositions();
  return dist;
}

template <bool safe,bool alEqDis>
double RMSD::optimalAlignmentWithCloseStructure(const  std::vector<double>  & align,
    const  std::vector<double>  & displace,
    const std::vector<Vector> & positions,
    const std::vector<Vector> & reference,
    std::vector<Vector>  & derivatives,
    Tensor & rotationPosClose,
    Tensor & rotationRefClose,
    std::array<std::array<Tensor,3>,3> & drotationPosCloseDrr01,
    bool squared) const {
  //initialize the data into the structure
  // typically the positions do not have the com neither calculated nor subtracted. This layer takes care of this business
  RMSDCoreData cd(align,displace,positions,reference);
  // transfer the settings for the center to let the CoreCalc deal with it
  cd.setPositionsCenterIsRemoved(positions_center_is_removed);
  if(positions_center_is_calculated) {cd.setPositionsCenter(positions_center);}
  else {cd.calcPositionsCenter();};

  cd.setReferenceCenterIsRemoved(reference_center_is_removed);
  if(!reference_center_is_calculated) {cd.calcReferenceCenter();}
  else {cd.setReferenceCenter(reference_center);}

  // instead of diagonalization, approximate with saved rotation matrix
  cd.doCoreCalcWithCloseStructure(safe,alEqDis, rotationPosClose, rotationRefClose, drotationPosCloseDrr01);
  // make the core calc distance
  double dist=cd.getDistance(squared);
  //  make the derivatives by using pieces calculated in coreCalc (probably the best is just to copy the vector...)
  derivatives=cd.getDDistanceDPositions();
  return dist;
}


template <bool safe,bool alEqDis>
double RMSD::optimalAlignment_PCA(const  std::vector<double>  & align,
                                  const  std::vector<double>  & displace,
                                  const std::vector<Vector> & positions,
                                  const std::vector<Vector> & reference,
                                  std::vector<Vector> & alignedpositions,
                                  std::vector<Vector> & centeredpositions,
                                  std::vector<Vector> & centeredreference,
                                  Tensor & Rotation,
                                  std::vector<Vector> & DDistDPos,
                                  Matrix<std::vector<Vector> > & DRotDPos,
                                  bool squared) const {
  //initialize the data into the structure
  // typically the positions do not have the com neither calculated nor subtracted. This layer takes care of this business
  RMSDCoreData cd(align,displace,positions,reference);
  // transfer the settings for the center to let the CoreCalc deal with it
  cd.setPositionsCenterIsRemoved(positions_center_is_removed);
  if(positions_center_is_calculated) {cd.setPositionsCenter(positions_center);}
  else {cd.calcPositionsCenter();};

  cd.setReferenceCenterIsRemoved(reference_center_is_removed);
  if(!reference_center_is_calculated) {cd.calcReferenceCenter();}
  else {cd.setReferenceCenter(reference_center);}

  // Perform the diagonalization and all the needed stuff
  cd.doCoreCalc(safe,alEqDis);
  // make the core calc distance
  double dist=cd.getDistance(squared);
  // make the derivatives by using pieces calculated in coreCalc (probably the best is just to copy the vector...)
  DDistDPos=cd.getDDistanceDPositions();
  // get the rotation matrix
  Rotation=cd.getRotationMatrixPositionsToReference();
  // get its derivative
  DRotDPos=cd.getDRotationDPositions(true); // this gives back the inverse
  // get aligned positions
  alignedpositions=cd.getAlignedPositionsToReference();
  // get centered positions
  centeredpositions=cd.getCenteredPositions();
  // get centered reference
  centeredreference=cd.getCenteredReference();
  return dist;
}


template <bool safe,bool alEqDis>
double RMSD::optimalAlignment_Fit(const  std::vector<double>  & align,
                                  const  std::vector<double>  & displace,
                                  const std::vector<Vector> & positions,
                                  const std::vector<Vector> & reference,
                                  Tensor & Rotation,
                                  Matrix<std::vector<Vector> > & DRotDPos,
                                  std::vector<Vector> & centeredpositions,
                                  Vector & center_positions,
                                  bool squared) {
  //initialize the data into the structure
  // typically the positions do not have the com neither calculated nor subtracted. This layer takes care of this business
  RMSDCoreData cd(align,displace,positions,reference);
  // transfer the settings for the center to let the CoreCalc deal with it
  cd.setPositionsCenterIsRemoved(positions_center_is_removed);
  if(positions_center_is_calculated) {cd.setPositionsCenter(positions_center);}
  else {cd.calcPositionsCenter();};

  cd.setReferenceCenterIsRemoved(reference_center_is_removed);
  if(!reference_center_is_calculated) {cd.calcReferenceCenter();}
  else {cd.setReferenceCenter(reference_center);}

  // Perform the diagonalization and all the needed stuff
  cd.doCoreCalc(safe,alEqDis);
  // make the core calc distance
  double dist=cd.getDistance(squared);
  // get the rotation matrix
  Rotation=cd.getRotationMatrixPositionsToReference();
  // get its derivative
  DRotDPos=cd.getDRotationDPositions(true); // this gives back the inverse
  // get centered positions
  centeredpositions=cd.getCenteredPositions();
  // get center
  center_positions=cd.getPositionsCenter();
  return dist;
}






/// This calculates the elements needed by the quaternion to calculate everything that is needed
/// additional calls retrieve different components
/// note that this considers that the centers of both reference and positions are already setted
/// but automatically should properly account for non removed components: if not removed then it
/// removes prior to calculation of the alignment
void RMSDCoreData::doCoreCalc(bool safe,bool alEqDis, bool only_rotation) {

  retrieve_only_rotation=only_rotation;
  const unsigned n=static_cast<unsigned int>(reference.size());

  plumed_massert(creference_is_calculated,"the center of the reference frame must be already provided at this stage");
  plumed_massert(cpositions_is_calculated,"the center of the positions frame must be already provided at this stage");

// This is the trace of positions*positions + reference*reference
  rr00=0.;
  rr11=0.;
// This is positions*reference
  Tensor rr01;
// center of mass managing: must subtract the center from the position or not?
  Vector cp; cp.zero(); if(!cpositions_is_removed)cp=cpositions;
  Vector cr; cr.zero(); if(!creference_is_removed)cr=creference;
// second expensive loop: compute second moments wrt centers
  for(unsigned iat=0; iat<n; iat++) {
    double w=align[iat];
    rr00+=dotProduct(positions[iat]-cp,positions[iat]-cp)*w;
    rr11+=dotProduct(reference[iat]-cr,reference[iat]-cr)*w;
    rr01+=Tensor(positions[iat]-cp,reference[iat]-cr)*w;
  }

// the quaternion matrix: this is internal
  Tensor4d m;

  m[0][0]=2.0*(-rr01[0][0]-rr01[1][1]-rr01[2][2]);
  m[1][1]=2.0*(-rr01[0][0]+rr01[1][1]+rr01[2][2]);
  m[2][2]=2.0*(+rr01[0][0]-rr01[1][1]+rr01[2][2]);
  m[3][3]=2.0*(+rr01[0][0]+rr01[1][1]-rr01[2][2]);
  m[0][1]=2.0*(-rr01[1][2]+rr01[2][1]);
  m[0][2]=2.0*(+rr01[0][2]-rr01[2][0]);
  m[0][3]=2.0*(-rr01[0][1]+rr01[1][0]);
  m[1][2]=2.0*(-rr01[0][1]-rr01[1][0]);
  m[1][3]=2.0*(-rr01[0][2]-rr01[2][0]);
  m[2][3]=2.0*(-rr01[1][2]-rr01[2][1]);
  m[1][0] = m[0][1];
  m[2][0] = m[0][2];
  m[2][1] = m[1][2];
  m[3][0] = m[0][3];
  m[3][1] = m[1][3];
  m[3][2] = m[2][3];


  Tensor dm_drr01[4][4];
  if(!alEqDis or !retrieve_only_rotation) {
    dm_drr01[0][0] = 2.0*Tensor(-1.0, 0.0, 0.0,  0.0,-1.0, 0.0,  0.0, 0.0,-1.0);
    dm_drr01[1][1] = 2.0*Tensor(-1.0, 0.0, 0.0,  0.0,+1.0, 0.0,  0.0, 0.0,+1.0);
    dm_drr01[2][2] = 2.0*Tensor(+1.0, 0.0, 0.0,  0.0,-1.0, 0.0,  0.0, 0.0,+1.0);
    dm_drr01[3][3] = 2.0*Tensor(+1.0, 0.0, 0.0,  0.0,+1.0, 0.0,  0.0, 0.0,-1.0);
    dm_drr01[0][1] = 2.0*Tensor( 0.0, 0.0, 0.0,  0.0, 0.0,-1.0,  0.0,+1.0, 0.0);
    dm_drr01[0][2] = 2.0*Tensor( 0.0, 0.0,+1.0,  0.0, 0.0, 0.0, -1.0, 0.0, 0.0);
    dm_drr01[0][3] = 2.0*Tensor( 0.0,-1.0, 0.0, +1.0, 0.0, 0.0,  0.0, 0.0, 0.0);
    dm_drr01[1][2] = 2.0*Tensor( 0.0,-1.0, 0.0, -1.0, 0.0, 0.0,  0.0, 0.0, 0.0);
    dm_drr01[1][3] = 2.0*Tensor( 0.0, 0.0,-1.0,  0.0, 0.0, 0.0, -1.0, 0.0, 0.0);
    dm_drr01[2][3] = 2.0*Tensor( 0.0, 0.0, 0.0,  0.0, 0.0,-1.0,  0.0,-1.0, 0.0);
    dm_drr01[1][0] = dm_drr01[0][1];
    dm_drr01[2][0] = dm_drr01[0][2];
    dm_drr01[2][1] = dm_drr01[1][2];
    dm_drr01[3][0] = dm_drr01[0][3];
    dm_drr01[3][1] = dm_drr01[1][3];
    dm_drr01[3][2] = dm_drr01[2][3];
  }


  Vector4d q;

  Tensor dq_drr01[4];
  if(!alEqDis or !only_rotation) {
    diagMatSym(m, eigenvals, eigenvecs );
    q=Vector4d(eigenvecs[0][0],eigenvecs[0][1],eigenvecs[0][2],eigenvecs[0][3]);
    double dq_dm[4][4][4];
    for(unsigned i=0; i<4; i++) for(unsigned j=0; j<4; j++) for(unsigned k=0; k<4; k++) {
          double tmp=0.0;
// perturbation theory for matrix m
          for(unsigned l=1; l<4; l++) tmp+=eigenvecs[l][j]*eigenvecs[l][i]/(eigenvals[0]-eigenvals[l])*eigenvecs[0][k];
          dq_dm[i][j][k]=tmp;
        }
// propagation to _drr01
    for(unsigned i=0; i<4; i++) {
      Tensor tmp;
      for(unsigned j=0; j<4; j++) for(unsigned k=0; k<4; k++) {
          tmp+=dq_dm[i][j][k]*dm_drr01[j][k];
        }
      dq_drr01[i]=tmp;
    }
  } else {
    TensorGeneric<1,4> here_eigenvecs;
    VectorGeneric<1> here_eigenvals;
    diagMatSym(m, here_eigenvals, here_eigenvecs );
    for(unsigned i=0; i<4; i++) eigenvecs[0][i]=here_eigenvecs[0][i];
    eigenvals[0]=here_eigenvals[0];
    q=Vector4d(eigenvecs[0][0],eigenvecs[0][1],eigenvecs[0][2],eigenvecs[0][3]);
  }

// This is the rotation matrix that brings reference to positions
// i.e. matmul(rotation,reference[iat])+shift is fitted to positions[iat]

  rotation[0][0]=q[0]*q[0]+q[1]*q[1]-q[2]*q[2]-q[3]*q[3];
  rotation[1][1]=q[0]*q[0]-q[1]*q[1]+q[2]*q[2]-q[3]*q[3];
  rotation[2][2]=q[0]*q[0]-q[1]*q[1]-q[2]*q[2]+q[3]*q[3];
  rotation[0][1]=2*(+q[0]*q[3]+q[1]*q[2]);
  rotation[0][2]=2*(-q[0]*q[2]+q[1]*q[3]);
  rotation[1][2]=2*(+q[0]*q[1]+q[2]*q[3]);
  rotation[1][0]=2*(-q[0]*q[3]+q[1]*q[2]);
  rotation[2][0]=2*(+q[0]*q[2]+q[1]*q[3]);
  rotation[2][1]=2*(-q[0]*q[1]+q[2]*q[3]);


  if(!alEqDis or !only_rotation) {
    drotation_drr01[0][0]=2*q[0]*dq_drr01[0]+2*q[1]*dq_drr01[1]-2*q[2]*dq_drr01[2]-2*q[3]*dq_drr01[3];
    drotation_drr01[1][1]=2*q[0]*dq_drr01[0]-2*q[1]*dq_drr01[1]+2*q[2]*dq_drr01[2]-2*q[3]*dq_drr01[3];
    drotation_drr01[2][2]=2*q[0]*dq_drr01[0]-2*q[1]*dq_drr01[1]-2*q[2]*dq_drr01[2]+2*q[3]*dq_drr01[3];
    drotation_drr01[0][1]=2*(+(q[0]*dq_drr01[3]+dq_drr01[0]*q[3])+(q[1]*dq_drr01[2]+dq_drr01[1]*q[2]));
    drotation_drr01[0][2]=2*(-(q[0]*dq_drr01[2]+dq_drr01[0]*q[2])+(q[1]*dq_drr01[3]+dq_drr01[1]*q[3]));
    drotation_drr01[1][2]=2*(+(q[0]*dq_drr01[1]+dq_drr01[0]*q[1])+(q[2]*dq_drr01[3]+dq_drr01[2]*q[3]));
    drotation_drr01[1][0]=2*(-(q[0]*dq_drr01[3]+dq_drr01[0]*q[3])+(q[1]*dq_drr01[2]+dq_drr01[1]*q[2]));
    drotation_drr01[2][0]=2*(+(q[0]*dq_drr01[2]+dq_drr01[0]*q[2])+(q[1]*dq_drr01[3]+dq_drr01[1]*q[3]));
    drotation_drr01[2][1]=2*(-(q[0]*dq_drr01[1]+dq_drr01[0]*q[1])+(q[2]*dq_drr01[3]+dq_drr01[2]*q[3]));
  }

  d.resize(n);

  // calculate rotation matrix derivatives and components distances needed for components only when align!=displacement
  if(!alEqDis)ddist_drotation.zero();
  #pragma omp simd
  for(unsigned iat=0; iat<n; iat++) {
    // components differences: this is useful externally
    d[iat]=positions[iat]-cp - matmul(rotation,reference[iat]-cr);
    //cerr<<"D "<<iat<<" "<<d[iat][0]<<" "<<d[iat][1]<<" "<<d[iat][2]<<"\n";
  }
  // ddist_drotation if needed
  if(!alEqDis or !only_rotation)
    for (unsigned iat=0; iat<n; iat++)
      ddist_drotation+=-2*displace[iat]*extProduct(d[iat],reference[iat]-cr);

  if(!alEqDis or !only_rotation) {
    ddist_drr01.zero();
    for(unsigned i=0; i<3; i++) for(unsigned j=0; j<3; j++) ddist_drr01+=ddist_drotation[i][j]*drotation_drr01[i][j];
  }
  // transfer this bools to the cd so that this settings will be reflected in the other calls
  this->alEqDis=alEqDis;
  this->safe=safe;
  isInitialized=true;

}
/// just retrieve the distance already calculated
double RMSDCoreData::getDistance( bool squared) {

  if(!isInitialized)plumed_merror("getDistance cannot calculate the distance without being initialized first by doCoreCalc ");

  double localDist=0.0;
  const unsigned n=static_cast<unsigned int>(reference.size());
  if(safe || !alEqDis) localDist=0.0;
  else
    localDist=eigenvals[0]+rr00+rr11;
  #pragma omp simd reduction(+:localDist)
  for(unsigned iat=0; iat<n; iat++) {
    if(alEqDis) {
      if(safe) localDist+=align[iat]*modulo2(d[iat]);
    } else {
      localDist+=displace[iat]*modulo2(d[iat]);
    }
  }
  if(!squared) {
    dist=std::sqrt(localDist);
    distanceIsMSD=false;
  } else {
    dist=localDist;
    distanceIsMSD=true;
  }
  hasDistance=true;
  return dist;
}

void RMSDCoreData::doCoreCalcWithCloseStructure(bool safe,bool alEqDis, Tensor & rotationPosClose, Tensor & rotationRefClose, std::array<std::array<Tensor,3>,3> & drotationPosCloseDrr01) {

  unsigned natoms = reference.size();
  Tensor ddist_drxy;
  ddist_drr01.zero();
  d.resize(natoms);

  // center of mass managing: must subtract the center from the position or not?
  Vector cp; cp.zero(); if(!cpositions_is_removed)cp=cpositions;
  Vector cr; cr.zero(); if(!creference_is_removed)cr=creference;
  //distance = \sum_{n=0}^{N} w_n(x_n-cpos-R_{XY} R_{AY} a_n)^2

  Tensor rotation = matmul(rotationPosClose, rotationRefClose);

  #pragma omp simd
  for (unsigned iat=0; iat<natoms; iat++) {
    d[iat] = positions[iat] - cp - matmul(rotation, reference[iat]-cr);
  }
  if (!alEqDis) {
    for (unsigned iat=0; iat<natoms; iat++) {
      //dist = \sum w_i(x_i - cpos - R_xy * R_ay * a_i)
      ddist_drxy += -2*displace[iat]*extProduct(matmul(d[iat], rotationRefClose), reference[iat]-cr);
    }
  }

  if (!alEqDis) {
    for(unsigned i=0; i<3; i++)
      for(unsigned j=0; j<3; j++)
        ddist_drr01+=ddist_drxy[i][j]*drotationPosCloseDrr01[i][j];
  }
  this->alEqDis=alEqDis;
  this->safe=safe;
  isInitialized=true;
}

std::vector<Vector> RMSDCoreData::getDDistanceDPositions() {
  std::vector<Vector>  derivatives;
  const unsigned n=static_cast<unsigned int>(reference.size());
  Vector ddist_dcpositions;
  derivatives.resize(n);
  double prefactor=1.0;
  if(!distanceIsMSD) prefactor*=0.5/dist;
  plumed_massert(!retrieve_only_rotation,"You used  only_rotation=true in doCoreCalc therefore you cannot retrieve this information now");
  if(!hasDistance)plumed_merror("getDPositionsDerivatives needs to calculate the distance via getDistance first !");
  if(!isInitialized)plumed_merror("getDPositionsDerivatives needs to initialize the coreData first!");
  Vector csum;
  for(unsigned iat=0; iat<n; iat++) {
    if(alEqDis) {
// there is no need for derivatives of rotation and shift here as it is by construction zero
// (similar to Hellman-Feynman forces)
      derivatives[iat]= 2*prefactor*align[iat]*d[iat];
    } else {
// these are the derivatives assuming the roto-translation as frozen
      Vector tmp1=2*displace[iat]*d[iat];
      derivatives[iat]=tmp1;
// derivative of cpositions
      ddist_dcpositions+=-tmp1;
      // these needed for com corrections
      Vector tmp2=matmul(ddist_drr01,reference[iat]-creference)*align[iat];
      derivatives[iat]+=tmp2;
      csum+=tmp2;
    }
  }

  if(!alEqDis)
    #pragma omp simd
    for(unsigned iat=0; iat<n; iat++) {derivatives[iat]= prefactor*(derivatives[iat]+(ddist_dcpositions-csum)*align[iat]); }

  return derivatives;
}

std::vector<Vector>  RMSDCoreData::getDDistanceDReference() {
  std::vector<Vector>  derivatives;
  const unsigned n=static_cast<unsigned int>(reference.size());
  Vector ddist_dcreference;
  derivatives.resize(n);
  double prefactor=1.0;
  if(!distanceIsMSD) prefactor*=0.5/dist;
  Vector csum;

  plumed_massert(!retrieve_only_rotation,"You used  only_rotation=true in doCoreCalc therefore you cannot retrieve this information now");
  if(!hasDistance)plumed_merror("getDDistanceDReference needs to calculate the distance via getDistance first !");
  if(!isInitialized)plumed_merror("getDDistanceDReference to initialize the coreData first!");
  // get the transpose rotation
  Tensor t_rotation=rotation.transpose();
  Tensor t_ddist_drr01=ddist_drr01.transpose();

// third expensive loop: derivatives
  for(unsigned iat=0; iat<n; iat++) {
    if(alEqDis) {
// there is no need for derivatives of rotation and shift here as it is by construction zero
// (similar to Hellman-Feynman forces)
      //TODO: check this derivative down here
      derivatives[iat]= -2*prefactor*align[iat]*matmul(t_rotation,d[iat]);
    } else {
// these are the derivatives assuming the roto-translation as frozen
      Vector tmp1=2*displace[iat]*matmul(t_rotation,d[iat]);
      derivatives[iat]= -tmp1;
// derivative of cpositions
      ddist_dcreference+=tmp1;
      // these below are needed for com correction
      Vector tmp2=matmul(t_ddist_drr01,positions[iat]-cpositions)*align[iat];
      derivatives[iat]+=tmp2;
      csum+=tmp2;
    }
  }

  if(!alEqDis)
    #pragma omp simd
    for(unsigned iat=0; iat<n; iat++) {derivatives[iat]= prefactor*(derivatives[iat]+(ddist_dcreference-csum)*align[iat]);}

  return derivatives;
}

/// this version does not calculate the derivative of rotation matrix as needed for SOMA
std::vector<Vector>  RMSDCoreData::getDDistanceDReferenceSOMA() {
  std::vector<Vector>  derivatives;
  const unsigned n=static_cast<unsigned int>(reference.size());
  Vector ddist_dcreference;
  derivatives.resize(n);
  double prefactor=1.0;
  if(!distanceIsMSD) prefactor*=0.5/dist;
  Vector csum,tmp1,tmp2;

  plumed_massert(!retrieve_only_rotation,"You used  only_rotation=true in doCoreCalc therefore you cannot retrieve this information now");
  if(!hasDistance)plumed_merror("getDDistanceDReference needs to calculate the distance via getDistance first !");
  if(!isInitialized)plumed_merror("getDDistanceDReference to initialize the coreData first!");
  // get the transpose rotation
  Tensor t_rotation=rotation.transpose();

// third expensive loop: derivatives
  for(unsigned iat=0; iat<n; iat++) {
    if(alEqDis) {
// there is no need for derivatives of rotation and shift here as it is by construction zero
// (similar to Hellman-Feynman forces)
      //TODO: check this derivative down here
      derivatives[iat]= -2*prefactor*align[iat]*matmul(t_rotation,d[iat]);
    } else {
// these are the derivatives assuming the roto-translation as frozen
      tmp1=2*displace[iat]*matmul(t_rotation,d[iat]);
      derivatives[iat]= -tmp1;
// derivative of cpositions
      ddist_dcreference+=tmp1;
    }
  }

  if(!alEqDis) for(unsigned iat=0; iat<n; iat++)derivatives[iat]=prefactor*(derivatives[iat]+ddist_dcreference*align[iat]);

  return derivatives;
}



/*
This below is the derivative of the rotation matrix that aligns the reference onto the positions
respect to positions
note that the this transformation overlap the  reference onto position
if inverseTransform=true then aligns the positions onto reference
*/
Matrix<std::vector<Vector> >  RMSDCoreData::getDRotationDPositions( bool inverseTransform ) {
  const unsigned n=static_cast<unsigned int>(reference.size());
  plumed_massert(!retrieve_only_rotation,"You used  only_rotation=true in doCoreCalc therefore you cannot retrieve this information now");
  if(!isInitialized)plumed_merror("getDRotationDPosition to initialize the coreData first!");
  Matrix<std::vector<Vector> > DRotDPos=Matrix<std::vector<Vector> >(3,3);
  // remember drotation_drr01 is Tensor drotation_drr01[3][3]
  //           (3x3 rot) (3x3 components of rr01)
  std::vector<Vector> v(n);
  Vector csum;
  // these below could probably be calculated in the main routine
  Vector cp; cp.zero(); if(!cpositions_is_removed)cp=cpositions;
  Vector cr; cr.zero(); if(!creference_is_removed)cr=creference;
  for(unsigned iat=0; iat<n; iat++) csum+=(reference[iat]-cr)*align[iat];
  for(unsigned iat=0; iat<n; iat++) v[iat]=(reference[iat]-cr-csum)*align[iat];
  for(unsigned a=0; a<3; a++) {
    for(unsigned b=0; b<3; b++) {
      if(inverseTransform) {
        DRotDPos[b][a].resize(n);
        for(unsigned iat=0; iat<n; iat++) {
          DRotDPos[b][a][iat]=matmul(drotation_drr01[a][b],v[iat]);
        }
      } else {
        DRotDPos[a][b].resize(n);
        for(unsigned iat=0; iat<n; iat++) {
          DRotDPos[a][b][iat]=matmul(drotation_drr01[a][b],v[iat]);
        }
      }
    }
  }
  return DRotDPos;
}

/*
This below is the derivative of the rotation matrix that aligns the reference onto the positions
respect to reference
note that the this transformation overlap the  reference onto position
if inverseTransform=true then aligns the positions onto reference
*/
Matrix<std::vector<Vector> >  RMSDCoreData::getDRotationDReference( bool inverseTransform ) {
  const unsigned n=static_cast<unsigned int>(reference.size());
  plumed_massert(!retrieve_only_rotation,"You used  only_rotation=true in doCoreCalc therefore you cannot retrieve this information now");
  if(!isInitialized)plumed_merror("getDRotationDPositions to initialize the coreData first!");
  Matrix<std::vector<Vector> > DRotDRef=Matrix<std::vector<Vector> >(3,3);
  // remember drotation_drr01 is Tensor drotation_drr01[3][3]
  //           (3x3 rot) (3x3 components of rr01)
  std::vector<Vector> v(n);
  Vector csum;
  // these below could probably be calculated in the main routine
  Vector cp; cp.zero(); if(!cpositions_is_removed)cp=cpositions;
  Vector cr; cr.zero(); if(!creference_is_removed)cr=creference;
  for(unsigned iat=0; iat<n; iat++) csum+=(positions[iat]-cp)*align[iat];
  for(unsigned iat=0; iat<n; iat++) v[iat]=(positions[iat]-cp-csum)*align[iat];

  for(unsigned a=0; a<3; a++) {
    for(unsigned b=0; b<3; b++) {
      Tensor t_drotation_drr01=drotation_drr01[a][b].transpose();
      if(inverseTransform) {
        DRotDRef[b][a].resize(n);
        for(unsigned iat=0; iat<n; iat++) {
          DRotDRef[b][a][iat]=matmul(t_drotation_drr01,v[iat]);
        }
      } else {
        DRotDRef[a][b].resize(n);
        for(unsigned iat=0; iat<n; iat++) {
          DRotDRef[a][b][iat]=matmul(t_drotation_drr01,v[iat]);
        }
      }
    }
  }
  return DRotDRef;
}


std::vector<Vector> RMSDCoreData::getAlignedReferenceToPositions() {
  std::vector<Vector> alignedref;
  const unsigned n=static_cast<unsigned int>(reference.size());
  alignedref.resize(n);
  if(!isInitialized)plumed_merror("getAlignedReferenceToPostions needs to initialize the coreData first!");
  // avoid to calculate matrix element but use the sum of what you have
  Vector cp; cp.zero(); if(!cpositions_is_removed)cp=cpositions;
  for(unsigned iat=0; iat<n; iat++)alignedref[iat]=-d[iat]+positions[iat]-cp;
  return alignedref;
}
std::vector<Vector> RMSDCoreData::getAlignedPositionsToReference() {
  std::vector<Vector> alignedpos;
  if(!isInitialized)plumed_merror("getAlignedPostionsToReference needs to initialize the coreData first!");
  const unsigned n=static_cast<unsigned int>(positions.size());
  alignedpos.resize(n);
  Vector cp; cp.zero(); if(!cpositions_is_removed)cp=cpositions;
  // avoid to calculate matrix element but use the sum of what you have
  for(unsigned iat=0; iat<n; iat++)alignedpos[iat]=matmul(rotation.transpose(),positions[iat]-cp);
  return alignedpos;
}


std::vector<Vector> RMSDCoreData::getCenteredPositions() {
  std::vector<Vector> centeredpos;
  const unsigned n=static_cast<unsigned int>(reference.size());
  centeredpos.resize(n);
  if(!isInitialized)plumed_merror("getCenteredPositions needs to initialize the coreData first!");
  // avoid to calculate matrix element but use the sum of what you have
  for(unsigned iat=0; iat<n; iat++)centeredpos[iat]=positions[iat]-cpositions;
  return centeredpos;
}

std::vector<Vector> RMSDCoreData::getCenteredReference() {
  std::vector<Vector> centeredref;
  const unsigned n=static_cast<unsigned int>(reference.size());
  centeredref.resize(n);
  if(!isInitialized)plumed_merror("getCenteredReference needs to initialize the coreData first!");
  // avoid to calculate matrix element but use the sum of what you have
  Vector cr; cr.zero(); if(!creference_is_removed)cr=creference;
  for(unsigned iat=0; iat<n; iat++)centeredref[iat]=reference[iat]-cr;
  return centeredref;
}


Vector RMSDCoreData::getPositionsCenter() {
  if(!isInitialized)plumed_merror("getCenteredPositions needs to initialize the coreData first!");
  return cpositions;
}

Vector RMSDCoreData::getReferenceCenter() {
  if(!isInitialized)plumed_merror("getCenteredPositions needs to initialize the coreData first!");
  return creference;
}

Tensor RMSDCoreData::getRotationMatrixReferenceToPositions() {
  if(!isInitialized)plumed_merror("getRotationMatrixReferenceToPositions needs to initialize the coreData first!");
  return rotation;
}

Tensor RMSDCoreData::getRotationMatrixPositionsToReference() {
  if(!isInitialized)plumed_merror("getRotationMatrixReferenceToPositions needs to initialize the coreData first!");
  return rotation.transpose();
}

const std::array<std::array<Tensor,3>,3> &  RMSDCoreData::getDRotationDRr01() const {
  if(!isInitialized)plumed_merror("getDRotationDRr01 needs to initialize the coreData first!");
  return drotation_drr01;
}



template double RMSD::optimalAlignment<true,true>(const  std::vector<double>  & align,
    const  std::vector<double>  & displace,
    const std::vector<Vector> & positions,
    const std::vector<Vector> & reference,
    std::vector<Vector>  & derivatives, bool squared)const;
template double RMSD::optimalAlignment<true,false>(const  std::vector<double>  & align,
    const  std::vector<double>  & displace,
    const std::vector<Vector> & positions,
    const std::vector<Vector> & reference,
    std::vector<Vector>  & derivatives, bool squared)const;
template double RMSD::optimalAlignment<false,true>(const  std::vector<double>  & align,
    const  std::vector<double>  & displace,
    const std::vector<Vector> & positions,
    const std::vector<Vector> & reference,
    std::vector<Vector>  & derivatives, bool squared)const;
template double RMSD::optimalAlignment<false,false>(const  std::vector<double>  & align,
    const  std::vector<double>  & displace,
    const std::vector<Vector> & positions,
    const std::vector<Vector> & reference,
    std::vector<Vector>  & derivatives, bool squared)const;



}
