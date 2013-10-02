/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013 The plumed team
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
#include "RMSD.h"
#include "PDB.h"
#include "Log.h"
#include "Exception.h"
#include <cmath>
#include <iostream>
#include "Matrix.h"
#include "Tools.h"
#include "Tensor.h"

using namespace std;
namespace PLMD{

RMSD::RMSD(Log & log ):
  alignmentMethod(SIMPLE),
  log(&log){}

RMSD& RMSD::operator=(const RMSD& v){
  alignmentMethod=v.alignmentMethod;
  reference=v.reference;
  align=v.align;
  displace=v.displace;

  log=v.log;
  return *this ;
}


RMSD::RMSD(const RMSD & oldrmsd):
  alignmentMethod(oldrmsd.alignmentMethod),
  reference(oldrmsd.reference),
  align(oldrmsd.align),
  displace(oldrmsd.align),
  log( oldrmsd.log )
  {  }

void RMSD::set(const PDB&pdb, string mytype ){

	setReference(pdb.getPositions());
	setAlign(pdb.getOccupancy());
	setDisplace(pdb.getBeta());
        setType(mytype);
}

void RMSD::setType(string mytype){

	alignmentMethod=SIMPLE; // initialize with the simplest case: no rotation
	if (mytype=="SIMPLE"){
		alignmentMethod=SIMPLE;
		log->printf("RMSD IS DONE WITH SIMPLE METHOD(NO ROTATION)\n");
	}
	else if (mytype=="OPTIMAL"){
		alignmentMethod=OPTIMAL;
		log->printf("RMSD IS DONE WITH OPTIMAL ALIGNMENT METHOD\n");
	}
	else if (mytype=="OPTIMAL-FAST"){
		alignmentMethod=OPTIMAL_FAST;
		log->printf("RMSD IS DONE WITH OPTIMAL-FAST ALIGNMENT METHOD (fast version, numerically less stable, only valid with align==displace)\n");
	}
	else plumed_merror("unknown RMSD type" + mytype);

}

void RMSD::clear(){
  reference.clear();
  align.clear();
  displace.clear();
}

RMSD::~RMSD(){
}

string RMSD::getMethod(){
	string mystring;
	switch(alignmentMethod){
		case SIMPLE: mystring.assign("SIMPLE");break; 
		case OPTIMAL: mystring.assign("OPTIMAL");break; 
		case OPTIMAL_FAST: mystring.assign("OPTIMAL-FAST");break; 
	}	
	return mystring;
}

void RMSD::setReference(const vector<Vector> & reference){
  unsigned n=reference.size();
  this->reference=reference;
  plumed_massert(align.empty(),"you should first clear() an RMSD object, then set a new referece");
  plumed_massert(displace.empty(),"you should first clear() an RMSD object, then set a new referece");
  align.resize(n,1.0);
  displace.resize(n,1.0);
}

void RMSD::setAlign(const vector<double> & align){
  plumed_massert(this->align.size()==align.size(),"mismatch in dimension of align/displace arrays");
  this->align=align;
}

void RMSD::setDisplace(const vector<double> & displace){
  plumed_massert(this->displace.size()==displace.size(),"mismatch in dimension of align/displace arrays");
  this->displace=displace;
}

double RMSD::calculate(const std::vector<Vector> & positions,std::vector<Vector> &derivatives, bool squared){

  double ret=0.;

  switch(alignmentMethod){
	case SIMPLE:
		//	do a simple alignment without rotation 
		ret=simpleAlignment(align,displace,positions,reference,log,derivatives,squared);
		break;	
	case OPTIMAL_FAST:
		// this is calling the fastest option:
                if(align==displace) ret=optimalAlignment<false,true>(align,displace,positions,reference,derivatives,squared); 
                else                ret=optimalAlignment<false,false>(align,displace,positions,reference,derivatives,squared); 
		break;
	case OPTIMAL:
		// this is the fast routine but in the "safe" mode, which gives less numerical error:
		if(align==displace) ret=optimalAlignment<true,true>(align,displace,positions,reference,derivatives,squared); 
		else ret=optimalAlignment<true,false>(align,displace,positions,reference,derivatives,squared); 
		break;	
  }	

  return ret;

}

double RMSD::simpleAlignment(const  std::vector<double>  & align,
		                     const  std::vector<double>  & displace,
		                     const std::vector<Vector> & positions,
		                     const std::vector<Vector> & reference ,
		                     Log* &log,
		                     std::vector<Vector>  & derivatives, bool squared) {
      double dist(0);
      double anorm(0);
      double dnorm(0);
      unsigned n=reference.size();

      Vector apositions;
      Vector areference;
      Vector dpositions;
      Vector dreference;

      for(unsigned i=0;i<n;i++){
        double aw=align[i];
        double dw=displace[i];
        anorm+=aw;
        dnorm+=dw;
        apositions+=positions[i]*aw;
        areference+=reference[i]*aw;
        dpositions+=positions[i]*dw;
        dreference+=reference[i]*dw;
      }

      double invdnorm=1.0/dnorm;
      double invanorm=1.0/anorm;

      apositions*=invanorm;
      areference*=invanorm;
      dpositions*=invdnorm;
      dreference*=invdnorm;

      Vector shift=((apositions-areference)-(dpositions-dreference));
      for(unsigned i=0;i<n;i++){
        Vector d=(positions[i]-apositions)-(reference[i]-areference);
        dist+=displace[i]*d.modulo2();
        derivatives[i]=2*(invdnorm*displace[i]*d+invanorm*align[i]*shift);
      }
     dist*=invdnorm;

     if(!squared){
	// sqrt and normalization
        dist=sqrt(dist);
	///// sqrt and normalization on derivatives
        for(unsigned i=0;i<n;i++){derivatives[i]*=(0.5/dist);}
      }
      return dist;
}

// notice that in the current implementation the safe argument only makes sense for
// align==displace
template <bool safe,bool alEqDis>
double RMSD::optimalAlignment(const  std::vector<double>  & align,
                                     const  std::vector<double>  & displace,
                                     const std::vector<Vector> & positions,
                                     const std::vector<Vector> & reference ,
                                     std::vector<Vector>  & derivatives, bool squared) {
  double dist(0);
  double norm(0);
  double dnorm(0);
  const unsigned n=reference.size();
// This is the trace of positions*positions + reference*reference
  double sum00w(0);
  double sum11w(0);
// This is positions*reference
  Tensor sum01w;

  derivatives.resize(n);

  Vector cpositions;
  Vector creference;

// first expensive loop: compute centers
  for(unsigned iat=0;iat<n;iat++){
    double w=align[iat];
    norm+=w;
    if(!alEqDis) dnorm+=displace[iat];
    cpositions+=positions[iat]*w;
    creference+=reference[iat]*w;
  }
  double invnorm=1.0/norm;
  double invdnorm;
  if(!alEqDis) invdnorm=1.0/dnorm;

  cpositions*=invnorm;
  creference*=invnorm;
  
// second expensive loop: compute second moments wrt centers
  for(unsigned iat=0;iat<n;iat++){
    double w=align[iat];
    sum00w+=dotProduct(positions[iat]-cpositions,positions[iat]-cpositions)*w;
    sum11w+=dotProduct(reference[iat]-creference,reference[iat]-creference)*w;
    sum01w+=Tensor(positions[iat]-cpositions,reference[iat]-creference)*w;
  }

  double rr00=sum00w*invnorm;
  Tensor rr01=sum01w*invnorm;
  double rr11=sum11w*invnorm;

  Matrix<double> m=Matrix<double>(4,4);
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
  if(!alEqDis){
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

  vector<double> eigenvals;
  Matrix<double> eigenvecs;
  int diagerror=diagMat(m, eigenvals, eigenvecs );

  if (diagerror!=0){
    string sdiagerror;
    Tools::convert(diagerror,sdiagerror);
    string msg="DIAGONALIZATION FAILED WITH ERROR CODE "+sdiagerror;
    plumed_merror(msg);
  }

  dist=eigenvals[0]+rr00+rr11;

  Matrix<double> ddist_dm(4,4);

  Vector4d q(eigenvecs[0][0],eigenvecs[0][1],eigenvecs[0][2],eigenvecs[0][3]);

  Tensor dq_drr01[4];
  if(!alEqDis){
    double dq_dm[4][4][4];
    for(unsigned i=0;i<4;i++) for(unsigned j=0;j<4;j++) for(unsigned k=0;k<4;k++){
      double tmp=0.0;
// perturbation theory for matrix m
      for(unsigned l=1;l<4;l++) tmp+=eigenvecs[l][j]*eigenvecs[l][i]/(eigenvals[0]-eigenvals[l])*eigenvecs[0][k];
      dq_dm[i][j][k]=tmp;
    }
// propagation to _drr01
    for(unsigned i=0;i<4;i++){
      Tensor tmp;
      for(unsigned j=0;j<4;j++) for(unsigned k=0;k<4;k++) {
        tmp+=dq_dm[i][j][k]*dm_drr01[j][k];
      }
      dq_drr01[i]=tmp;
    }
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

  
  Tensor drotation_drr01[3][3];
  if(!alEqDis){
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

  double prefactor=2.0*invnorm;

  if(!squared && alEqDis) prefactor*=0.5/sqrt(dist);

// if "safe", recompute dist here to a better accuracy
  if(safe || !alEqDis) dist=0.0;

// If safe is set to "false", MSD is taken from the eigenvalue of the M matrix
// If safe is set to "true", MSD is recomputed from the rotational matrix
// For some reason, this last approach leads to less numerical noise but adds an overhead

  Tensor ddist_drotation;
  Vector ddist_dcpositions;

// third expensive loop: derivatives
  for(unsigned iat=0;iat<n;iat++){
    Vector d(positions[iat]-cpositions - matmul(rotation,reference[iat]-creference));
    if(alEqDis){
// there is no need for derivatives of rotation and shift here as it is by construction zero
// (similar to Hellman-Feynman forces)
      derivatives[iat]= prefactor*align[iat]*d;
       if(safe) dist+=align[iat]*invnorm*modulo2(d);
    } else {
// the case for align != displace is different, sob:
      dist+=displace[iat]*invdnorm*modulo2(d);
// these are the derivatives assuming the roto-translation as frozen
      derivatives[iat]=2*displace[iat]*invdnorm*d;
// here I accumulate derivatives wrt rotation matrix ..
      ddist_drotation+=-2*displace[iat]*invdnorm*extProduct(d,reference[iat]-creference);
// .. and cpositions
      ddist_dcpositions+=-2*displace[iat]*invdnorm*d;
    }
  }

  if(!alEqDis){
    Tensor ddist_drr01;
    for(unsigned i=0;i<3;i++) for(unsigned j=0;j<3;j++) ddist_drr01+=ddist_drotation[i][j]*drotation_drr01[i][j];
    for(unsigned iat=0;iat<n;iat++){
// this is propagating to positions.
// I am implicitly using the derivative of rr01 wrt positions here
      derivatives[iat]+=matmul(ddist_drr01,(reference[iat]-creference))*align[iat]*invnorm;
      derivatives[iat]+=ddist_dcpositions*align[iat]*invnorm;
    }
  }
  if(!squared){
    dist=sqrt(dist);
    if(!alEqDis){
      double xx=0.5/dist;
      for(unsigned iat=0;iat<n;iat++) derivatives[iat]*=xx;
    }
  }

  return dist;
}

}
