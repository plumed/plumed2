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
#include "RMSD.h"
#include "PDB.h"
#include "Log.h"
#include "OptimalAlignment.h"
#include "Exception.h"
#include <cmath>
#include <iostream>
#include "Matrix.h"
#include "Tools.h"

using namespace std;
namespace PLMD{

RMSD::RMSD(Log & log ):
  alignmentMethod(SIMPLE),
  myoptimalalignment(NULL),
  log(&log){}

RMSD& RMSD::operator=(const RMSD& v){
  alignmentMethod=v.alignmentMethod;
  reference=v.reference;
  align=v.align;
  displace=v.displace;
// in this manner the new RMSD is built empty and will just allocate its own
// myoptimalalignment when used (in calculate())
  myoptimalalignment=NULL;

  log=v.log;
  return *this ;
}


RMSD::RMSD(const RMSD & oldrmsd):
  alignmentMethod(oldrmsd.alignmentMethod),
  reference(oldrmsd.reference),
  align(oldrmsd.align),
  displace(oldrmsd.align),
// in this manner the new RMSD is built empty and will just allocate its own
// myoptimalalignment when used (in calculate())
  myoptimalalignment(NULL),
  log( oldrmsd.log )
  {  }

void RMSD::set(const PDB&pdb, string mytype ){

	setReference(pdb.getPositions());
	setAlign(pdb.getOccupancy());
	setDisplace(pdb.getBeta());
        setType(mytype);
}

void RMSD::setType(string mytype){
	myoptimalalignment=NULL;

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
	if(myoptimalalignment!=NULL) delete myoptimalalignment;
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
		ret=optimalAlignment(align,displace,positions,reference,derivatives,squared); 
		break;
	case OPTIMAL:
		if (myoptimalalignment==NULL){ // do full initialization	
		//
		// I create the object only here
		// since the alignment object require to know both position and reference
		// and it is possible only at calculate time
		//
			myoptimalalignment=new OptimalAlignment(align,displace,positions,reference,log);
        	}
		// this changes the P0 according the running frame
		(*myoptimalalignment).assignP0(positions);
		ret=(*myoptimalalignment).calculate(squared, derivatives);
		//(*myoptimalalignment).weightedFindiffTest(false);
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
	  double norm(0);
	  unsigned n=reference.size();
	  for(unsigned i=0;i<n;i++){
	      Vector d=delta(reference[i],positions[i]);
	      derivatives[i]=2.0*displace[i]*d;
	      dist+=displace[i]*d.modulo2();
	      norm+=displace[i];
      }

     double ret; 
     if(!squared){
	// sqrt and normalization
        ret=sqrt(dist/norm);
	///// sqrt and normalization on derivatives
        for(unsigned i=0;i<n;i++){derivatives[i]*=(0.5/ret/norm);}
      }else{
        ret=dist/norm;
      }
      return ret;
}

double RMSD::optimalAlignment(const  std::vector<double>  & align,
                                     const  std::vector<double>  & displace,
                                     const std::vector<Vector> & positions,
                                     const std::vector<Vector> & reference ,
                                     std::vector<Vector>  & derivatives, bool squared) {
  plumed_massert(displace==align,"OPTIMAL_FAST version of RMSD can only be used when displace weights are same as align weights");

  double dist(0);
  double norm(0);
  const unsigned n=reference.size();
// for these tensors I directly accumulate the trace:
  double sum00w(0);
  double sum11w(0);
// for these
  Tensor sum01w;

  derivatives.resize(n);

  Vector cpositions;
  Vector creference;

// first expensive loop: compute centers
  for(unsigned iat=0;iat<n;iat++){
    double w=align[iat];
    norm+=w;
    cpositions+=positions[iat]*w;
    creference+=reference[iat]*w;
  }
  cpositions/=norm;
  creference/=norm;
  
// second expensive loop: compute second moments wrt centers
  for(unsigned iat=0;iat<n;iat++){
    double w=align[iat];
    sum00w+=dotProduct(positions[iat]-cpositions,positions[iat]-cpositions)*w;
    sum01w+=Tensor(positions[iat]-cpositions,reference[iat]-creference)*w;
    sum11w+=dotProduct(reference[iat]-creference,reference[iat]-creference)*w;
  }

  double rr00=sum00w/norm;
  Tensor rr01=sum01w/norm;
  double rr11=sum11w/norm;

  Matrix<double> m=Matrix<double>(4,4);
  m[0][0]=rr00+rr11+2.0*(-rr01[0][0]-rr01[1][1]-rr01[2][2]);
  m[1][1]=rr00+rr11+2.0*(-rr01[0][0]+rr01[1][1]+rr01[2][2]);
  m[2][2]=rr00+rr11+2.0*(+rr01[0][0]-rr01[1][1]+rr01[2][2]);
  m[3][3]=rr00+rr11+2.0*(+rr01[0][0]+rr01[1][1]-rr01[2][2]);
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

  vector<double> eigenvals;
  Matrix<double> eigenvecs;
  int diagerror=diagMat(m, eigenvals, eigenvecs );

  if (diagerror!=0){
    string sdiagerror;
    Tools::convert(diagerror,sdiagerror);
    string msg="DIAGONALIZATION FAILED WITH ERROR CODE "+sdiagerror;
    plumed_merror(msg);
  }

  dist=eigenvals[0];

  Matrix<double> ddist_dm(4,4);

  Vector4d q(eigenvecs[0][0],eigenvecs[0][1],eigenvecs[0][2],eigenvecs[0][3]);

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

  double prefactor=2.0/norm;
  Vector shift=cpositions-matmul(rotation,creference);

  if(!squared){
    dist=sqrt(dist);
    prefactor*=0.5/dist;
  }

// third expensive loop: derivatives
  for(unsigned iat=0;iat<n;iat++){
// there is no need for derivatives of rotation and shift here as it is by construction zero
// (similar to Hellman-Feynman forces)
    derivatives[iat]= prefactor*align[iat]*(positions[iat]-shift - matmul(rotation,reference[iat]));
  }

  return dist;
}

}
