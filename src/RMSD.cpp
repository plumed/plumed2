#include "RMSD.h"
#include "PDB.h"
#include "Log.h"
#include <cassert>
#include <cmath>
#include <iostream>

using namespace std;
using namespace PLMD;

void RMSD::setFromPDB(const PDB&pdb, string mytype ){
  myoptimalalignment=NULL; 
  alignmentMethod=SIMPLE; // initialize with the simplest case: no rotation
  if (mytype=="SIMPLE"){ 	alignmentMethod=SIMPLE; log.printf("RMSD IS DONE WITH SIMPLE METHOD(NO ROTATION)\n")
;}
  else if (mytype=="OPTIMAL"){ 	alignmentMethod=OPTIMAL; log.printf("RMSD IS DONE WITH OPTIMAL ALIGNMENT METHOD\n"); }
  else assert(0);
  setReference(pdb.getPositions());
  setAlign(pdb.getOccupancy());
  setDisplace(pdb.getBeta());
}

void RMSD::clear(){
  reference.clear();
  align.clear();
  displace.clear();
}

string RMSD::getMethod(){
	string mystring;
	switch(alignmentMethod){
		case SIMPLE: mystring.assign("SIMPLE");break; 
		case OPTIMAL: mystring.assign("OPTIMAL");break; 
	}	
	return mystring;
}

void RMSD::setReference(const vector<Vector> & reference){
  unsigned n=reference.size();
  this->reference=reference;
  assert(align.size()==0);
  assert(displace.size()==0);
  align.resize(n,1.0);
  displace.resize(n,1.0);
}

void RMSD::setAlign(const vector<double> & align){
  assert(this->align.size()==align.size());
  this->align=align;
}

void RMSD::setDisplace(const vector<double> & displace){
  assert(this->displace.size()==displace.size());
  this->displace=displace;
}

double RMSD::calculate(const std::vector<Vector> & positions,std::vector<Vector> &derivatives){

  const unsigned n=reference.size();

  double ret=0.;

  switch(alignmentMethod){
	case SIMPLE:
		//	do a simple alignment without rotation 
		ret=simpleAlignment(align,displace,positions,reference,log,derivatives);
		break;	
	case OPTIMAL:
		if (myoptimalalignment==NULL){ // do full initialization	
			myoptimalalignment=new OptimalAlignment(align,displace,positions,reference,log);
        }
		// this changes the P1 according the running frame
		(*myoptimalalignment).assignP0(positions);

		ret=(*myoptimalalignment).calculate(true, derivatives);
		//(*myoptimalalignment).weightedFindiffTest(false);

		break;	
  }	

  return ret;

}

double RMSD::simpleAlignment(const  std::vector<double>  & align,
		                     const  std::vector<double>  & displace,
		                     const std::vector<Vector> & positions,
		                     const std::vector<Vector> & reference ,
		                     Log &log,
		                     std::vector<Vector>  & derivatives) {
	  double dist(0);
	  double norm(0);
	  int n=reference.size();
	  for(unsigned i=0;i<n;i++){
	      Vector d=delta(reference[i],positions[i]);
	      derivatives[i]=2.0*d;
	      dist+=displace[i]*d.modulo2();
	      norm+=displace[i];
      }

	// sqrt and normalization
     double ret=sqrt(dist/norm);
	///// sqrt and normalization on derivatives
	  for(unsigned i=0;i<n;i++){derivatives[i]*=(0.5/ret/norm);}
	  return ret;
}
