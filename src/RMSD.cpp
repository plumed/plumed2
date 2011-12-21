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
  log.printf("RMSD: begin calculated\n");
  switch(alignmentMethod){
	case SIMPLE:
		//	do a simple alignment without rotation 
		break;	
	case OPTIMAL:
		if (myoptimalalignment==NULL){ // do full initialization	
			myoptimalalignment=new OptimalAlignment(align,displace,positions,reference,log);
        }
		// this changes the P1 according the running frame
		(*myoptimalalignment).assignP1(positions);
		ret=(*myoptimalalignment).calculate(derivatives);
		break;	
  }	
  log.printf("RMSD: done!\n");

  return ret;
///  bool simple,trivial;
///  simple=true;
///  trivial=true; // means: no alignment!!
///  for(unsigned i=0;i<n;i++) if(align[i]!=1.0) simple=false;
///  for(unsigned i=0;i<n;i++) if(align[i]!=0.0) trivial=false;
///  for(unsigned i=0;i<n;i++) if(displace[i]!=1.0) simple=false;
///
///  double dist(0);
///  double norm(0);
///  if(trivial){
///    for(unsigned i=0;i<n;i++){
///      Vector d=delta(reference[i],positions[i]);
///      derivatives[i]=2.0*d;
///      dist+=displace[i]*d.modulo2();
///      norm+=displace[i];
///    }
///  } else {
///// TODO
///    assert(trivial);
///  }
///
///// sqrt and normalization
///  double ret=sqrt(dist/norm);
///// sqrt and normalization on derivatives
///  for(unsigned i=0;i<n;i++){derivatives[i].scale(0.5/ret/norm);}

}


