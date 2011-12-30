#include "OptimalAlignment.h"
#include <cassert>
#include <cmath>
#include <iostream>

using namespace std;
using namespace PLMD;

OptimalAlignment::OptimalAlignment( const  std::vector<double>  & align, const  std::vector<double>  & displace, const std::vector<Vector> & p0, const std::vector<Vector> & p1 , Log &log )
:log(log){
	// copy the structure into place
	this->p0=p0;
	this->p1=p1;
	this->align=align;
	this->displace=displace;
	// kearsley init to null
	mykearsley=NULL;
	// basic check
	if(p0.size() != p1.size()){
		log.printf("THE SIZE OF THE TWO FRAMES TO BE ALIGNED ARE DIFFERENT\n");
    }
	// fast behaviour: if all the alignment and displacement are 1.0 then go for fast 
	fast=true;
	for (unsigned i=0;i<align.size();i++ ){
		if(align[i]!=displace[i])fast=false;
	}
	if (mykearsley==NULL) {
		mykearsley=new Kearsley(p0,p1,align,log);
	}

//	if(fast){
//		log.printf("USING FAST ALIGNMENT\n");
//	}else{
//		log.printf("USING SLOW ALIGNMENT\n");
//	}

};

void OptimalAlignment::assignP0(  const std::vector<Vector> & p0 ){
	this->p0=p1;
	if(mykearsley!=NULL){mykearsley->assignP0(p0);}else{cerr<<"kearsley is not initialized"<<endl; exit(0);}
}

void OptimalAlignment::assignP1(  const std::vector<Vector> & p1 ){
	this->p1=p1;
	if(mykearsley!=NULL){mykearsley->assignP1(p1);}else{cerr<<"kearsley is not initialized"<<endl; exit(0);}
}

double OptimalAlignment::calculate( std::vector<Vector> & derivatives){
	bool rmsd=true;
	double err;

	// at this point everything should be already in place for calculating the alignment (p1,p2,align)
	// here everything is done with kearsley algorithm. Extension to other optimal alignment algos is
	// possible here below with a switch

	mykearsley->calculate(rmsd);  // this calculates the MSD: transform into RMSD

	// check findiff alignment
	// mykearsley->finiteDifferenceInterface(rmsd);

	if(fast){
		err=mykearsley->err;
		derrdp0=mykearsley->derrdp0;
		derrdp1=mykearsley->derrdp1;
		derivatives=derrdp0;
	}else{
		std::cerr<<"USING SLOW"<<std::endl;
		log.printf("SLOW ALIGNMENT NOT YET IMPLEMENTED");exit(0);
	}
	// destroy the kearsley object?

	return err;
}


