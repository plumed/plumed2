#include "OptimalAlignment.h"
#include <cassert>
#include <cmath>
#include <iostream>

using namespace std;
using namespace PLMD;

OptimalAlignment::OptimalAlignment( const  std::vector<double>  & align, const  std::vector<double>  & displace, const std::vector<Vector> & p1, const std::vector<Vector> & p2 , Log &log )
:log(log){
	// copy the structure into place
	this->p1=p1;
	this->p2=p2;
	this->align=align;
	this->displace=displace;
	// kearsley init to null
	mykearsley=NULL;
	// basic check
	if(p1.size() != p2.size()){
		log.printf("THE SIZE OF THE TWO FRAMES TO BE ALIGNED ARE DIFFERENT\n");
    }
	// fast behaviour: if all the alignment and displacement are 1.0 then go for fast 
	fast=true;
	for (unsigned i=0;i<align.size();i++ ){
		if(align[i]!=displace[i])fast=false;
	}
	if (mykearsley==NULL) {
		mykearsley=new Kearsley(p1,p2,align,log);
	}

//	if(fast){
//		log.printf("USING FAST ALIGNMENT\n");
//	}else{
//		log.printf("USING SLOW ALIGNMENT\n");
//	}

};

void OptimalAlignment::assignP1(  const std::vector<Vector> & p1 ){
	this->p1=p1;
	if(mykearsley!=NULL){mykearsley->assignP1(p1);}else{cerr<<"kearsley is not initialized"<<endl; exit(0);}
} 

void OptimalAlignment::assignP2(  const std::vector<Vector> & p2 ){
	this->p2=p2;
	if(mykearsley!=NULL){mykearsley->assignP2(p2);}else{cerr<<"kearsley is not initialized"<<endl; exit(0);}
}

double OptimalAlignment::calculate(const std::vector<Vector> & derivatives){

	// at this point everything should be already in place for calculating the alignment (p1,p2,align)
	mykearsley->calculate();
	if(fast){

		//std::cerr<<"USING FAST"<<std::endl;
	}else{
		std::cerr<<"USING SLOW"<<std::endl;
	}
	// destroy the kearsley object

	return 0.;
}


