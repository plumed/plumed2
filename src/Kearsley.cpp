#include "Kearsley.h"
#include <cassert>
#include <cmath>
#include <iostream>

using namespace std;
using namespace PLMD;

Kearsley::Kearsley(const vector<Vector> &p1, const vector<Vector> &p2, const vector<double> &align, Log &log):log(log) {
	// copy the structure
	this->p1=p1;
	this->p2=p2;
	com1_is_removed=false;
	com2_is_removed=false;
	this->align=align;
	// now make an initial allocation
	int n=p1.size();
	// eventually here one should make a "hard" resize of all the structures
//	log.printf("Reallocating a size of %d atoms to kearsley structure\n",n);

//	try{
//		diff.resize(n);
//	}catch(bad_alloc&) {
//		cerr<<"Cannot allocate the vector in Kearsley"<<endl;
//	}
//	try{
//		diff.resize(n);
//	}catch(bad_alloc&) {
//		cerr<<"Cannot allocate the vector in Kearsley"<<endl;
//	}

};
// do the alignment
double Kearsley::calculate() {
	// basic sanity check
	if(p1.size()!=p2.size() || p2.size()!=align.size()){
			cerr<<"Kearsley: looks like you have not properly allocated the vectors: the two frames have different size"<<endl;
			cerr<<"size of p1 is :"<<p1.size()<<endl;
			cerr<<"size of p2 is :"<<p2.size()<<endl;
			cerr<<"size of align is :"<<align.size()<<endl;
			exit(0);
	}
	if(p1.size()==0 || p2.size()==0  ){
		cerr<<"Kearsley: looks like you have not properly allocated the vectors: they do not contain anything"<<endl;
		exit(0);
	}
	double m[4][4],rr1[4],rr0[4],q[4],lambda[4],dddq[4][4][4],gamma[3][3][3],rrsq;
	double dm_r1[4][4][3],dm_r0[4][4][3];
	double pi1[3][3],pi0[3][3];
	double xx[3], totalign=0., totdisplace, s, tmp1, fact1, fact2;
	int i,j,k,l,ii,ll,jj,mm,n,nn,iii;

	vector<int> alignmap; // on the fly map done for optimization

	int natoms=p1.size();


	// calculate coms

	for(i=0;i<natoms;i++){
		totalign+=align[i];
		if (align[i]!=0.)alignmap.push_back(i);
	}


	bool do_center=true; // keep it for legacy code compatibility

	if(com1_is_removed==false){

		for(j=0;j<3;j++)xx[j]=0.;

		if(do_center) {// if you dont need to center no prob...
			for(j=0;j<3;j++){
				for(i=0;i<alignmap.size();i++){
					xx[j]+=p1[i][j]*align[alignmap[i]];
				}
			}
			for(j=0;j<3;j++)xx[j]=xx[j]/(totalign);
		}

		for(j=0;j<3;j++)com1[j]=xx[j];

		if (p1reset.size()==0){p1reset.resize(natoms);}

		for(i=0;i<natoms;i++){
			for(j=0;j<3;j++)p1reset[i][j]=p1[i][j]-xx[j];
		}
		com1_is_removed=true;
	}

	if(com2_is_removed==false){

		for(j=0;j<3;j++)xx[j]=0.;

		if(do_center) {// if you dont need to center no prob...
			for(j=0;j<3;j++){
				for(i=0;i<alignmap.size();i++){
					xx[j]+=p2[i][j]*align[alignmap[i]];
				}
			}
			for(j=0;j<3;j++)xx[j]=xx[j]/(totalign);
		}

		for(j=0;j<3;j++)com2[j]=xx[j];

		if (p2reset.size()==0){p2reset.resize(natoms);}

		for(i=0;i<natoms;i++){
			for(j=0;j<3;j++)p2reset[i][j]=p2[i][j]-xx[j];
		}
		com2_is_removed=true;

	}

	//
	// CLEAN M MATRIX
	for(i=0;i<4;i++){
		for(j=0;j<4;j++){
			m[i][j]=0.;
		}
	}

	// ASSIGN MATRIX ELEMENTS USING ONLY THE ATOMS INVOLVED IN ALIGNMENT

	for(i=0;i<alignmap.size();i++){

		k=alignmap[i];
        tmp1=align[k];

		// adopt scaled coordinates

		rr1[0]=p2reset[k][0]*tmp1;
        rr1[1]=p2reset[k][1]*tmp1;
        rr1[2]=p2reset[k][2]*tmp1;
        rr0[0]=p1reset[k][0]*tmp1;
        rr0[1]=p1reset[k][1]*tmp1;
        rr0[2]=p1reset[k][2]*tmp1;

        rrsq=(pow(rr0[0],2)+pow(rr0[1],2)+pow(rr0[2],2)+pow(rr1[0],2)+pow(rr1[1],2)+pow(rr1[2],2));

        m[0][0] +=  rrsq+2.*(-rr0[0]*rr1[0]-rr0[1]*rr1[1]-rr0[2]*rr1[2]);
        m[1][1] +=  rrsq+2.*(-rr0[0]*rr1[0]+rr0[1]*rr1[1]+rr0[2]*rr1[2]);
        m[2][2] +=  rrsq+2.*(+rr0[0]*rr1[0]-rr0[1]*rr1[1]+rr0[2]*rr1[2]);
        m[3][3] +=  rrsq+2.*(+rr0[0]*rr1[0]+rr0[1]*rr1[1]-rr0[2]*rr1[2]);
        m[0][1] += 2.*(-rr0[1]*rr1[2]+rr0[2]*rr1[1]);
        m[0][2] += 2.*( rr0[0]*rr1[2]-rr0[2]*rr1[0]);
        m[0][3] += 2.*(-rr0[0]*rr1[1]+rr0[1]*rr1[0]);
        m[1][2] -= 2.*( rr0[0]*rr1[1]+rr0[1]*rr1[0]);
        m[1][3] -= 2.*( rr0[0]*rr1[2]+rr0[2]*rr1[0]);
        m[2][3] -= 2.*( rr0[1]*rr1[2]+rr0[2]*rr1[1]);

	};
	m[1][0] = m[0][1];
	m[2][0] = m[0][2];
	m[2][1] = m[1][2];
	m[3][0] = m[0][3];
	m[3][1] = m[1][3];
	m[3][2] = m[2][3];

	// diagonalize the 4x4 matrix


	return 0.;
}
;
void Kearsley::assignP1(const std::vector<Vector> & p1) {
	this->p1=p1;
	com1_is_removed=false;
}
void Kearsley::assignP2(const std::vector<Vector> & p2) {
	this->p2=p2;
	com2_is_removed=false;
}

