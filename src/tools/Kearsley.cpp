/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2016 The plumed team
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
#include "Kearsley.h"
#include <cmath>
#include <iostream>
#include <cstdlib>
#include "Matrix.h"
#include "Tensor.h"
#include "Log.h"
#include "Matrix.h"
#include "Random.h"

using namespace std;
namespace PLMD{

// put some notes

Kearsley::Kearsley(const vector<Vector> &p0, const vector<Vector> &p1, const vector<double> &align, Log* &log):
	log(log),
	p0(p0),
	p1(p1),
	align(align),
	com0_is_removed(false),
	com1_is_removed(false),
	err(0.0)
{
	// now make an initial allocation
//	int n=p0.size();
	// eventually here one should make a "hard" resize of all the structures
//	log->printf("Reallocating a size of %d atoms to kearsley structure\n",n);

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

}
// do the alignment

double Kearsley::calculate(bool rmsd) {
	// just an ad-hoc scaling factor
	double myscale=1.;
	// basic sanity check
	if(p0.size()!=p1.size() || p1.size()!=align.size()){
			cerr<<"Kearsley: looks like you have not properly allocated the vectors: the two frames have different size"<<endl;
			cerr<<"size of p0 is :"<<p0.size()<<endl;
			cerr<<"size of p1 is :"<<p1.size()<<endl;
			cerr<<"size of align is :"<<align.size()<<endl;
			exit(0);
	}
	if(p0.empty() || p1.empty()  ){
		cerr<<"Kearsley: looks like you have not properly allocated the vectors: they do not contain anything"<<endl;
		exit(0);
	}
        Vector rr1,rr0;
        Vector4d q;
	double dddq[4][4][4],gamma[3][3][3],rrsq;
	Matrix<double> m=Matrix<double>(4,4);
//	double dm_r1[4][4][3],dm_r0[4][4][3];
        Vector dm_r1[4][4];
        Vector dm_r0[4][4];
        Tensor pi1,pi0;
	Tensor d;
	Tensor dinv;
        Vector xx;
	double totalign=0.,  s, tmp1, err;
	unsigned  i,j,k,l,ii,ll,jj,mm,n,nn,iii;
	bool verbose=false;

	vector<int> alignmap; // on the fly map done for optimization

	unsigned natoms=p0.size();


	// calculate coms

	for(i=0;i<natoms;i++){
		if (align[i]>0.){
			alignmap.push_back(i);
			totalign+=align[i];
		}
		if (align[i]<0.){cerr<<"FOUND ALIGNMENT WEIGHT NEGATIVE!"<<endl;exit(0);};

	}


	// later will be implemented something for optimizing this piece of crap

	com0_is_removed=false;
	com1_is_removed=false;

	bool do_center=true; // keep it for legacy code compatibility

	if(com0_is_removed==false){

		xx.zero();

		if(do_center) {// if you dont need to center no prob...
			for(i=0;i<alignmap.size();i++){
				xx+=p0[alignmap[i]]*align[alignmap[i]];
			}
			xx/=(totalign);
		}

		com0=xx;

		if (p0reset.empty()){p0reset.resize(natoms);}
		for(i=0;i<natoms;i++){
			p0reset[i]=p0[i]-xx;
		}
		com0_is_removed=true;

		if (verbose){
			log->printf("P0 RESET\n");
			for(i=0;i<natoms;i++){
				log->printf("ATOM %6u  C   ALA     1    %8.3f%8.3f%8.3f  1.00  1.00\n",i+1,p0reset[i][0]/myscale,p0reset[i][1]/myscale,p0reset[i][2]/myscale);
			}
			log->printf("END\n");
		}

	}

	if(com1_is_removed==false){

		xx.zero();

		if(do_center) {// if you dont need to center no prob...
			for(i=0;i<alignmap.size();i++){
				xx+=p1[alignmap[i]]*align[alignmap[i]];
			}
			xx/=(totalign);
		}

		com1=xx;

		if (p1reset.empty()){p1reset.resize(natoms);}

		for(i=0;i<natoms;i++){
				p1reset[i]=p1[i]-xx;
		}
		com1_is_removed=true;

		if(verbose){
			log->printf("P1 RESET\n");
			for(i=0;i<natoms;i++){
				log->printf("ATOM %6u  C   ALA     1    %8.3f%8.3f%8.3f  1.00  1.00\n",i+1,p1reset[i][0]/myscale,p1reset[i][1]/myscale,p1reset[i][2]/myscale);
			}
			log->printf("END\n");
		}

	}

	bool fake=false;
	if(fake){
	// case of trivial alignment
         // set rotmat
                rotmat0on1=Tensor::identity();
		if (p0reset.size()==0){p0reset.resize(natoms);}
		if (p1reset.size()==0){p1reset.resize(natoms);}

		derrdp0.resize(natoms);
		derrdp1.resize(natoms);
		dmatdp0.resize(3*3*3*natoms);for(i=0;i<dmatdp0.size();i++)dmatdp0[i]=0.;
		dmatdp1.resize(3*3*3*natoms);for(i=0;i<dmatdp1.size();i++)dmatdp1[i]=0.;

		err=0.;
		for(i=0;i<natoms;i++){
			if(align[i]>0.)err+=align[i]*modulo2(p0reset[i]-p1reset[i]);
		}

		return 0.;
	}
	//
	// CLEAN M MATRIX

	m=0.;

	// ASSIGN MATRIX ELEMENTS USING ONLY THE ATOMS INVOLVED IN ALIGNMENT

	for(i=0;i<alignmap.size();i++){

		k=alignmap[i];
		tmp1=sqrt(align[k]/totalign);

		// adopt scaled coordinates

		rr1=p1reset[k]*tmp1;
		rr0=p0reset[k]*tmp1;

		rrsq=modulo2(rr0)+modulo2(rr1);

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

	vector<double> eigenvals;
	Matrix<double> eigenvecs;

	int diagerror=diagMat(m, eigenvals, eigenvecs );

	if (diagerror!=0){cerr<<"DIAGONALIZATION FAILED WITH ERROR CODE "<<diagerror<<endl;exit(0);}

	s=1.0;
	if(eigenvecs(0,0)<0.)s=-1.;//correct for negative values (?)
	// eigenvecs are in rows!!

	q[0]=s*eigenvecs(0,0);
	q[1]=s*eigenvecs(0,1);
	q[2]=s*eigenvecs(0,2);
	q[3]=s*eigenvecs(0,3);
	err=eigenvals[0];

	//log->printf(" ERR: %20.10f \n",err);

	if(verbose){
		log->printf(" ERR: %f \n",err);
		for (i=0;i<4;i++){
			log->printf(" EIGENVALS: %f \n",eigenvals[i]);
		}
	}

	if(abs(eigenvals[0]-eigenvals[1])<1.e-8){
		cerr<<"DIAGONALIZATION: NON UNIQUE SOLUTION"<<endl;exit(0);
	}

	/*
	 * the ROTATION matrix
	 */

	d[0][0]=q[0]*q[0]+q[1]*q[1]-q[2]*q[2]-q[3]*q[3]       ;
	d[1][0]=2.0*(q[1]*q[2]-q[0]*q[3]);
	d[2][0]=2.0*(q[1]*q[3]+q[0]*q[2]);
	d[0][1]=2.0*(q[1]*q[2]+q[0]*q[3]);
	d[1][1]=q[0]*q[0]+q[2]*q[2]-q[1]*q[1]-q[3]*q[3];
	d[2][1]=2.0*(q[2]*q[3]-q[0]*q[1]);
	d[0][2]=2.0*(q[1]*q[3]-q[0]*q[2]);
	d[1][2]=2.0*(q[2]*q[3]+q[0]*q[1]);
	d[2][2]=q[0]*q[0]+q[3]*q[3]-q[1]*q[1]-q[2]*q[2];

	/*
	 * first derivative in perturbation theory : derivative of the rotation matrix respect to the
	 * quternion vectors
	 */

	dddq[0][0][0]= 2.0*q[0];
	dddq[1][0][0]=-2.0*q[3];
	dddq[2][0][0]= 2.0*q[2];
	dddq[0][1][0]= 2.0*q[3];
	dddq[1][1][0]= 2.0*q[0];
	dddq[2][1][0]=-2.0*q[1];
	dddq[0][2][0]=-2.0*q[2];
	dddq[1][2][0]= 2.0*q[1];
	dddq[2][2][0]= 2.0*q[0];

	dddq[0][0][1]= 2.0*q[1];
	dddq[1][0][1]= 2.0*q[2];
	dddq[2][0][1]= 2.0*q[3];
	dddq[0][1][1]= 2.0*q[2];
	dddq[1][1][1]=-2.0*q[1];
	dddq[2][1][1]=-2.0*q[0];
	dddq[0][2][1]= 2.0*q[3];
	dddq[1][2][1]= 2.0*q[0];
	dddq[2][2][1]=-2.0*q[1];

	dddq[0][0][2]=-2.0*q[2];
	dddq[1][0][2]= 2.0*q[1];
	dddq[2][0][2]= 2.0*q[0];
	dddq[0][1][2]= 2.0*q[1];
	dddq[1][1][2]= 2.0*q[2];
	dddq[2][1][2]= 2.0*q[3];
	dddq[0][2][2]=-2.0*q[0];
	dddq[1][2][2]= 2.0*q[3];
	dddq[2][2][2]=-2.0*q[2];

	dddq[0][0][3]=-2.0*q[3];
	dddq[1][0][3]=-2.0*q[0];
	dddq[2][0][3]= 2.0*q[1];
	dddq[0][1][3]= 2.0*q[0];
	dddq[1][1][3]=-2.0*q[3];
	dddq[2][1][3]= 2.0*q[2];
	dddq[0][2][3]= 2.0*q[1];
	dddq[1][2][3]= 2.0*q[2];
	dddq[2][2][3]= 2.0*q[3];

	/*
	 * Build gamma 3x3x3 matrix
	 */
	for(i=0;i<3;i++){     //direction
		for(j=0;j<3;j++){     //direction
			for(k=0;k<3;k++){     //eigenvector number
				gamma[i][j][k]=0.0;
				for(l=0;l<4;l++){   //components of each eigenvector in pert. series
					if(abs(eigenvals[0]-eigenvals[k+1])<1.e-8){
						log->printf(" FOUND DEGENERACY IN RMSD_ESS ROUTINE \n");
						log->printf(" I'm DYING....\n");
						log->printf(" COPYING STACK HERE \n");
						log->printf(" P0\n");
						for(ll=0;ll<natoms;ll++)log->printf(" %f %f %f \n",p0reset[ll][0],p0reset[ll][1],p0reset[ll][2]);
						log->printf(" P1\n");
						for(ll=0;ll<natoms;ll++)log->printf(" %f %f %f \n",p1reset[ll][0],p1reset[ll][1],p1reset[ll][2]);
						exit(0);
					}
					else{
						gamma[i][j][k]  +=  dddq[i][j][l]*eigenvecs(k+1,l)/(eigenvals[0]-eigenvals[k+1]);
					}
				}
                //log->printf("GAMMA %2d %2d %2d V %12.6f\n",i,j,k,gamma[i][j][k]);
			}
		}
	}

	// allocate various arrays

	dmatdp1.resize(3*3*3*natoms);
	for(i=0;i<dmatdp1.size();i++)dmatdp1[i]=0.;
	dmatdp0.resize(3*3*3*natoms);
	for(i=0;i<dmatdp0.size();i++)dmatdp0[i]=0.;

	vector<double> dd_dr_temp;dd_dr_temp.resize(natoms);

	vector<Vector> derr_dr1;
	derr_dr1.resize(natoms);
	vector<Vector> derr_dr0;
	derr_dr0.resize(natoms);
	vector<Vector> array_3_n;
	array_3_n.resize(natoms);


	/*
	 * Table of Derivative of the quaternion matrix respect to atom position: needed only if simple
	 * alignment is required and no correction respect to the rotation matrix is wanted
	 */

	for(iii=0;iii<alignmap.size();iii++){

		i=alignmap[iii];
		tmp1=sqrt(align[i]/totalign);

		// once again: derivative respect to scaled distance

		rr1=2.*p1reset[i]*tmp1;
		rr0=2.*p0reset[i]*tmp1;


		dm_r1 [0][0][0]=(rr1[0]-rr0[0]);
		dm_r1 [0][0][1]=(rr1[1]-rr0[1]);
		dm_r1 [0][0][2]=(rr1[2]-rr0[2]);

		dm_r1 [0][1][0]=0.;
		dm_r1 [0][1][1]= rr0[2];
		dm_r1 [0][1][2]=-rr0[1];

		dm_r1 [0][2][0]=-rr0[2];
		dm_r1 [0][2][1]= 0.;
		dm_r1 [0][2][2]= rr0[0];

		dm_r1 [0][3][0]= rr0[1];
		dm_r1 [0][3][1]=-rr0[0];
		dm_r1 [0][3][2]= 0.;

		dm_r1 [1][1][0]=(rr1[0]-rr0[0]);
		dm_r1 [1][1][1]=(rr1[1]+rr0[1]);
		dm_r1 [1][1][2]=(rr1[2]+rr0[2]);

		dm_r1 [1][2][0]=-rr0[1];
		dm_r1 [1][2][1]=-rr0[0];
		dm_r1 [1][2][2]= 0.;

		dm_r1 [1][3][0]=-rr0[2];
		dm_r1 [1][3][1]= 0.;
		dm_r1 [1][3][2]=-rr0[0];

		dm_r1 [2][2][0]=(rr1[0]+rr0[0]);
		dm_r1 [2][2][1]=(rr1[1]-rr0[1]);
		dm_r1 [2][2][2]=(rr1[2]+rr0[2]);

		dm_r1 [2][3][0]=0.;
		dm_r1 [2][3][1]=-rr0[2];
		dm_r1 [2][3][2]=-rr0[1];

		dm_r1 [3][3][0]=(rr1[0]+rr0[0]);
		dm_r1 [3][3][1]=(rr1[1]+rr0[1]);
		dm_r1 [3][3][2]=(rr1[2]-rr0[2]);
		/*
   derivative respec to to the other vector
		 */
		dm_r0 [0][0][0]=-(rr1[0]-rr0[0]);
		dm_r0 [0][0][1]=-(rr1[1]-rr0[1]);
		dm_r0 [0][0][2]=-(rr1[2]-rr0[2]);

		dm_r0 [0][1][0]=0.       ;
		dm_r0 [0][1][1]=-rr1[2];
		dm_r0 [0][1][2]=rr1[1];

		dm_r0 [0][2][0]= rr1[2];
		dm_r0 [0][2][1]= 0.;
		dm_r0 [0][2][2]=-rr1[0];

		dm_r0 [0][3][0]=-rr1[1] ;
		dm_r0 [0][3][1]= rr1[0];
		dm_r0 [0][3][2]= 0.;

		dm_r0 [1][1][0]=-(rr1[0]-rr0[0]);
		dm_r0 [1][1][1]=(rr1[1]+rr0[1]);
		dm_r0 [1][1][2]=(rr1[2]+rr0[2]);

		dm_r0 [1][2][0]=-rr1[1];
		dm_r0 [1][2][1]=-rr1[0];
		dm_r0 [1][2][2]= 0.;

		dm_r0 [1][3][0]=-rr1[2];
		dm_r0 [1][3][1]= 0.;
		dm_r0 [1][3][2]=-rr1[0];

		dm_r0 [2][2][0]=(rr1[0]+rr0[0]);
		dm_r0 [2][2][1]=-(rr1[1]-rr0[1]);
		dm_r0 [2][2][2]=(rr1[2]+rr0[2]);

		dm_r0 [2][3][0]=0.;
		dm_r0 [2][3][1]=-rr1[2];
		dm_r0 [2][3][2]=-rr1[1];

		dm_r0 [3][3][0]=(rr1[0]+rr0[0]);
		dm_r0 [3][3][1]=(rr1[1]+rr0[1]);
		dm_r0 [3][3][2]=-(rr1[2]-rr0[2]);
		/*
		 * write the diagonal
		 */

		dm_r1[1][0]=dm_r1[0][1];
		dm_r1[2][0]=dm_r1[0][2];
		dm_r1[3][0]=dm_r1[0][3];
		dm_r1[2][1]=dm_r1[1][2];
		dm_r1[3][1]=dm_r1[1][3];
		dm_r1[3][2]=dm_r1[2][3];

		dm_r0[1][0]=dm_r0[0][1];
		dm_r0[2][0]=dm_r0[0][2];
		dm_r0[3][0]=dm_r0[0][3];
		dm_r0[2][1]=dm_r0[1][2];
		dm_r0[3][1]=dm_r0[1][3];
		dm_r0[3][2]=dm_r0[2][3];


		//log->printf("DMDR0 ALIGN %f AT %d VAL %f\n",align[i],i,dm_r0[0][0][j]);


		/*
		 * pi matrix : coefficents in per theory
		 */

        pi0.zero();
		pi1.zero();
		derr_dr1[i].zero();
		derr_dr0[i].zero();

		for(k=0;k<4;k++){
			for(l=0;l<4;l++){
				derr_dr1[i]+=(q[k]*q[l])*dm_r1[l][k];
				derr_dr0[i]+=(q[k]*q[l])*dm_r0[l][k];
				for(mm=0;mm<3;mm++)for(j=0;j<3;j++){
					pi0[mm][j]+=eigenvecs(mm+1,k)*dm_r0[l][k][j]*q[l];
					pi1[mm][j]+=eigenvecs(mm+1,k)*dm_r1[l][k][j]*q[l];
				};
			};
		};
/*
		derr_dr1[i]/=totalign;
		derr_dr0[i]/=totalign;


*/

		for(j=0;j<3;j++){
			for (k=0;k<3;k++){
				for(l=0;l<3;l++){
					int ind=j*3*3*natoms+k*3*natoms+l*natoms+i;
					dmatdp0[ind]=0.;
					dmatdp1[ind]=0.;
					for(ii=0;ii<3;ii++){
						dmatdp1[ind]+=gamma[j][k][ii]*pi1[ii][l];
						dmatdp0[ind]+=gamma[j][k][ii]*pi0[ii][l];
					}

				}

			}

		}
	}




	// end of the calculation of the derivative of the rotation matrix

	/*
	 * Now correct for center of mass: only if needed
	 *
	 */

	bool comcorr_r1=true;

	if(comcorr_r1){
		for(k=0;k<alignmap.size();k++){
			i=alignmap[k];
			tmp1=sqrt(align[i]/totalign);
			array_3_n[i]=tmp1*derr_dr1[i];
			if(do_center){
				for(jj=0;jj<alignmap.size();jj++){
					j=alignmap[jj];
					array_3_n[i]-=tmp1*(align[j]/totalign)*derr_dr1[j];
				}
			}
		}
		for(k=0;k<alignmap.size();k++){
			i=alignmap[k];
			derr_dr1[i]=array_3_n[i];
		}
	}


	bool do_comcorr_r0=true;
	//
	// correction for r0 frame
	//
	if(do_comcorr_r0){
		for(k=0;k<alignmap.size();k++){
			i=alignmap[k];
			tmp1=sqrt(align[i]/totalign);
			array_3_n[i]=tmp1*derr_dr0[i];
			if(do_center){
				for(jj=0;jj<alignmap.size();jj++){
					j=alignmap[jj];
					array_3_n[i]-=tmp1*(align[j]/totalign)*derr_dr0[j];
				}
			}
		}
		for(k=0;k<alignmap.size();k++){
			i=alignmap[k];
			derr_dr0[i]=array_3_n[i];
		}
	}


	bool do_der_r1=true;
	bool do_der_r0=true;
	bool do_der_rotmat=true;

	if(do_der_r1 && do_der_rotmat){
		for(i=0;i<3;i++){
			for(j=0;j<3;j++){
				for(k=0;k<3;k++){
					for(ll=0;ll<alignmap.size();ll++){
						l=alignmap[ll];
						int ind=i*3*3*natoms+j*3*natoms+k*natoms+l;
						tmp1=sqrt(align[l]/totalign);
						dd_dr_temp[l]=tmp1*dmatdp1[ind];
						if(do_center){
							for(nn=0;nn<alignmap.size();nn++){
								n=alignmap[nn];
								dd_dr_temp[l]-=dmatdp1[ind-l+n]*tmp1*align[n]/totalign;
							}
						}

					}
					for(ll=0;ll<alignmap.size();ll++){
						l=alignmap[ll];
						int ind=i*3*3*natoms+j*3*natoms+k*natoms+l;
						dmatdp1[ind]=dd_dr_temp[l];
					}

				}
			}
		}

	}

	if(do_der_r0 && do_der_rotmat){
		for(i=0;i<3;i++){
			for(j=0;j<3;j++){
				for(k=0;k<3;k++){
					for(ll=0;ll<alignmap.size();ll++){
						l=alignmap[ll];
						int ind=i*3*3*natoms+j*3*natoms+k*natoms+l;
						tmp1=sqrt(align[l]/totalign);
						dd_dr_temp[l]=tmp1*dmatdp0[ind];
						if(do_center){
							for(nn=0;nn<alignmap.size();nn++){
								n=alignmap[nn];
								dd_dr_temp[l]-=dmatdp0[ind-l+n]*tmp1*align[n]/totalign;
							}
						}
					}
					for(ll=0;ll<alignmap.size();ll++){
						l=alignmap[ll];
						int ind=i*3*3*natoms+j*3*natoms+k*natoms+l;
						dmatdp0[ind]=dd_dr_temp[l];
					}
				}
			}
		}
	}


	bool do_p1rotated=true;
	if (do_p1rotated){
		// resize if not allocated

		if(p1.size()!=p1rotated.size())p1rotated.resize(p1.size());

//		exit(0);

		for(i=0;i<natoms;i++) p1rotated[i]=matmul(d,p1reset[i]);

		// reallocate difference vectors
		if(p1.size()!=diff1on0.size())diff1on0.resize(p1.size());
		for(i=0;i<natoms;i++){
			diff1on0[i]=p1rotated[i]-p0reset[i];
		}

		if(verbose){
			log->printf("P1-RESET-AND-ROTATED\n");
			for(i=0;i<natoms;i++){
				log->printf("ATOM %6u  C   ALA     2    %8.3f%8.3f%8.3f  1.00  1.00\n",i+1,p1rotated[i][0]/myscale,p1rotated[i][1]/myscale,p1rotated[i][2]/myscale);
			}
			log->printf("END\n");
			log->printf("P0-RESET\n");
			for(i=0;i<natoms;i++){
				log->printf("ATOM %6u  C   ALA     2    %8.3f%8.3f%8.3f  1.00  1.00\n",i+1,p0reset[i][0]/myscale,p0reset[i][1]/myscale,p0reset[i][2]/myscale);
			}
			log->printf("END\n");
		}

	}



        dinv=inverse(d);

	bool do_p0rotated=true;
	if (do_p0rotated){
		if(p0.size()!=p0rotated.size())p0rotated.resize(p0.size());
		for(i=0;i<natoms;i++) p0rotated[i]=matmul(dinv,p0reset[i]);
		if(p1.size()!=diff0on1.size())diff0on1.resize(p1.size());
		for(i=0;i<natoms;i++) diff0on1[i]=p0rotated[i]-p1reset[i];
		if(verbose){
			log->printf("P0-RESET AND INVERSE ROTATED\n");
			for(i=0;i<natoms;i++){
				log->printf("ATOM %6u  C   ALA     1    %8.3f%8.3f%8.3f  1.00  1.00\n",i+1,p0rotated[i][0]/myscale,p0rotated[i][1]/myscale,p0rotated[i][2]/myscale);
			}
			log->printf("END\n");
			log->printf("P1-RESET\n");
			for(i=0;i<natoms;i++){
				log->printf("ATOM %6u  C   ALA     2    %8.3f%8.3f%8.3f  1.00  1.00\n",i+1,p1reset[i][0]/myscale,p1reset[i][1]/myscale,p1reset[i][2]/myscale);
			}
			log->printf("END\n");
		}
	}
	// copy on the official vectors:
	rotmat0on1=d;
	rotmat1on0=dinv;
	derrdp0.resize(natoms);
	derrdp0=derr_dr0;
	derrdp1.resize(natoms);
	derrdp1=derr_dr1;

	// now rescale accordingly for rmsd instead of msd
	if(rmsd){
		err=sqrt(err);
		double tmp=0.5/err;
		for(ii=0;ii<alignmap.size();ii++){
				i=alignmap[ii];
				derrdp0[i]*=tmp;
				derrdp1[i]*=tmp;
		}

	}

	return err;

}
void Kearsley::assignP1(const std::vector<Vector> & p1) {
	this->p1=p1;
	com1_is_removed=false;
}
void Kearsley::assignP0(const std::vector<Vector> & p0) {
	this->p0=p0;
	com0_is_removed=false;
}

void Kearsley::assignAlign(const std::vector<double> & align) {
	this->align=align;
}

void Kearsley::finiteDifferenceInterface(bool rmsd){
log->printf("Entering rmsd finite difference test system for kearsley\n");
log->printf("-------------------------------------------\n");
log->printf("TEST1: derivative of the value (derr_dr0/derr_dr1)\n");
//// test 1
unsigned i,j,l,m;
double step=1.e-6,olderr,delta;
// messing up a bit with align weights
double delta1;
vector<double> align1;
align1.resize(p0.size());
Random rnd;
for (i=0;i<p0.size();i++){
		// draw a random number
	    delta=rnd.RandU01();
	    delta1=rnd.RandU01();
	    if(delta>delta1){
	    //if(delta>0.3){
	    	align1[i]=delta;
	    }else{align1[i]=0.;};
	   // log->printf("ALIGN %d IS %8.3f\n",i,align1[i]);
}
assignAlign(align1);
//// get initial value of the error and derivative of it
olderr=calculate(rmsd);
log->printf("INITIAL ERROR VALUE: %e\n",olderr);
// store the matrix
Tensor old_rotmat0on1=rotmat0on1;

//// get initial value of the error and derivative of it

log->printf("TESTING: derrdp1 \n");
for(unsigned j=0;j<3;j++){
   for(unsigned i=0;i<derrdp1.size();i++){
       // random displacement
       delta=(rnd.RandU01()-0.5)*2*step;
       p1[i][j]+=delta;
	   com1_is_removed=false; // this is required whenever the assignment is not done with the methods
       com0_is_removed=false; // this is required whenever the assignment is not done with the methods
       err=calculate(rmsd);
       //log->printf("INITIAL ERROR VALUE: %e NEW ERROR %e DELTA %e ELEM %d %d \n",olderr,err,delta,i,j );
       p1[i][j]-=delta;
       switch(j){
         case 0:
             log->printf("TESTING: X  %4u ANAL %18.9f NUMER %18.9f DELTA %18.9f ALIGN %6.2f\n",i,derrdp1[i][j],(err-olderr)/delta,derrdp1[i][j]-(err-olderr)/delta,align[i]);break;
         case 1:
             log->printf("TESTING: Y  %4u ANAL %18.9f NUMER %18.9f DELTA %18.9f ALIGN %6.2f\n",i,derrdp1[i][j],(err-olderr)/delta,derrdp1[i][j]-(err-olderr)/delta,align[i]);break;
         case 2:
             log->printf("TESTING: Z  %4u ANAL %18.9f NUMER %18.9f DELTA %18.9f ALIGN %6.2f\n",i,derrdp1[i][j],(err-olderr)/delta,derrdp1[i][j]-(err-olderr)/delta,align[i]);break;

       }
   }
}
//exit(0);
log->printf("TESTING: derrdp0 \n");
for(unsigned j=0;j<3;j++){
   for(unsigned i=0;i<derrdp0.size();i++){
       // random displacement
       delta=(rnd.RandU01()-0.5)*2*step;
       p0[i][j]+=delta;
       com0_is_removed=false; // this is required whenever the assignment is not done with the methods
       com1_is_removed=false; // this is required whenever the assignment is not done with the methods

       err=calculate(rmsd);
       p0[i][j]-=delta;
       switch(j){
         case 0:
             log->printf("TESTING: X  %4u ANAL %18.9f NUMER %18.9f DELTA %18.9f ALIGN %6.2f\n",i,derrdp0[i][j],(err-olderr)/delta,derrdp0[i][j]-(err-olderr)/delta,align[i]);break;
         case 1:
             log->printf("TESTING: Y  %4u ANAL %18.9f NUMER %18.9f DELTA %18.9f ALIGN %6.2f\n",i,derrdp0[i][j],(err-olderr)/delta,derrdp0[i][j]-(err-olderr)/delta,align[i]);break;
         case 2:
             log->printf("TESTING: Z  %4u ANAL %18.9f NUMER %18.9f DELTA %18.9f ALIGN %6.2f\n",i,derrdp0[i][j],(err-olderr)/delta,derrdp0[i][j]-(err-olderr)/delta,align[i]);break;
       }
   }
}

log->printf("TESTING: dmatdp0 \n");
for(l=0;l<3;l++){
  for(m=0;m<3;m++){
    for(j=0;j<3;j++){
       for(i=0;i<p0.size();i++){
           // random displacement
           delta=(rnd.RandU01()-0.5)*2*step;
           p0[i][j]+=delta;
           com0_is_removed=false;
           com1_is_removed=false;
           calculate(rmsd);
           p0[i][j]-=delta;
       	   int ind=l*3*3*p0.size()+m*3*p0.size()+j*p0.size()+i;
           switch(j){
           	 case 0:
                log->printf("TESTING: DMATDP0 [ %u ][ %u ]:  X %u ANAL %18.9f NUMER %18.9f DELTA %18.9f ALIGN %6.2f\n",l,m,i,dmatdp0[ind],(rotmat0on1[l][m]- old_rotmat0on1[l][m])/delta,dmatdp0[ind]-(rotmat0on1[l][m]- old_rotmat0on1[l][m])/delta,align[i]);break;

           	 case 1:
                log->printf("TESTING: DMATDP0 [ %u ][ %u ]:  Y %u ANAL %18.9f NUMER %18.9f DELTA %18.9f ALIGN %6.2f\n",l,m,i,dmatdp0[ind],(rotmat0on1[l][m]- old_rotmat0on1[l][m])/delta,dmatdp0[ind]-(rotmat0on1[l][m]- old_rotmat0on1[l][m])/delta,align[i]);break;

           	 case 2:
                log->printf("TESTING: DMATDP0 [ %u ][ %u ]:  Z %u ANAL %18.9f NUMER %18.9f DELTA %18.9f ALIGN %6.2f\n",l,m,i,dmatdp0[ind],(rotmat0on1[l][m]- old_rotmat0on1[l][m])/delta,dmatdp0[ind]-(rotmat0on1[l][m]- old_rotmat0on1[l][m])/delta,align[i]);break;

           }
       }
    }
  }
}
log->printf("TESTING: dmatdp1 \n");
for(l=0;l<3;l++){
  for(m=0;m<3;m++){
    for(j=0;j<3;j++){
       for(i=0;i<p1.size();i++){
           // random displacement
           delta=(rnd.RandU01()-0.5)*2*step;
           p1[i][j]+=delta;
           com0_is_removed=false;
           com1_is_removed=false;
           calculate(rmsd);
           p1[i][j]-=delta;
       	   int ind=l*3*3*p1.size()+m*3*p1.size()+j*p1.size()+i;
           switch(j){

             case 0:
                log->printf("TESTING: DMATDP1 [ %u ][ %u ]:  X %u ANAL %18.9f NUMER %18.9f DELTA %18.9f ALIGN %6.2f\n",l,m,i,dmatdp1[ind],(rotmat0on1[l][m]- old_rotmat0on1[l][m])/delta,dmatdp1[ind]-(rotmat0on1[l][m]- old_rotmat0on1[l][m])/delta,align[i]);break;

             case 1:
                log->printf("TESTING: DMATDP1 [ %u ][ %u ]:  Y %u ANAL %18.9f NUMER %18.9f DELTA %18.9f ALIGN %6.2f\n",l,m,i,dmatdp1[ind],(rotmat0on1[l][m]- old_rotmat0on1[l][m])/delta,dmatdp1[ind]-(rotmat0on1[l][m]- old_rotmat0on1[l][m])/delta,align[i]);break;

             case 2:
                log->printf("TESTING: DMATDP1 [ %u ][ %u ]:  Z %u ANAL %18.9f NUMER %18.9f DELTA %18.9f ALIGN %6.2f\n",l,m,i,dmatdp1[ind],(rotmat0on1[l][m]- old_rotmat0on1[l][m])/delta,dmatdp1[ind]-(rotmat0on1[l][m]- old_rotmat0on1[l][m])/delta,align[i]);break;

           }
       }
    }
  }
}

	exit(0);
}
}
