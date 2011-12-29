#include "Kearsley.h"
#include <cassert>
#include <cmath>
#include <iostream>
#include "Matrix.h"
#include "Tensor.h"

using namespace std;
using namespace PLMD;

// put some notes

Kearsley::Kearsley(const vector<Vector> &p0, const vector<Vector> &p1, const vector<double> &align, Log &log):log(log) {
	// copy the structure
	this->p0=p0;
	this->p1=p1;
	com0_is_removed=false;
	com1_is_removed=false;
	this->align=align;
	// now make an initial allocation
	int n=p0.size();
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
	// just an ad-hoc scaling factor
	double myscale=0.1;
	// basic sanity check
	if(p0.size()!=p1.size() || p1.size()!=align.size()){
			cerr<<"Kearsley: looks like you have not properly allocated the vectors: the two frames have different size"<<endl;
			cerr<<"size of p0 is :"<<p0.size()<<endl;
			cerr<<"size of p1 is :"<<p1.size()<<endl;
			cerr<<"size of align is :"<<align.size()<<endl;
			exit(0);
	}
	if(p0.size()==0 || p1.size()==0  ){
		cerr<<"Kearsley: looks like you have not properly allocated the vectors: they do not contain anything"<<endl;
		exit(0);
	}
	double rr1[4],rr0[4],q[4],dddq[4][4][4],gamma[3][3][3],rrsq;
	Matrix<double> m=Matrix<double>(4,4);
	double dm_r1[4][4][3],dm_r0[4][4][3];
	double pi1[3][3],pi0[3][3];
	Tensor d;
	Tensor dinv;
	double xx[3], totalign=0.,  s, tmp1, err;
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


	bool do_center=true; // keep it for legacy code compatibility

	if(com0_is_removed==false){

		for(j=0;j<3;j++)xx[j]=0.;

		if(do_center) {// if you dont need to center no prob...
			for(j=0;j<3;j++){
				for(i=0;i<alignmap.size();i++){
					xx[j]+=p0[alignmap[i]][j]*align[alignmap[i]];
				}
			}
			for(j=0;j<3;j++)xx[j]=xx[j]/(totalign);
		}

		for(j=0;j<3;j++)com0[j]=xx[j];

		if (p0reset.size()==0){p0reset.resize(natoms);}
		for(i=0;i<natoms;i++){
			for(j=0;j<3;j++)p0reset[i][j]=p0[i][j]-xx[j];
		}
		com0_is_removed=true;

		if (verbose){
			log.printf("P0 RESET\n");
			for(i=0;i<natoms;i++){
				log.printf("ATOM %6d  C   ALA     1    %8.3f%8.3f%8.3f  1.00  1.00\n",i+1,p0reset[i][0]/myscale,p0reset[i][1]/myscale,p0reset[i][2]/myscale);
			}
			log.printf("END\n");
		}

	}

	if(com1_is_removed==false){

		for(j=0;j<3;j++)xx[j]=0.;

		if(do_center) {// if you dont need to center no prob...
			for(j=0;j<3;j++){
				for(i=0;i<alignmap.size();i++){
					xx[j]+=p1[alignmap[i]][j]*align[alignmap[i]];
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

		if(verbose){
			log.printf("P1 RESET\n");
			for(i=0;i<natoms;i++){
				log.printf("ATOM %6d  C   ALA     1    %8.3f%8.3f%8.3f  1.00  1.00\n",i+1,p1reset[i][0]/myscale,p1reset[i][1]/myscale,p1reset[i][2]/myscale);
			}
			log.printf("END\n");
		}

	}

	//
	// CLEAN M MATRIX

	m=0.;

	// ASSIGN MATRIX ELEMENTS USING ONLY THE ATOMS INVOLVED IN ALIGNMENT

	for(i=0;i<alignmap.size();i++){

		k=alignmap[i];
        tmp1=align[k];

		// adopt scaled coordinates

		rr1[0]=p1reset[k][0]*tmp1;
        rr1[1]=p1reset[k][1]*tmp1;
        rr1[2]=p1reset[k][2]*tmp1;
        rr0[0]=p0reset[k][0]*tmp1;
        rr0[1]=p0reset[k][1]*tmp1;
        rr0[2]=p0reset[k][2]*tmp1;

        rrsq=(pow(rr0[0],2)+pow(rr0[1],2)+pow(rr0[2],2)+pow(rr1[0],2)+pow(rr1[1],2)+pow(rr1[2],2));

        m(0,0) +=  rrsq+2.*(-rr0[0]*rr1[0]-rr0[1]*rr1[1]-rr0[2]*rr1[2]);
        m(1,1) +=  rrsq+2.*(-rr0[0]*rr1[0]+rr0[1]*rr1[1]+rr0[2]*rr1[2]);
        m(2,2) +=  rrsq+2.*(+rr0[0]*rr1[0]-rr0[1]*rr1[1]+rr0[2]*rr1[2]);
        m(3,3) +=  rrsq+2.*(+rr0[0]*rr1[0]+rr0[1]*rr1[1]-rr0[2]*rr1[2]);
        m(0,1) += 2.*(-rr0[1]*rr1[2]+rr0[2]*rr1[1]);
        m(0,2) += 2.*( rr0[0]*rr1[2]-rr0[2]*rr1[0]);
        m(0,3) += 2.*(-rr0[0]*rr1[1]+rr0[1]*rr1[0]);
        m(1,2) -= 2.*( rr0[0]*rr1[1]+rr0[1]*rr1[0]);
        m(1,3) -= 2.*( rr0[0]*rr1[2]+rr0[2]*rr1[0]);
        m(2,3) -= 2.*( rr0[1]*rr1[2]+rr0[2]*rr1[1]);

	};
	m(1,0) = m(0,1);
	m(2,0) = m(0,2);
	m(2,1) = m(1,2);
	m(3,0) = m(0,3);
	m(3,1) = m(1,3);
	m(3,2) = m(2,3);

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
	err=eigenvals[0]/totalign;
	if(verbose){
		log.printf(" ERR: %f \n",err);
		for (i=0;i<4;i++){
			log.printf(" EIGENVALS: %f \n",eigenvals[i]);
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
						log.printf(" FOUND DEGENERACY IN RMSD_ESS ROUTINE \n");
						log.printf(" I'm DYING....\n");
						log.printf(" COPYING STACK HERE \n");
						log.printf(" P0\n");
						for(ll=0;ll<natoms;ll++)log.printf(" %f %f %f \n",p0reset[ll][0],p0reset[ll][1],p0reset[ll][2]);
						log.printf(" P1\n");
						for(ll=0;ll<natoms;ll++)log.printf(" %f %f %f \n",p1reset[ll][0],p1reset[ll][1],p1reset[ll][2]);
						exit(0);
					}
					else{
						gamma[i][j][k]  +=  dddq[i][j][l]*m(l,k+1)/(eigenvals[0]-eigenvals[k+1]);
					}
				}
			}
		}
	}

	// allocate various arrays
	Matrix4d<double> dd_dr1=Matrix4d<double>(3,3,3,natoms);
	Matrix4d<double> dd_dr0=Matrix4d<double>(3,3,3,natoms);
	Matrix4d<double> dd_dr_temp=Matrix4d<double>(3,3,3,natoms);

	Matrix4d<double> dm_r0_store=Matrix4d<double>(3,3,3,natoms);
	Matrix4d<double> dm_r1_store=Matrix4d<double>(3,3,3,natoms);


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
		tmp1=align[i];

		// once again: derivative respect to scaled distance

		rr1[0]=2.*p1reset[i][0]*tmp1;
		rr1[1]=2.*p1reset[i][1]*tmp1;
		rr1[2]=2.*p1reset[i][2]*tmp1;
		rr0[0]=2.*p0reset[i][0]*tmp1;
		rr0[1]=2.*p0reset[i][1]*tmp1;
		rr0[2]=2.*p0reset[i][2]*tmp1;


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

		for(j=0;j<3;j++){

			dm_r1[1][0][j]=dm_r1[0][1][j];
			dm_r1[2][0][j]=dm_r1[0][2][j];
			dm_r1[3][0][j]=dm_r1[0][3][j];
			dm_r1[2][1][j]=dm_r1[1][2][j];
			dm_r1[3][1][j]=dm_r1[1][3][j];
			dm_r1[3][2][j]=dm_r1[2][3][j];

			dm_r0[1][0][j]=dm_r0[0][1][j];
			dm_r0[2][0][j]=dm_r0[0][2][j];
			dm_r0[3][0][j]=dm_r0[0][3][j];
			dm_r0[2][1][j]=dm_r0[1][2][j];
			dm_r0[3][1][j]=dm_r0[1][3][j];
			dm_r0[3][2][j]=dm_r0[2][3][j];

			for(ll=0;ll<4;ll++){
				for(mm=0;mm<4;mm++){
					dm_r0_store(ll,mm,j,i)=dm_r0[ll][mm][j];
					dm_r1_store(ll,mm,j,i)=dm_r1[ll][mm][j];
				};
			};
		};

		/*
		 * pi matrix : coefficents in per theory
		 */

		for(j=0;j<3;j++){
			pi1[0][j]=0.;
			pi1[1][j]=0.;
			pi1[2][j]=0.;
			pi0[0][j]=0.;
			pi0[1][j]=0.;
			pi0[2][j]=0.;

			derr_dr1[i][j]=0.;
			derr_dr0[i][j]=0.;

			for(k=0;k<4;k++){
				for(l=0;l<4;l++){
					derr_dr1[i][j]=derr_dr1[i][j]+q[k]*q[l]*dm_r1[l][k][j];
					derr_dr0[i][j]=derr_dr0[i][j]+q[k]*q[l]*dm_r0[l][k][j];
					for(mm=0;mm<3;mm++){
						pi0[mm][j]+=m(k,mm+1)*dm_r0[l][k][j]*q[l];
						pi1[mm][j]+=m(k,mm+1)*dm_r1[l][k][j]*q[l];
					};
				};
			};
			derr_dr1[i][j]=derr_dr1[i][j]/totalign;
			derr_dr0[i][j]=derr_dr0[i][j]/totalign;

		};


		for(j=0;j<3;j++){
			for (k=0;k<3;k++){
				for(l=0;l<3;l++){
					dd_dr1(j,k,l,i)=0.;
					dd_dr0(j,k,l,i)=0.;
					for(ii=0;ii<3;ii++){
						dd_dr1(j,k,l,i)+=gamma[j][k][ii]*pi1[ii][l];
						dd_dr0(j,k,l,i)+=gamma[j][k][ii]*pi0[ii][l];
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
		for(l=0;l<3;l++){
			for(k=0;k<alignmap.size();k++){
				i=alignmap[k];
				array_3_n[i][l]=align[i]*derr_dr1[i][l];
				tmp1=align[i]/totalign;
				if(do_center){
					for(jj=0;jj<alignmap.size();jj++){
						j=alignmap[jj];
						array_3_n[i][l]-=tmp1*align[j]*derr_dr1[j][l];
					}
				}

			}
		}
		for(l=0;l<3;l++){
			for(k=0;k<alignmap.size();k++){
				i=alignmap[k];
				derr_dr1[i][l]=array_3_n[i][l];
			}
		}
	}


	bool do_comcorr_r0=true;
	//
	// correction for r0 frame
	//
	if(do_comcorr_r0){
		for(l=0;l<3;l++){
			for(k=0;k<alignmap.size();k++){
				i=alignmap[k];
				array_3_n[i][l]=align[i]*derr_dr0[i][l];
				tmp1=align[i]/totalign;
				if(do_center){
					for(jj=0;jj<alignmap.size();jj++){
						j=alignmap[jj];
						array_3_n[i][l]-=tmp1*align[j]*derr_dr0[j][l];
					}
				}
			}
		}
		for(l=0;l<3;l++){
			for(k=0;k<alignmap.size();k++){
				i=alignmap[k];
				derr_dr0[i][l]=array_3_n[i][l];
			}
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
						dd_dr_temp(i,j,k,l)=align[l]*dd_dr1(i,j,k,l);
						tmp1=align[l]/totalign;
						if(do_center){
							for(nn=0;nn<alignmap.size();nn++){
								n=align[nn];
								dd_dr_temp(i,j,k,l)-=dd_dr1(i,j,k,n)*tmp1*align[n];
							}
						}

					}
				}
			}
		}
		for(i=0;i<3;i++){
			for(j=0;j<3;j++){
				for(k=0;k<3;k++){
					for(ll=0;ll<align.size();ll++){
						l=align[ll];
						dd_dr1(i,j,k,l)=dd_dr_temp(i,j,k,l);
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
						dd_dr_temp(i,j,k,l)=align[l]*dd_dr0(i,j,k,l);
						tmp1=align[l]/totalign;
						if(do_center){
							for(nn=0;nn<alignmap.size();nn++){
								n=alignmap[nn];
								dd_dr_temp(i,j,k,l)-=dd_dr0(i,j,k,n)*tmp1*align[n];
							}
						}
					}
				}
			}
		}
		for(i=0;i<3;i++){
			for(j=0;j<3;j++){
				for(k=0;k<3;k++){
					for(ll=0;ll<alignmap.size();ll++){
						l=alignmap[ll];
						dd_dr0(i,j,k,l)=dd_dr_temp(i,j,k,l);
					}
				}
			}
		}
	}

	bool do_p1rotated=true;
	if (do_p1rotated){
		// resize if not allocated
		if(p1.size()!=p1rotated.size())p1rotated.resize(p1.size());

		for(i=0;i<natoms;i++){
			p1rotated[i][0]=d[0][0]*p1reset[i][0]+
					        d[0][1]*p1reset[i][1]+
					        d[0][2]*p1reset[i][2];
			p1rotated[i][1]=d[1][0]*p1reset[i][0]+
					        d[1][1]*p1reset[i][1]+
					        d[1][2]*p1reset[i][2];
			p1rotated[i][2]=d[2][0]*p1reset[i][0]+
					        d[2][1]*p1reset[i][1]+
					        d[2][2]*p1reset[i][2];
		}

		// reallocate difference vectors
		if(p1.size()!=diff.size())diff.resize(p1.size());
		for(i=0;i<natoms;i++){
			diff[i][0]=p1rotated[i][0]-p0reset[i][0];
			diff[i][1]=p1rotated[i][1]-p0reset[i][1];
			diff[i][2]=p1rotated[i][2]-p0reset[i][2];
		}

		if(verbose){
			log.printf("P1-RESET-AND-ROTATED\n");
			for(i=0;i<natoms;i++){
				log.printf("ATOM %6d  C   ALA     2    %8.3f%8.3f%8.3f  1.00  1.00\n",i+1,p1rotated[i][0]/myscale,p1rotated[i][1]/myscale,p1rotated[i][2]/myscale);
			}
			log.printf("END\n");
			log.printf("P0-RESET\n");
			for(i=0;i<natoms;i++){
				log.printf("ATOM %6d  C   ALA     2    %8.3f%8.3f%8.3f  1.00  1.00\n",i+1,p0reset[i][0]/myscale,p0reset[i][1]/myscale,p0reset[i][2]/myscale);
			}
			log.printf("END\n");
		}

	}

	// wonderful supersoviet hardcoded 3x3 matrix inversion ;)
	double det;
	det=         d[0][0]*(d[2][2]*d[1][1]-d[2][1]*d[1][2])
				-d[1][0]*(d[2][2]*d[0][1]-d[2][1]*d[0][2])
				+d[2][0]*(d[1][2]*d[0][1]-d[1][1]*d[0][2]);

	dinv[0][0]= (d[2][2]*d[1][1]-d[2][1]*d[1][2])/det;//a22a11-a21a12
	dinv[0][1]=-(d[2][2]*d[0][1]-d[2][1]*d[0][2])/det;//-(a22a01-a21a02)
	dinv[0][2]= (d[1][2]*d[0][1]-d[1][1]*d[0][2])/det;//a12a01-a11a02

	dinv[1][0]=-(d[2][2]*d[1][0]-d[2][0]*d[1][2])/det;  //-(a22a10-a20a12)
	dinv[1][1]= (d[2][2]*d[0][0]-d[2][0]*d[0][2])/det;	// 	a22a00-a20a02
	dinv[1][2]=-(d[1][2]*d[0][0]-d[1][0]*d[0][2])/det;	// 	-(a12a00-a10a02)

	dinv[2][0]= (d[2][1]*d[1][0]-d[2][0]*d[1][1])/det;  //  a21a10-a20a11
	dinv[2][1]=-(d[2][1]*d[0][0]-d[2][0]*d[0][1])/det;	//	-(a21a00-a20a01)
	dinv[2][2]= (d[1][1]*d[0][0]-d[1][0]*d[0][1])/det;  //a11a00-a10a01

	bool do_p0rotated=true;
	if (do_p0rotated){
		if(p0.size()!=p0rotated.size())p0rotated.resize(p0.size());
		for(i=0;i<natoms;i++){
			p0rotated[i][0]=dinv[0][0]*p0reset[i][0]+
							dinv[0][1]*p0reset[i][1]+
							dinv[0][2]*p0reset[i][2];
			p0rotated[i][1]=dinv[1][0]*p0reset[i][0]+
							dinv[1][1]*p0reset[i][1]+
							dinv[1][2]*p0reset[i][2];
			p0rotated[i][2]=dinv[2][0]*p0reset[i][0]+
							dinv[2][1]*p0reset[i][1]+
							dinv[2][2]*p0reset[i][2];
		}
		if(verbose){
			log.printf("P0-RESET AND INVERSE ROTATED\n");
			for(i=0;i<natoms;i++){
				log.printf("ATOM %6d  C   ALA     1    %8.3f%8.3f%8.3f  1.00  1.00\n",i+1,p0rotated[i][0]/myscale,p0rotated[i][1]/myscale,p0rotated[i][2]/myscale);
			}
			log.printf("END\n");
			log.printf("P1-RESET\n");
			for(i=0;i<natoms;i++){
				log.printf("ATOM %6d  C   ALA     2    %8.3f%8.3f%8.3f  1.00  1.00\n",i+1,p1reset[i][0]/myscale,p1reset[i][1]/myscale,p1reset[i][2]/myscale);
			}
			log.printf("END\n");
		}
	}
	// copy on the official vectors:
	rotmat0on1=d;
	rotmat1on0=dinv;
	derrdp0.resize(natoms);
	derrdp0=derr_dr0;
	derrdp1.resize(natoms);
	derrdp1=derr_dr1;
	dmatdp0.resize(3,3,3,natoms);
	dmatdp0=dd_dr0;
	dmatdp1.resize(3,3,3,natoms);
	dmatdp1=dd_dr1;

	return err;

}
;
void Kearsley::assignP1(const std::vector<Vector> & p1) {
	this->p1=p1;
	com1_is_removed=false;
}
void Kearsley::assignP0(const std::vector<Vector> & p0) {
	this->p0=p0;
	com0_is_removed=false;
}

void Kearsley::finiteDifferenceInterface(){
log.printf("Entering rmsd finite difference test system\n");
log.printf("-------------------------------------------\n");
log.printf("TEST1: derivative of the value (derr_dr0/derr_dr1)\n");
//// test 1
unsigned i,j,k,l,m;
double step=1.e-9,olderr,delta;

//// get initial value of the error and derivative of it
olderr=calculate();
log.printf("INITIAL ERROR VALUE: %e\n",olderr);
//// store the derivative
vector<Vector> old_derrdp0=derrdp0 ;
vector<Vector> old_derrdp1=derrdp1 ;
// store the matrix
Tensor old_rotmat0on1=rotmat0on1,old_rotmat1on0=rotmat1on0;
// store the deriv of matrix respect to atoms
Matrix4d<double>  old_dmatdp0=dmatdp0;
Matrix4d<double>  old_dmatdp1=dmatdp1;

log.printf("TESTING: derrdp1 \n");
for(unsigned j=0;j<3;j++){
   for(unsigned i=0;i<derrdp1.size();i++){
       // random displacement
       delta=(drand48()-0.5)*2*step;
       p1[i][j]+=delta;
	   com1_is_removed=false; // this is required whenever the assignment is not done with the methods
       com0_is_removed=false; // this is required whenever the assignment is not done with the methods
       err=calculate();
       //log.printf("INITIAL ERROR VALUE: %e NEW ERROR %e DELTA %e ELEM %d %d \n",olderr,err,delta,i,j );
       p1[i][j]-=delta;
       switch(j){
         case 0:
             log.printf("TESTING: X  %4d ANAL %18.9f NUMER %18.9f DELTA %18.9f\n",i,derrdp1[i][j],(err-olderr)/delta,derrdp1[i][j]-(err-olderr)/delta);break;
         case 1:
             log.printf("TESTING: Y  %4d ANAL %18.9f NUMER %18.9f DELTA %18.9f\n",i,derrdp1[i][j],(err-olderr)/delta,derrdp1[i][j]-(err-olderr)/delta);break;
         case 2:
             log.printf("TESTING: Z  %4d ANAL %18.9f NUMER %18.9f DELTA %18.9f\n",i,derrdp1[i][j],(err-olderr)/delta,derrdp1[i][j]-(err-olderr)/delta);break;

       }
   }
}
log.printf("TESTING: derrdp0 \n");
for(unsigned j=0;j<3;j++){
   for(unsigned i=0;i<derrdp0.size();i++){
       // random displacement
       delta=(drand48()-0.5)*2*step;
       p0[i][j]+=delta;
       com0_is_removed=false; // this is required whenever the assignment is not done with the methods
       com1_is_removed=false; // this is required whenever the assignment is not done with the methods

       err=calculate();
       p0[i][j]-=delta;
       switch(j){
         case 0:
             log.printf("TESTING: X  %4d ANAL %18.9f NUMER %18.9f DELTA %18.9f\n",i,derrdp0[i][j],(err-olderr)/delta,derrdp0[i][j]-(err-olderr)/delta);break;
         case 1:
             log.printf("TESTING: Y  %4d ANAL %18.9f NUMER %18.9f DELTA %18.9f\n",i,derrdp0[i][j],(err-olderr)/delta,derrdp0[i][j]-(err-olderr)/delta);break;
         case 2:
             log.printf("TESTING: Z  %4d ANAL %18.9f NUMER %18.9f DELTA %18.9f\n",i,derrdp0[i][j],(err-olderr)/delta,derrdp0[i][j]-(err-olderr)/delta);break;

       }
   }
}

log.printf("TESTING: dmatdp0 \n");
for(l=0;l<3;l++){
  for(m=0;m<3;m++){
    for(j=0;j<3;j++){
       for(i=0;i<p0.size();i++){
           // random displacement
           delta=(drand48()-0.5)*2*step;
           p0[i][j]+=delta;
           com0_is_removed=false;
           com1_is_removed=false;
           calculate();
           p0[i][j]-=delta;
           switch(j){
             case 0:
                log.printf("TESTING: DMATDP0 [ %d ][ %d ]:  X %d ANAL %18.9f NUMER %18.9f DELTA %18.9f\n",l,m,i,dmatdp0(l,m,j,i),(rotmat0on1[l][m]- old_rotmat0on1[l][m])/delta,dmatdp0(l,m,j,i)-(rotmat0on1[l][m]- old_rotmat0on1[l][m])/delta);break;
           }
       }
    }
  }
}
log.printf("TESTING: dmatdp1 \n");
for(l=0;l<3;l++){
  for(m=0;m<3;m++){
    for(j=0;j<3;j++){
       for(i=0;i<p1.size();i++){
           // random displacement
           delta=(drand48()-0.5)*2*step;
           p1[i][j]+=delta;
           com0_is_removed=false;
           com1_is_removed=false;
           calculate();
           p1[i][j]-=delta;
           switch(j){
             case 0:
                log.printf("TESTING: DMATDP1 [ %d ][ %d ]:  X %d ANAL %18.9f NUMER %18.9f DELTA %18.9f\n",l,m,i,dmatdp1(l,m,j,i),(rotmat0on1[l][m]- old_rotmat0on1[l][m])/delta,dmatdp1(l,m,j,i)-(rotmat0on1[l][m]- old_rotmat0on1[l][m])/delta);break;
           }
       }
    }
  }
}



//fprintf(mtd_data.fplog,"TESTING: dd_dr1 \n");
//for(l=0;l<3;l++){
//  for(m=0;m<3;m++){
//    for(j=0;j<3;j++){
//       for(i=0;i<inpack.natoms;i++){
//           // random displacement
//           delta=(drand48()-0.5)*2*step;
//           inpack.r1[j][i]+=delta;
//           rmsd_mini_pack(inpack,work,7,2,0);
//           inpack.r1[j][i]-=delta;
//           switch(j){
//             case 0:
//                fprintf(mtd_data.fplog,"TESTING: DD_DR1 [ %d ][ %d ]:  X %d ANAL %18.9f NUMER %18.9f DELTA %18.9f\n",l,m,i,dd_dr1[l][m][j][i],(work->d[l][m]- oldd[l][m])/delta,dd_dr1[l][m][j][i]-(work->d[l][m]- oldd[l][m])/delta);break;
//             case 1:
//                fprintf(mtd_data.fplog,"TESTING: DD_DR1 [ %d ][ %d ]:  Y %d ANAL %18.9f NUMER %18.9f DELTA %18.9f\n",l,m,i,dd_dr1[l][m][j][i],(work->d[l][m]- oldd[l][m])/delta,dd_dr1[l][m][j][i]-(work->d[l][m]- oldd[l][m])/delta);break;
//             case 2:
//                fprintf(mtd_data.fplog,"TESTING: DD_DR1 [ %d ][ %d ]:  Z %d ANAL %18.9f NUMER %18.9f DELTA %18.9f\n",l,m,i,dd_dr1[l][m][j][i],(work->d[l][m]- oldd[l][m])/delta,dd_dr1[l][m][j][i]-(work->d[l][m]- oldd[l][m])/delta);break;
//
//
//           }
//       }
//    }
//  }
//}
	exit(0);
};
