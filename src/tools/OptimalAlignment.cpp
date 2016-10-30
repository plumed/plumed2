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
#include "OptimalAlignment.h"
#include "Kearsley.h"
#include "Log.h"
#include <cmath>
#include <iostream>
#include <cstdlib>
#include "Random.h"

using namespace std;
namespace PLMD{

OptimalAlignment::OptimalAlignment( const  std::vector<double>  & align, const  std::vector<double>  & displace, const std::vector<Vector> & p0, const std::vector<Vector> & p1 , Log* &log )
:log(log){

	// kearsley init to null
	mykearsley=NULL;
	if (mykearsley==NULL) {
		mykearsley=new Kearsley(p0,p1,align,log);
	}
	// copy the structure into place
	assignP0(p0);
	assignP1(p1);
	assignAlign(align);
	assignDisplace(displace);

	// basic check
	if(p0.size() != p1.size()){
		(this->log)->printf("THE SIZE OF THE TWO FRAMES TO BE ALIGNED ARE DIFFERENT\n");
    }
	// fast behaviour: if all the alignment and displacement are 1.0 then go for fast 
	fast=true;
	for (unsigned i=0;i<align.size();i++ ){
                if(align[i]!=displace[i] || align[i]!=1.0)fast=false;
	}

}

OptimalAlignment::~OptimalAlignment(){
	if (mykearsley!=NULL)delete mykearsley;
}

void OptimalAlignment::assignP0(  const std::vector<Vector> & p0 ){
	this->p0=p0;
	if(mykearsley!=NULL){mykearsley->assignP0(p0);}else{cerr<<"kearsley is not initialized"<<endl; exit(0);}
}

void OptimalAlignment::assignP1(  const std::vector<Vector> & p1 ){
	this->p1=p1;
	if(mykearsley!=NULL){mykearsley->assignP1(p1);}else{cerr<<"kearsley is not initialized"<<endl; exit(0);}
}

void OptimalAlignment::assignAlign(  const std::vector<double> & align ){
	this->align=align;
	if(mykearsley!=NULL){mykearsley->assignAlign(align);}else{cerr<<"kearsley is not initialized"<<endl; exit(0);}
}

void OptimalAlignment::assignDisplace(  const std::vector<double> & displace ){
	this->displace=displace;
}


double OptimalAlignment::calculate(bool squared, std::vector<Vector> & derivatives){

	bool rmsd=!squared ;  

	double err;

	// at this point everything should be already in place for calculating the alignment (p1,p2,align)
	// here everything is done with kearsley algorithm. Extension to other optimal alignment algos is
	// possible here below with a switch

	err=mykearsley->calculate(rmsd);  // this calculates the MSD: transform into RMSD

	// check findiff alignment
	 //mykearsley->finiteDifferenceInterface(rmsd);

	if(fast){
		//log->printf("Doing fast: ERR %12.6f \n",err);
		derrdp0=mykearsley->derrdp0;
		derrdp1=mykearsley->derrdp1;
		derivatives=derrdp0;
	}else{
		/// TODO this interface really sucks since is strongly asymmetric should be re-engineered.
		err=weightedAlignment(rmsd);
		//log->printf("Doing slow: ERR %12.6f \n",err);
		derivatives=derrdp0;
	}
	// destroy the kearsley object?

	return err;
}

#ifdef __INTEL_COMPILER
#pragma intel optimization_level 2
#endif
/// this does the weighed alignment if the vector of alignment is different from displacement
double OptimalAlignment::weightedAlignment( bool rmsd){
	double tmp0,tmp1,walign,wdisplace,const1,ret;
	unsigned  i,k,l,m,n,o,oo,mm;

	unsigned natoms=p0.size();

	Kearsley *myk=mykearsley;  /// just a shortcut

	/// TODO : these two blocks can be calculated once forever after the initialization (exception for certain methods?)

	/// clear derivatives
	if (derrdp0.size()!=natoms)derrdp0.resize(natoms);
	if (derrdp1.size()!=natoms)derrdp1.resize(natoms);

	// clear the container
	for(i=0;i<natoms;i++){
		derrdp0[i][0]=derrdp0[i][1]=derrdp0[i][2]=0.;
		derrdp1[i][0]=derrdp1[i][1]=derrdp1[i][2]=0.;
	}

	walign=0.;
	vector<int> alignmap;
	for(i=0;i<natoms;i++){
		if (align[i]>0.){
			alignmap.push_back(i);
			walign+=align[i];
		}
		if (align[i]<0.){cerr<<"FOUND ALIGNMENT WEIGHT NEGATIVE!"<<endl;exit(0);};
	}

	wdisplace=0.;
	vector<int> displacemap;
	for(i=0;i<natoms;i++){
		if (displace[i]>0.){
			displacemap.push_back(i);
			wdisplace+=displace[i];
		}
		if (displace[i]<0.){cerr<<"FOUND ALIGNMENT WEIGHT NEGATIVE!"<<endl;exit(0);};
	}


	tmp0=0.;

	vector<Vector> array_n_3;
	array_n_3.resize(natoms);
	for(i=0;i<array_n_3.size();i++)array_n_3[i][0]=array_n_3[i][1]=array_n_3[i][2]=0.;

	//    err= (1/totdisplace) sum_k_l   displace_k*((p0reset_k_l- sum_k_m rot_l_m*p1reset_k_m )**2)

	//for(kk=0;kk<displacemap.size();kk++){
	//	k=displacemap[kk];
	for(k=0;k<natoms;k++){

		for(l=0;l<3;l++){

			tmp1=0.;
			// contribution from rotated reference frame //
			for(m=0;m<3;m++){
				tmp1-=myk->rotmat0on1[l][m]*myk->p1reset[k][m];
			}

			// contribution from running centered frame //
			tmp1+= myk->p0reset[k][l];

			array_n_3[k][l]=tmp1; // store coefficents for derivative usage//
			tmp0+=tmp1*tmp1*displace[k]; //squared distance added//
		}

	}

	tmp0=tmp0/wdisplace;

  // log->printf(" ERRR NEW %f \n",tmp0);

	ret=tmp0;

	/* DERIVATIVE CALCULATION:respect to running frame */
	for(k=0;k<natoms;k++){
		for(l=0;l<3;l++){

			tmp1 =2.*array_n_3[k][l]*displace[k]/wdisplace ; //ok

			const1=2.*align[k]/(walign*wdisplace);

			if(const1>0.){
				for(oo=0;oo<displacemap.size();oo++){
					o=displacemap[oo];
					tmp1 -=const1*array_n_3[o][l]*displace[o];   //ok
				}
			}

			for(mm=0;mm<displacemap.size();mm++){
				m=displacemap[mm];
				const1=2.* displace[m]/wdisplace ;
				for(n=0;n<3;n++){
					tmp0=0.;
					for(o=0;o<3;o++){
						int ind=n*3*3*natoms+o*3*natoms+l*natoms+k;  //ok
						tmp0+=myk->dmatdp0[ind]*myk->p1reset[m][o];
					}
					tmp0*=-const1*array_n_3[m][n];

					tmp1+=tmp0;
				}
			}

			derrdp0[k][l]=tmp1;

		}
	}
	//exit(0);

	//return ret;
	bool do_frameref_der=true;

	// derivatives of
	//    err= (1/totdisplace) sum_k_l   displace_k*((p0reset_k_l- sum_m rot_l_m*p1reset_k_m )**2)
	// respect p1:
	// derr_dp1=(1/totdisplace) sum_k_l 2*displace_k*(p0reset_k_l- sum_m rot_l_m*p1reset_k_m )
	// 									 *d/dp1 ( p0reset_k_l- sum_m rot_l_m*p1reset_k_m)
	// =
	// (1/totdisplace) sum_k_l 2*displace_k*(p0reset_k_l- sum_m rot_l_m*p1reset_k_m )*
	//							*(d/dp1 p0reset_k_l
	//								- sum_m (d/dp1  rot_l_m)*p1reset_k_m
	//								- sum_m rot_l_m*(d/dp1 p1reset_k_m ) )
	// =
	// 				sum_k_l 2*displace_k/totdisplace* array_n_3_k_l
	//							*(- sum_m (d/dp1  rot_l_m)*p1reset_k_m
	//								- sum_m rot_l_m*(d/dp1 p1reset_k_m ) )

	if(do_frameref_der){
		for(k=0;k<natoms;k++){
//			for(kk=0;kk<displacemap.size();kk++){
//				k=displacemap[kk];

			for(l=0;l<3;l++){

				tmp1=0.;
				for(mm=0;mm<displacemap.size();mm++){
					m=displacemap[mm];
					const1=2.* displace[m]/wdisplace ;
					for(n=0;n<3;n++){
						tmp0=0.;
						for(o=0;o<3;o++){
							int ind=n*3*3*natoms+o*3*natoms+l*natoms+k;
							tmp0+=myk->dmatdp1[ind]*myk->p1reset[m][o];
						}
						tmp0*=-const1*array_n_3[m][n];
						tmp1+= tmp0;
					}

				}

				tmp0=0.;
				for(o=0;o<3;o++){
							tmp0+=array_n_3[k][o]*myk->rotmat0on1[o][l];
				}
				tmp1+=-tmp0*2.*displace[k]/wdisplace;

				tmp0=0.;

				for(mm=0;mm<displacemap.size();mm++){
					m=displacemap[mm];
					for(o=0;o<3;o++){
						tmp0+=array_n_3[m][o]*myk->rotmat0on1[o][l]*displace[m];
					}
				}
				tmp1 += tmp0*2.*align[k]/(walign*wdisplace);

				derrdp1[k][l]=tmp1;

			}
		}
	}

	/// probably not really the way it should be
	if (rmsd){
		ret=sqrt(ret);
			double tmp=0.5/ret;
			for(unsigned i=0;i<natoms;i++){
					derrdp0[i][0]=derrdp0[i][0]*tmp;
					derrdp0[i][1]=derrdp0[i][1]*tmp;
					derrdp0[i][2]=derrdp0[i][2]*tmp;
					derrdp1[i][0]=derrdp1[i][0]*tmp;
					derrdp1[i][1]=derrdp1[i][1]*tmp;
					derrdp1[i][2]=derrdp1[i][2]*tmp;

			}
	}

	return ret;
}

double OptimalAlignment::weightedFindiffTest( bool rmsd){

        Random rnd;

	log->printf("Entering rmsd finite difference test system\n ");
	log->printf("RMSD OR MSD: %s\n",(rmsd)?"rmsd":"msd");
	log->printf("-------------------------------------------\n");
	log->printf("TEST1: derivative of the value (derr_dr0/derr_dr1)\n");
	//// test 1
	double step=1.e-8,olderr,delta,err;
	vector<Vector> fakederivatives;
	fakederivatives.resize(p0.size());
	fast=false;
	// randomizing alignments and  displacements
/*	for (i=0;i<p0.size();i++){
		// draw a random number
	    delta=drand48();
	    delta1=drand48();
	    if(delta>delta1){
	    	align[i]=delta;
	    }else{align[i]=0.;};
	    delta=drand48();
	    delta1=drand48();
	    if(delta>delta1){
	    	displace[i]=delta;
	    }else{displace[i]=0.;}
	}*/
	//// get initial value of the error and derivative of it
	assignAlign(align);
	assignDisplace(displace);
	olderr=calculate(rmsd, fakederivatives);

	log->printf("INITIAL ERROR VALUE: %e\n",olderr);

	// randomizing alignments and  displacements
	log->printf("TESTING: derrdp0 \n");

	for(unsigned j=0;j<3;j++){
	   for(unsigned i=0;i<derrdp0.size();i++){
	       // random displacement
	       delta=(rnd.RandU01()-0.5)*2*step;
	       p0[i][j]+=delta;
	       assignP0( p0 );
	       err=calculate(rmsd, fakederivatives);
	       p0[i][j]-=delta;
	       assignP0( p0 );
	       switch(j){
	         case 0:
	             log->printf("TESTING: X  %4u ANAL %18.9f NUMER %18.9f DELTA %18.9f DISP %6.2f ALIGN %6.2f \n",i,derrdp0[i][j],(err-olderr)/delta,derrdp0[i][j]-(err-olderr)/delta,displace[i],align[i]);break;
	         case 1:
	             log->printf("TESTING: Y  %4u ANAL %18.9f NUMER %18.9f DELTA %18.9f DISP %6.2f ALIGN %6.2f \n",i,derrdp0[i][j],(err-olderr)/delta,derrdp0[i][j]-(err-olderr)/delta,displace[i],align[i]);break;
	         case 2:
	             log->printf("TESTING: Z  %4u ANAL %18.9f NUMER %18.9f DELTA %18.9f DISP %6.2f ALIGN %6.2f \n",i,derrdp0[i][j],(err-olderr)/delta,derrdp0[i][j]-(err-olderr)/delta,displace[i],align[i]);break;

	       }
	   }
	}

	log->printf("TESTING: derrdp1 \n");
	for(unsigned j=0;j<3;j++){
	   for(unsigned i=0;i<derrdp1.size();i++){
	       // random displacement
	       delta=(rnd.RandU01()-0.5)*2*step;
	       p1[i][j]+=delta;
	       assignP1( p1 );
	       err=calculate(rmsd, fakederivatives);
	       p1[i][j]-=delta;
	       assignP1( p1 );
	       switch(j){
	         case 0:
	             log->printf("TESTING: X  %4u ANAL %18.9f NUMER %18.9f DELTA %18.9f DISP %6.2f ALIGN %6.2f \n",i,derrdp1[i][j],(err-olderr)/delta,derrdp1[i][j]-(err-olderr)/delta,displace[i],align[i]);break;
	         case 1:
	             log->printf("TESTING: Y  %4u ANAL %18.9f NUMER %18.9f DELTA %18.9f DISP %6.2f ALIGN %6.2f \n",i,derrdp1[i][j],(err-olderr)/delta,derrdp1[i][j]-(err-olderr)/delta,displace[i],align[i]);break;
	         case 2:
	             log->printf("TESTING: Z  %4u ANAL %18.9f NUMER %18.9f DELTA %18.9f DISP %6.2f ALIGN %6.2f \n",i,derrdp1[i][j],(err-olderr)/delta,derrdp1[i][j]-(err-olderr)/delta,displace[i],align[i]);break;

	       }
	   }
	}
	exit(0);

// This is to avoid warnings:
  return 0.0;

}



}
