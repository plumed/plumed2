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


}

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

	err=mykearsley->calculate(rmsd);  // this calculates the MSD: transform into RMSD

	// check findiff alignment
	// mykearsley->finiteDifferenceInterface(rmsd);

	if(fast){
		log.printf("Doing fast: ERR %12.6f \n",err);
		derrdp0=mykearsley->derrdp0;
		derrdp1=mykearsley->derrdp1;
		derivatives=derrdp0;
	}else{
		/// TODO this interface really sucks since is strongly asymmetric should be re-engineered.
		err=weightedAlignment(rmsd);
		log.printf("Doing slow: ERR %12.6f \n",err);
		derivatives=derrdp0;
	}
	// destroy the kearsley object?

	return err;
};
/// this does the weighed alignment if the vector of alignment is different from displacement
double OptimalAlignment::weightedAlignment( bool rmsd){
	double tmp0,tmp1,walign,wdisplace,ndisplace,const1,ret;
	unsigned  i,j,k,l,m,n,o,oo,mm,kk;

	unsigned natoms=p0.size();

	Kearsley *myk=mykearsley;  /// just a shortcut

	/// TODO : these two blocks can be calculated once forever after the initialization (exception for certain methods?)

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

	for(kk=0;kk<displacemap.size();kk++){
		k=displacemap[kk];
		for(l=0;l<3;l++){

			tmp1=0.;
			// contribution from rotated reference frame //
			for(m=0;m<3;m++){
				tmp1-=myk->rotmat0on1[l][m]*myk->p1reset[k][m];
			}

			// contribution from running centered frame //
			tmp1+= myk->p0reset[k][l];
	//		log.printf("WEIGHTED %3d COMP %1d VAL %f  XX %f CUMUL %f \n",k,l,tmp1*sqrt(displace[k]/wdisplace),tmp1, tmp0);

			array_n_3[k][l]=tmp1; // store coefficents for derivative usage//
			tmp0+=tmp1*tmp1*displace[k]; //squared distance added//
		}
	}
  //  log.printf(" ERRR NEW1 %f %f \n",tmp0,wdisplace);

	tmp0=tmp0/wdisplace;

   // log.printf(" ERRR NEW %f \n",tmp0);

	ret=tmp0;

	/* DERIVATIVE CALCULATION:respect to running frame */
	for(k=0;k<natoms;k++){
		for(l=0;l<3;l++){

			tmp1 =2.*array_n_3[k][l]*displace[k]/wdisplace ;

			const1=2.*align[k]/(walign*wdisplace);

			if(const1!=0.){
				for(oo=0;oo<displacemap.size();oo++){
					o=displacemap[oo];
					tmp1 -=const1*array_n_3[o][l]*displace[o];
				}
			}

				for(mm=0;mm<displacemap.size();mm++){
					m=displacemap[mm];
					const1=2.* displacemap[m]/wdisplace ;
					for(n=0;n<3;n++){
						tmp0=0.;
						for(o=0;o<3;o++){
							int ind=n*3*3*natoms+o*3*natoms+l*natoms+k;
							tmp0+=myk->dmatdp0[ind]*myk->p1reset[m][o];
						}
						tmp0*=-const1*array_n_3[m][n];

						tmp1+=tmp0;
					}
				}

			derrdp1[k][l]=tmp1;

		}
	}


	bool do_frameref_der;

	if(do_frameref_der){
		for(k=0;k<natoms;k++){
			for(l=0;l<3;l++){
				/////////////////////////////////////
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

				/////////////////////////////////////
			}
		}
	}

	/// probably not really the way it should be
	if (rmsd){
		ret=sqrt(ret);
			double tmp=0.5/ret;
			for(int ii=0;ii<alignmap.size();ii++){
					i=alignmap[ii];
					derrdp0[i][0]=derrdp0[i][0]*tmp;
					derrdp0[i][1]=derrdp0[i][1]*tmp;
					derrdp0[i][2]=derrdp0[i][2]*tmp;
					derrdp1[i][0]=derrdp1[i][0]*tmp;
					derrdp1[i][1]=derrdp1[i][1]*tmp;
					derrdp1[i][2]=derrdp1[i][2]*tmp;

			}
	}

	//log.printf(" ERRR NEW %f \n",ret);
	return ret;
};

double OptimalAlignment::weightedFindiffTest( bool rmsd){

	log.printf("Entering rmsd finite difference test system\n");
	log.printf("-------------------------------------------\n");
	log.printf("TEST1: derivative of the value (derr_dr0/derr_dr1)\n");
	//// test 1
	unsigned i,j,l,m;
	double step=1.e-7,olderr,delta,delta1,err;
	vector<Vector> fakederivatives;
	fakederivatives.resize(p0.size());
	fast=false;
	// randomizing alignments and  displacements
	for (i=0;i<p0.size();i++){
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
	}
	//// get initial value of the error and derivative of it
	olderr=calculate(fakederivatives);
	log.printf("INITIAL ERROR VALUE: %e\n",olderr);

	// randomizing alignments and  displacements
	log.printf("TESTING: derrdp1 \n");
	for(unsigned j=0;j<3;j++){
	   for(unsigned i=0;i<derrdp1.size();i++){
	       // random displacement
	       delta=(drand48()-0.5)*2*step;
	       p1[i][j]+=delta;
	       assignP1( p1 );
	       err=calculate(fakederivatives);
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
	exit(0);
}



