#include "FlexibleBin.h"
#include "ActionWithArguments.h"
#include <cmath>
#include <iostream>
#include <vector>
#include <Matrix.h>

using namespace std;
using namespace PLMD;


FlexibleBin::FlexibleBin(int type, ActionWithArguments *paction,  double const &d ):type(type),paction(paction),sigma(d){
	// initialize the averages and the variance matrices
	if(type==diffusion){
		unsigned ncv=paction->getNumberOfArguments();	
		vector<double> average(ncv*(ncv+1)/2);
		vector<double> variance(ncv*(ncv+1)/2);
	}

}
/// Update the flexible bin 
/// in case of diffusion based: update at every step
/// in case of gradient based: update only when you add the hill 
void FlexibleBin::update(bool nowAddAHill){
	unsigned ncv=paction->getNumberOfArguments();
	unsigned dimension=ncv*(ncv+1)/2;	
	// this is done all the times from scratch. It is not an accumulator 
	unsigned  k=0;
	unsigned i,j;
	vector<double> cv;
	vector<double> delta; 
	double decay=paction->getTimeStep()/sigma;
	// here update the flexible bin according to the needs	
	switch (type){
		// This should be called every time
		case diffusion: 
			cerr<< "Doing diffusion "<<endl; 
                        //
                        // THE AVERAGE VALUE
                        //
			// beware: the pbc
			delta.resize(ncv);      
			for(i=0;i<ncv;i++)cv.push_back(paction->getArgument(i));
			if(average.size()==0){ // initial time: just set the initial vector
				average.resize(ncv);
				for(i=0;i<ncv;i++)average[i]=cv[i]; 	
			}else{ // accumulate  
				for(i=0;i<ncv;i++){
					delta[i]=paction->difference(i,average[i],cv[i]); 
					average[i]+=decay*delta[i]; 
					average[i]=paction->bringBackInPbc(i,average[i]); // equation 8 of "Metadynamics with adaptive Gaussians"
				}
					
			}
                        //
                        // THE VARIANCE
                        //
			if(variance.size()==0){ 
				variance.resize(dimension,0.); // nonredundant members dimension=ncv*(ncv+1)/2;
			}else{
				k=0;
				for(i=0;i<ncv;i++){
					for(j=i;j<ncv;j++){ // upper diagonal loop 
						variance[k]+=decay*(delta[i]*delta[j]-variance[k]);	
			              		//cerr<<"MM "<<"I "<<i<<" J "<<j<<" K "<<k<<" "<<variance[k]<<endl;		
						k++;
					}
				}
			}
			break;
		case geometry: 
                        //
                        //this calculates in variance the \nabla CV_i \dot \nabla CV_j
                        //
			variance.resize(dimension);
			//cerr<< "Doing geometry "<<endl; 
			// now the signal for retrieving the gradients should be already given by checkNeedsGradients.
			// here just do the projections
			// note that the call  checkNeedsGradients() in BiasMetaD takes care of switching on the call to gradients 
			if (nowAddAHill){// geometry is sync with hill deposition 
			        //cerr<< "add a hill "<<endl; 
				k=0;
				for(unsigned i=0;i<ncv;i++){
   				       for(unsigned j=i;j<ncv;j++){
 				              // eq 12 of "Metadynamics with adaptive Gaussians" 		
				      	      variance[k]=sigma*sigma*(paction->getProjection(i,j));
				              k++;		
  				       }
				};				
			};
			break;
		default:
			cerr<< "This flexible bin is not recognized  "<<endl; 
			exit(1)	;	
	}

}

vector<double> FlexibleBin::getMatrix() const{
	return variance;
};

///
/// Calculate the matrix of  (dcv_i/dx)*(dcv_j/dx)^-1 
/// that is needed for the metrics in metadynamics
///
///
vector<double> FlexibleBin::getInverseMatrix() const{
	unsigned ncv=paction->getNumberOfArguments();
	unsigned dimension=ncv*(ncv+1)/2;	
	Matrix<double> matrix(ncv,ncv); 
	unsigned i,j,k;	
	k=0;
	//paction->log<<"------------ GET INVERSE MATRIX ---------------\n";
        // place the matrix in a complete matrix for compatibility 
	for (i=0;i<ncv;i++){
		for (j=i;j<ncv;j++){
			matrix(j,i)=matrix(i,j)=variance[k];
			k++;	
		}
	}
	//paction->log<<"MATRIX\n";
	//matrixOut(paction->log,matrix);	
        // get the inverted matrix
	Matrix<double> invmatrix(ncv,ncv);	
	Invert(matrix,invmatrix);
        vector<double> uppervec(ncv*(ncv+1)/2); 
        // upper diagonal of the inverted matrix (that is symmetric) 
	k=0;
	for (i=0;i<ncv;i++){
		for (j=i;j<ncv;j++){
			uppervec[k]=invmatrix(i,j);
			//paction->log<<"VV "<<i<<" "<<j<<" "<<uppervec[k]<<"\n";	
			k++;
		}
	}
	//paction->log<<"------------ END GET INVERSE MATRIX ---------------\n";

	return uppervec;  
};
