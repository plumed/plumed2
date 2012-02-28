#include "FlexibleBin.h"
#include "ActionWithArguments.h"
#include <cmath>
#include <iostream>
#include <vector>

using namespace std;
using namespace PLMD;


FlexibleBin::FlexibleBin(int type, ActionWithArguments *paction ):type(type),paction(paction){
	// initialize the averages and the variance matrices
	if(type==diffusion){
		unsigned ncv=paction->getNumberOfArguments();	
		vector<double> average(ncv*(ncv+1)/2);
		vector<double> variance(ncv*(ncv+1)/2);
	}

}

void FlexibleBin::update(double const sigma){
	unsigned ncv=paction->getNumberOfArguments();
	unsigned dimension=ncv*(ncv+1)/2;	
	// this is done all the times from scratch. It is not an accumulator 
	vector<double> matrix(dimension);
	int k=0;
	// here update the flexible bin according to the needs	
	switch (type){
		case diffusion: 
			// This should be called every time
			cerr<< "Doing diffusion "<<endl; 
			plumed_massert(average.size()==dimension,"Attempting to change the diffusion matrix on the fly. Too bad! " );
			plumed_massert(variance.size()==dimension,"Attempting to change the diffusion matrix on the fly. Too bad! " );
			break;
		case geometry: 
			cerr<< "Doing geometry "<<endl; 
			// now the signal for retrieving the gradients should be already given by checkNeedsGradients.
			// here just do the projections
			/// TODO: call this only when needed 
			for(unsigned i=0;i<ncv;i++){
   			       for(unsigned j=0;j<ncv;j++){
			      	      //fprintf(fp,fmt.c_str(),getProjection(i,j));
			      	      matrix[k]=paction->getProjection(i,j);
			              cerr<<"MM "<<k<<" "<<matrix[k]<<endl;		
			              k++;		
  			       }
			};				
			break;
		default:
			cerr<< "This flexible bin is not recognized  "<<endl; 
			exit(1)	;	
	}

}
