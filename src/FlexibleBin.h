#ifndef __PLUMED_FlexibleBin_h
#define __PLUMED_FlexibleBin_h

#include<vector>
using namespace std;

namespace PLMD{

class ActionWithArguments;

class FlexibleBin{
	private:
		double sigma;	
		const int type;
		// this contains all the infos about the CVs including periodicity
		ActionWithArguments *paction;
		// variance is the matrix that really matters
		vector<double> variance;	
		// this is only there
		vector<double> average;
	public:
		/// a constructor that takes the pointer of the action that contains it
		FlexibleBin(int type,ActionWithArguments *paction, double const &d);
		~FlexibleBin();
		/// update the average (always for diffusion) or calculate the geom covariance (  only when do_when_zero is zero)
		void update(bool nowAddAHill );
		vector<double> getMatrix() const;
		vector<double> getInverseMatrix() const;
		enum AdaptiveHillsType { none, diffusion, geometry }; 
};



}
#endif
