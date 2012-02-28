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
		// for diffusion based metrics: do allocate the accumulators 
		vector<double> average;
		vector<double> variance;	
		ActionWithArguments *paction;
	public:
		FlexibleBin(int type,ActionWithArguments *paction);
		~FlexibleBin();
		void update(double const s);
		enum AdaptiveHillsType { none, diffusion, geometry }; 
};



}
#endif
