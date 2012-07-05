#ifndef __PLUMED_HistogramBead_h
#define __PLUMED_HistogramBead_h

#include <cassert>
#include <cmath>
#include <string>
#include <vector>

namespace PLMD {

class Log;

/**
\ingroup TOOLBOX
A class for calculating whether or not values are within a given range using : \f$ \sum_i \int_a^b G( s_i, \sigma*(b-a) ) \f$    
*/      

class HistogramBead{
private:	
	bool init;
	double lowb;
	double highb;
	double width;
public:
        static std::string documentation( bool dir );
        static std::string histodocs();
        static void generateBins( const std::string& params, const std::string& dd, std::vector<std::string>& bins );  
	HistogramBead();
        std::string description() const ;
        bool hasBeenSet() const;
        void set(const std::string& params, const std::string& dd, std::string& errormsg);
	void set(double l, double h, double w);
	double calculate(double x, double&df) const;
	double getlowb() const ;
	double getbigb() const ;
	void printKeywords(Log& log) const;
};	

inline
HistogramBead::HistogramBead():
init(false)
{		
}

inline
bool HistogramBead::hasBeenSet() const {
  return init;
}

inline
double HistogramBead::getlowb() const { return lowb; }
	
inline
double HistogramBead::getbigb() const { return highb; }
	
inline
double HistogramBead::calculate( double x, double& df ) const {
	const double pi=3.141592653589793238462643383279502884197169399375105820974944592307;
	assert(init); double lowB, upperB;
	lowB = ( lowb - x ) / ( sqrt(2.0) * width );
	upperB = ( highb - x ) / ( sqrt(2.0) * width ) ;
	df = ( exp( -lowB*lowB ) - exp( -upperB*upperB ) ) / ( sqrt(2*pi)*width );
	return 0.5*( erf( upperB ) - erf( lowB ) );
}

}

#endif
