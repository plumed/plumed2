#ifndef __PLUMED_HistogramBead_h
#define __PLUMED_HistogramBead_h

#include <cassert>
#include <cmath>

namespace PLMD {

class HistogramBead{
private:	
	bool init;
	double lowb;
	double highb;
	double width;
public:
	HistogramBead();
	void set(double l, double h, double w);
	double calculate(double x, double&df) const;
	double getlowb() const ;
	double getbigb() const ;
};	

inline
HistogramBead::HistogramBead():
init(false)
{		
} 	

inline
void HistogramBead::set( double l, double h, double w){
	init=true; lowb=l; highb=h; width=w*(h-l); 
	assert( highb>lowb && width>0 );
}

inline
double HistogramBead::getlowb() const { return lowb; }
	
inline
double HistogramBead::getbigb() const { return highb; }
	
inline
double HistogramBead::calculate( double x, double& df ) const {
	const double pi=3.141592653589793238462643383279502884197169399375105820974944592307;
	assert(init); double lowB, upperB, val=0;
	lowB = ( lowb - x ) / ( sqrt(2) * width );
	upperB = ( highb - x ) / ( sqrt(2) * width ) ;
	df = ( exp( -lowB*lowB ) - exp( -upperB*upperB ) ) / ( sqrt(2*pi)*width );
	return 0.5*( erf( upperB ) - erf( lowB ) );
}

}

#endif
