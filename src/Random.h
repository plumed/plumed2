#ifndef __RANDOM_H
#define __RANDOM_H

#include <iostream>

namespace PLMD{

class Random{
	static const int IA=16807,IM=2147483647,IQ=127773,IR=2836,NTAB=32;
	static const int NDIV=(1+(IM-1)/NTAB);
	static const double EPS;
	static const double AM;
	static const double RNMX;
	static const double fact;
	bool incPrec;
        bool switchGaussian;
        double saveGaussian;
	int iy;
	int iv[NTAB];
	int idum;
	std::string name;
public:
	Random(const std::string & name="");
	void setSeed(int idum);
	double RandU01();
	double U01();
	double U01d();
	void WriteStateFull(std::ostream &)const;
	void ReadStateFull (std::istream &);
	void fromString(const std::string & str);
	void toString(std::string & str)const;
	friend std::ostream & operator<<(std::ostream & out,const Random & rng){
		rng.WriteStateFull(out); return out;
	};
	friend std::istream & operator>>(std::istream & in,Random & rng){
		rng.ReadStateFull(in); return in;
	};
	double Gaussian();
	void IncreasedPrecis(bool i){incPrec=i;};
};

}

#endif

