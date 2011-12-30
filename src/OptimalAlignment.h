#ifndef __PLUMED_OptimalAlignment_h
#define __PLUMED_OptimalAlignment_h

#include "Vector.h"
#include "Tensor.h"
#include "Kearsley.h"
#include "Log.h"
#include <vector>

using namespace std;

namespace PLMD{


/// A class that implements Kearsley's calculation 

class	OptimalAlignment 
{
  // kearsley should contain all the derivatives respect the rotation matrix and respect to the alignment 
  Kearsley *mykearsley;	
  // the optimal alignment should only contain the derivative respect the distance 
   std::vector<double> displace;
   std::vector<double> align;
   std::vector<Vector> p0;
   std::vector<Vector> p1;
   std::vector<Vector> derrdp0;
   std::vector<Vector> derrdp1;

   Log &log;
   bool fast;	

public:
  OptimalAlignment( const  std::vector<double>  & align,  const std::vector<double>   & displace, const std::vector<Vector> & p0, const std::vector<Vector> & p1 , Log &log );
  void assignP0(  const std::vector<Vector> & p0 );
  void assignP1(  const std::vector<Vector> & p1 );
  double calculate( std::vector<Vector> & derivatives);
};

}

#endif

