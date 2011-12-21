#ifndef __PLUMED_Kearsley_h
#define __PLUMED_Kearsley_h

#include "Vector.h"
#include "Tensor.h"
#include "Log.h"
#include <vector>

using namespace std;

namespace PLMD{

/// A class that implements Kearsley's calculation 

class Kearsley 
{
  // general log
  Log &log;
  // position of atoms 
  vector<Vector> p1;
  vector<Vector> p2;
  // alignment and displace vectors
  vector<double> align;
  // error  
  double err;
  // displacement
  vector<Vector> diff;  
  // center of mass
  Vector com1;
  Vector com2;
  bool com1_is_removed;
  bool com2_is_removed;
  // position resetted wrt coms 
  vector<Vector> p1reset;
  vector<Vector> p2reset;
  // rotation matrices 
  Tensor rotmat1on2,rotmat2on1;
  // derivatives
  vector<Vector> derrdp1;
  vector<Vector> derrdp2;
  //derivative of the rotation matrix 9 x 3N 
  vector< vector<Vector> > dmatdp1;
  vector< vector<Vector> > dmatdp2;

public:
  // initialize the structure
  Kearsley(  const vector<Vector> &p1, const vector<Vector> &p2,  const vector<double> &align , Log &log);
  // switch the assignment of the structures
  void assignP1(const std::vector<Vector> & p1);
  void assignP2(const std::vector<Vector> & p2);
  // this makes the real buisness
  double calculate( );
};

}

#endif

