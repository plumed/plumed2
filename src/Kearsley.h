#ifndef __PLUMED_Kearsley_h
#define __PLUMED_Kearsley_h

#include "Vector.h"
#include "Tensor.h"
#include "Log.h"
#include "Matrix.h"
#include <vector>

using namespace std;

namespace PLMD{

/// A class that implements Kearsley's calculation 

class Kearsley
{
  // general log
  Log &log;
  // position of atoms 
  vector<Vector> p0;
  vector<Vector> p1;
  // alignment and displace vectors
  vector<double> align;

  bool com0_is_removed;
  bool com1_is_removed;

public:
  // these are all the possible data that one might have from Kearsley
  // error  
  double err;
  // displacement
  vector<Vector> diff0on1;
  vector<Vector> diff1on0;

  // center of mass
  Vector com0;
  Vector com1;
  // position resetted wrt coms 
  vector<Vector> p0reset;
  vector<Vector> p1reset;
  // position rotated
  vector<Vector> p0rotated;
  vector<Vector> p1rotated;
  // rotation matrices 
  Tensor rotmat0on1,rotmat1on0;
  // derivatives
  vector<Vector> derrdp0;
  vector<Vector> derrdp1;
  //derivative of the rotation matrix 9 x 3N 
  vector<double> dmatdp0; //(3*3*3*natoms)
  vector<double> dmatdp1;

  // initialize the structure
  Kearsley(  const vector<Vector> &p0, const vector<Vector> &p1,  const vector<double> &align , Log &log);
  // switch the assignment of the structures
  void assignP0(const std::vector<Vector> & p0);
  void assignP1(const std::vector<Vector> & p1);
  // just for regularity test
  void finiteDifferenceInterface(bool rmsd);
  // this makes the real buisness
  double calculate( bool rmsd );
};

}

#endif

