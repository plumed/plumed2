#ifndef __PLUMED_Kearsley_h
#define __PLUMED_Kearsley_h

#include "Vector.h"
#include "Tensor.h"
#include "Log.h"
#include "Matrix.h"
#include <vector>

namespace PLMD{

/// A class that implements Kearsley's calculation 

class Kearsley
{
  // general log
  Log &log;
  // position of atoms 
  std::vector<Vector> p0;
  std::vector<Vector> p1;
  // alignment and displace vectors
  std::vector<double> align;

  bool com0_is_removed;
  bool com1_is_removed;

public:
  // these are all the possible data that one might have from Kearsley
  // error  
  double err;
  // displacement
  std::vector<Vector> diff0on1;
  std::vector<Vector> diff1on0;

  // center of mass
  Vector com0;
  Vector com1;
  // position resetted wrt coms 
  std::vector<Vector> p0reset;
  std::vector<Vector> p1reset;
  // position rotated
  std::vector<Vector> p0rotated;
  std::vector<Vector> p1rotated;
  // rotation matrices 
  Tensor rotmat0on1,rotmat1on0;
  // derivatives
  std::vector<Vector> derrdp0;
  std::vector<Vector> derrdp1;
  //derivative of the rotation matrix 9 x 3N 
  std::vector<double> dmatdp0; //(3*3*3*natoms)
  std::vector<double> dmatdp1;

  // initialize the structure
  Kearsley(  const std::vector<Vector> &p0, const std::vector<Vector> &p1,  const std::vector<double> &align , Log &log);
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

