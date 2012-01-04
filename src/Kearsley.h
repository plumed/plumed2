#ifndef __PLUMED_Kearsley_h
#define __PLUMED_Kearsley_h

#include "Vector.h"
#include "Tensor.h"
#include "Log.h"
#include "Matrix.h"
#include <vector>

namespace PLMD{

/// A class that implements Kearsley's calculation 
/// which is optimal alignment via quaternion and
/// analytical derivatives via perturbation theory

class Kearsley
{
  /// general log reference that needs to be initialized when constructed
  Log &log;
  /// position of atoms (first frame. In md is the running frame)
  std::vector<Vector> p0;
  /// position of atoms (second frame. In md is the  reference frame)
  std::vector<Vector> p1;
  /// alignment weight: the rmsd/msd that it provides is only based on this scalar
  std::vector<double> align;

  bool com0_is_removed;
  bool com1_is_removed;

public:
  /// error: the distance between two frames (might be rmsd/msd. See below)
  double err;
  /// displacement: the vector that goes from the p0 onto p1
  std::vector<Vector> diff0on1;
  /// displacement: the vector that goes from the p1 onto p0 (via inverse rotation)
  std::vector<Vector> diff1on0;

  /// center of mass of p0
  Vector com0;
  /// center of mass of p1
  Vector com1;
  /// position resetted wrt coms p0
  std::vector<Vector> p0reset;
  /// position resetted wrt coms p1
  std::vector<Vector> p1reset;
  /// position rotated: p0
  std::vector<Vector> p0rotated;
  /// position rotated: p1
  std::vector<Vector> p1rotated;
  /// rotation matrices p0 on p1 and reverse (p1 over p0)
  Tensor rotmat0on1,rotmat1on0;
  /// derivatives: derivative of the error respect p0
  std::vector<Vector> derrdp0;
  /// derivatives: derivative of the error respect p1
  std::vector<Vector> derrdp1;
  /// derivative of the rotation matrix
  /// note the dimension 3x3 x 3 x N
  std::vector<double> dmatdp0;
  std::vector<double> dmatdp1;

  /// constructor: need the two structure, the alignment vector and  the log reference
  Kearsley(  const std::vector<Vector> &p0, const std::vector<Vector> &p1,  const std::vector<double> &align , Log &log);
  /// switch the assignment of the structure p0 (e.g. at each md step)
  void assignP0(const std::vector<Vector> & p0);
  /// derivatives: derivative of the error respect p1
  void assignP1(const std::vector<Vector> & p1);
  /// transfer the alignment vector
  void assignAlign(const std::vector<double> & align);
  /// finite differences of all the relevant quantities: takes a bool which decides if giving back rmsd or not (msd in this case)
  void finiteDifferenceInterface(bool rmsd);
  // this makes the real calculation: the rmsd bool decides wether doing rmsd or msd
  double calculate( bool rmsd );
};

}

#endif

