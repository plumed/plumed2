#ifndef __PLUMED_Torsion_h
#define __PLUMED_Torsion_h

#include "Vector.h"
#include "Tools.h"

namespace PLMD{

/// Class to compute torsional angles.
/// I define it as a class even if it does not contain anything. The reason
/// is that in the future I would like to extend it to contain options about
/// how the calculation should be done. So, for now use it as
/// Torsion t;
/// double angle=t.compute(v1,v2,v3);
/// I know it is a bit misleading. If we really do not need to store "options"
/// inside the Torsion class, we can remove it later and write compute as
/// a static function.
class Torsion{
// still empty, but may accomodate some options in the future
public:
/// Compute the angle between the porjections of v1 and v3 on the plane orthogonal
/// to v2. To have a "normal" definition (= plumed1), use it as
/// compute(r01,r12,r23);
/// See ColvarTorsion for a practical usage...
  double compute(const Vector& v1,const Vector& v2,const Vector& v3)const;
/// This should be the version which also computes the derivatives wrt the arguments.
/// Still it is not implemented!! So, you can only do numerical derivatives
///  double compute(const Vector& v1,const Vector& v2,const Vector& v3,Vector& d1,Vector& d2,Vector& d3)const;
};

}

#endif
