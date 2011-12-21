#ifndef __PLUMED_Angle_h
#define __PLUMED_Angle_h

#include "Vector.h"

namespace PLMD{

/// Class to compute angles.
/// I define it as a class even if it does not contain anything. The reason
/// is that in the future I would like to extend it to contain options about
/// how the calculation should be done. So, for now use it as
/// Angle a;
/// double angle=a.compute(v1,v2);
/// I know it is a bit misleading. If we really do not need to store "options"
/// inside the Angle class, we can remove it later and write compute as
/// a static function.
class Angle{
// still empty, but may accomodate some options in the future
public:
/// Compute the angle between vectors v1 and v2
  double compute(const Vector& v1,const Vector& v2)const;
/// Compute the angle between vectors v1 and v2 and its derivatives wrt v1 and v2
  double compute(const Vector& v1,const Vector& v2,Vector& d1,Vector& d2)const;
};

}

#endif
