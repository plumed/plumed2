#ifndef __PLUMED_Units_h
#define __PLUMED_Units_h

namespace PLMD{

class Units{
public:
  double energy,length,time;
  Units();
};

inline
Units::Units():
  energy(1.0),
  length(1.0),
  time(1.0)
{ 
}

}

#endif
