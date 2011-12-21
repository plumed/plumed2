#ifndef __PLUMED_Units_h
#define __PLUMED_Units_h

namespace PLMD{

/// Small utility class.
/// It has no implemented methods, and all its data are public.
/// It just simplify the syntax of functions which should pass the
/// value of all the units.
class Units{
public:
/// Units for energy, expressed in kj/mol (e.g. 4.184 means kcal/mol)
  double energy;
/// Units for length, expressed in nm (e.g. 0.1 means A)
  double length;
/// Units for time, expressed in ps (e.g. 0.001 means fs)
  double time;
// Constructor, setting default values (1.0)
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
