#ifndef __PLUMED_WithCmd_h
#define __PLUMED_WithCmd_h

#include <string>

namespace PLMD{

/// Base for classes with cmd() method.
/// This is an abstract base class for classes with
/// cmd() method. It takes care of "const" cast, and
/// in the future it may be used to enforce some sort
/// of type checking on passed arguments.
class WithCmd{
public:
  virtual ~WithCmd();
/// Const val version, which indeed just overrides the const and call the virtual method.
  void cmd(const std::string& key,const void*val);
/// This has to be implemented in daughter classes.
  virtual void cmd(const std::string& key,void*val=NULL)=0;
};

inline
void WithCmd::cmd(const std::string& key,const void*val){
// this is nasty trick:
  cmd(key,const_cast<void*>(val));
// in this manner, a const pointer can be used for val, allowing the user to pass
// arguments such as cmd("pippo","pluto")
// but here we override the const
}

inline
WithCmd::~WithCmd(){
// do nothing
// here just to allow inheriting from this class properly
}

}

#endif
