#ifndef __PLUMED_CLToolRegister_h
#define __PLUMED_CLToolRegister_h

#include <string>
#include <map>
#include <set>
#include <iostream>
#include <vector>

namespace PLMD{

class CLTool;
class CLToolOptions;

/// Same as ActionRegister, but for CLTools
class CLToolRegister{
/// Write on a stream the list of registered directives
  friend std::ostream &operator<<(std::ostream &,const CLToolRegister&);
/// Pointer to a function which, given the options, create an CLTool
  typedef CLTool*(*creator_pointer)(const CLToolOptions&);
/// Map cltool to a function which creates the related object
  std::map<std::string,creator_pointer> m;
/// Iterator over the map
  typedef std::map<std::string,creator_pointer>::iterator mIterator;
/// Iterator over the map
  typedef std::map<std::string,creator_pointer>::const_iterator const_mIterator;
/// Set of disabled cltools (which were registered more than once)
  std::set<std::string> disabled;
public:
/// Register a new class.
/// \param key The name of the directive to be used in the input file
/// \param key A pointer to a function which creates an object of that class
  void add(std::string key,creator_pointer);
/// Verify if a directive is present in the register
  bool check(std::string cltool);
/// Create an CLTool of the type indicated in the options
/// \param ao object containing information for initialization, such as the full input line, a pointer to PlumedMain, etc
  CLTool* create(const CLToolOptions&ao);
  void remove(creator_pointer);
  ~CLToolRegister();
/// Returns a list of the allowed CLTools
  std::vector<std::string> list()const;
};

/// Function returning a reference to the CLToolRegister.
/// \relates CLToolRegister
/// To avoid problems with order of initialization, this function contains
/// a static CLToolRegister which is built the first time the function is called.
/// In this manner, it is always initialized before it's used
CLToolRegister& cltoolRegister();

std::ostream & operator<<(std::ostream &log,const CLToolRegister&ar);

}


/// Shortcut for CLTool registration
/// \relates PLMD::CLToolRegister
/// For easier registration, this file also provides a macro PLUMED_REGISTER_CLTOOL.
/// \param classname the name of the class to be registered
/// \param directive a string containing the corresponding directive
/// This macro should be used in the .cpp file of the corresponding class
#define PLUMED_REGISTER_CLTOOL(classname,directive) \
  static class classname##RegisterMe{ \
    static PLMD::CLTool* create(const PLMD::CLToolOptions&ao){(void)ao;return new classname;} \
  public: \
    classname##RegisterMe(){PLMD::cltoolRegister().add(directive,create);}; \
    ~classname##RegisterMe(){PLMD::cltoolRegister().remove(create);}; \
  } classname##RegisterMeObject;


#endif

