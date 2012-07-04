#ifndef __PLUMED_Citations_h
#define __PLUMED_Citations_h

#include <vector>
#include <string>
#include <iosfwd>

namespace PLMD{

/**
Class taking care of bibliography.
*/

class Citations{
  std::vector<std::string> items;
public:
/// Add a citation.
/// It returns a string containing the reference number, something like "[10]"
  std::string cite(const std::string &);
/// Dumps the bibliography.
/// It writes on the ostream the list of all the bibliographic itemrs
/// prefixed with their reference number
  friend std::ostream &operator<<(std::ostream &,const Citations&);
};

std::ostream & operator<<(std::ostream &,const Citations&);

}

#endif
