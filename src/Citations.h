#ifndef __PLUMED_Citations_h
#define __PLUMED_Citations_h

#include <vector>
#include <string>
#include <iosfwd>

namespace PLMD{

/**
Class taking care of bibliography.

This class contains a vector of citations. To add a new citations, use cite(). To print
the entire bibliography, just dump on a ostream. Everytime cite is used, a string
containing the number of the citation is returned. If the same citation is added twice,
the same string is returned, so that this example will produce only two bibliographic items:
\verbatim
#include "Citations.h"
#include <iostream>
int main(int argc,char**argv){
  PLMD::Citations citations;
  std::cout << citations.cite("Pinco e Pallino, Il Piccolo 33, 444 (2012)") << "\n";
  std::cout << citations.cite("Other cite") << "\n";
  std::cout << citations.cite("Pinco e Pallino, Il Piccolo 33, 444 (2012)") << "\n";

  std::cout << "Bibliography\n"<< citations;
  return 0;
}
\endverbatim
*/

class Citations{
  std::vector<std::string> items;
public:
/// Add a citation.
/// It returns a string containing the reference number, something like "[10]"
  std::string cite(const std::string &);
/// Dumps the bibliography.
/// It writes on the ostream the list of all the bibliographic items
/// prefixed with their reference number
  friend std::ostream &operator<<(std::ostream &,const Citations&);
};

std::ostream & operator<<(std::ostream &,const Citations&);

}

#endif
