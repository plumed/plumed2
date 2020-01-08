/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2020 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed.org for more information.

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#ifndef __PLUMED_tools_Citations_h
#define __PLUMED_tools_Citations_h

#include <vector>
#include <string>
#include <iosfwd>

namespace PLMD {

/**
\ingroup TOOLBOX
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

class Citations {
  std::vector<std::string> items;
public:
/// Add a citation.
/// It returns a string containing the reference number, something like "[10]"
  std::string cite(const std::string &);
/// Dumps the bibliography.
/// It writes on the ostream the list of all the bibliographic items
/// prefixed with their reference number
  friend std::ostream &operator<<(std::ostream &,const Citations&);
/// Delete all references
  void clear();
/// Check if bibliography is empty
  bool empty()const;
};

std::ostream & operator<<(std::ostream &,const Citations&);

}

#endif
