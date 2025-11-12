/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2025 The plumed team
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
#ifndef __PLUMED_TEST_MACROS
#define __PLUMED_TEST_MACROS
#include <iostream>
#include <fstream>
#include <algorithm>
#include <iterator>
#include <vector>
#include <string>

///A very barebone tee implementation
///Useful for runnig make test on the fly with
class tee {
  std::ofstream ofs;
public:
  tee(std::string filename) : ofs(filename) {}
  template<typename T>
  tee& operator<<(const T& t) {
    ofs<<t;
    std::cout <<t;
    return *this;
  }
  template<typename Iterable>
  tee&  dump(const Iterable& v) {
    using value_type=typename Iterable::value_type;
    std::copy(v.begin(),v.end(),std::ostream_iterator<value_type>(ofs," "));
    std::copy(v.begin(),v.end(),std::ostream_iterator<value_type>(std::cout," "));
    return *this;
  }
};

//a single letter change to tee for having the output reported in the report.txt, for conveniece
class etee {
  std::ofstream ofs;
public:
  etee(std::string filename) : ofs(filename) {}
  template<typename T>
  etee& operator<<(const T& t) {
    ofs<<t;
    std::cerr <<t;
    return *this;
  }
  template<typename Iterable>
  etee&  dump(const Iterable& v) {
    using value_type=typename Iterable::value_type;
    std::copy(v.begin(),v.end(),std::ostream_iterator<value_type>(ofs," "));
    std::copy(v.begin(),v.end(),std::ostream_iterator<value_type>(std::cerr," "));
    return *this;
  }
};

namespace PLMDTests {

/** A simple structure that encloses a title and a value
*
* to use it store it in a vector (use the shortcut `PLMDTests::testcases`)
* and then you can "navigate" the list of testcase with a rangebased loop
* `for(auto&[title,value] : listOfTestcases){...}`
* in this way you can output the title of the testcase with relative simplicity
*/
template <typename T>
struct namedVal {
  std::string name;
  T val;
};

template <typename T>
using testcases=std::vector<namedVal<T>>;

}
#endif //__PLUMED_TEST_MACROS
