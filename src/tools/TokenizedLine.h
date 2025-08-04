
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
#ifndef __PLUMED_tools_TokenizedLine_h
#define __PLUMED_tools_TokenizedLine_h
#include "Tools.h"
#include <map>
#include <vector>
namespace PLMD {
///This class abstracts the input line in tokens.
///
///The underlying container accelerates the lookup for the keywords
///ISSUE: In case of a vector keyword and the use of the '@replicas:' idiom:
/// If the user wants to chante only the FIRS element of the vector with @replica
/// the parse will interpret the whole string as a replica string
/// There are no problem with @replicas: applied on other components:
///  - {@replicas:{1,3,4},1,3} won't work (the result ideally is equivalent to "@replicas:{{1,1,3},{3,1,3},{4,1,3}}"
///  - {1,@replicas:{1,3,4},3} will work (the result is equivalent to "@replicas:{{1,1,3},{1,3,3},{1,4,3}}"
class TokenizedLine {
public:
  // the first element is the whole line, if contains @replica: will get also divided
  using line=std::vector<std::string>;
  using mapType=std::map<std::string,line,std::less<void>>;
private:
  using vectorIt = typename std::vector<std::string>::iterator;
  using const_vectorIt = typename std::vector<std::string>::const_iterator;
  mapType tokens;
public:
  static std::string_view replica(const line&,int rep=-1);
  struct presentAndFound {
    bool present;
    bool found;
  };
/// Initializer from vector iterators
  TokenizedLine(vectorIt begin, vectorIt end);
/// Initializer from vector iterators
  TokenizedLine(const_vectorIt begin, const_vectorIt end);
/// Initializer from a vector of strings (ideally the output of getWords)
  TokenizedLine(const std::vector<std::string>&);
///return a plain string with the all the current KEY=value combinations, it is possible to clear the tokens after that
  std::string convertToString(bool alsoClear);
///returns the list of the keys:
  std::string keyList(std::string_view sep = ", ");
///return a keyword and its argument
  std::string getKeyword(std::string_view key) const;
///return the value of the asked key, "" otherwise;
  std::string getValue(std::string_view key, int rep=-1) const;
/// Return the size of the underlying container;
  std::size_t size() const;
/// Returns true if the underlying container is empty
  bool empty() const;
/// Read a value from the tokens and remove it from the list
  template<typename T>
  presentAndFound readAndRemove(std::string_view key,
                                T& value,
                                int rep=-1);
/// Read a list of values from the tokens and remove it from the list
  template<typename T>
  presentAndFound readAndRemoveVector(std::string_view key,
                                      std::vector<T>& value,
                                      int rep=-1);
///return true if the flag is present and removes it from the tokens
  bool readAndRemoveFlag(std::string_view key);
};



template<typename T>
TokenizedLine::presentAndFound TokenizedLine::readAndRemove(std::string_view key,
    T& value,
    const int replica_index ) {
  auto keytext = tokens.find(key);
  bool present = keytext != tokens.end();
  bool found=false;
  if(present) {
    found = Tools::parse(TokenizedLine::replica(keytext->second,replica_index),
                         value);
    tokens.erase(keytext);
  }
  return {present, found};
}


template<typename T>
TokenizedLine::presentAndFound TokenizedLine::readAndRemoveVector(std::string_view key,
    std::vector<T>& value,
    const int replica_index ) {

  auto keytext = tokens.find(key);
  bool present = keytext != tokens.end();
  bool found=false;
  if(present) {
    found = Tools::parseVector(TokenizedLine::replica(keytext->second,replica_index),
                               value,
                               replica_index);
    tokens.erase(keytext);
  }
  return {present, found};
}


} //namespace PLMD
#endif //__PLUMED_tools_TokenizedLine_h

