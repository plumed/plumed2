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
#ifndef __PLUMED_core_CLTool_h
#define __PLUMED_core_CLTool_h
#include <cstdio>
#include <vector>
#include <string>
#include <cstdio>
#include "tools/Tools.h"
#include "tools/Keywords.h"

namespace PLMD {

class Communicator;

class CLToolOptions {
  friend class CLTool;
  friend class CLToolRegister;
private:
  std::vector<std::string> line;
/// The documentation for this command line tool
  const Keywords& keys;
  static Keywords emptyKeys;
public:
  explicit CLToolOptions(const std::string &name);
  CLToolOptions(const CLToolOptions& co, const Keywords& k);
};

/**
\ingroup INHERIT
This is the abstract base class to use for implementing new command line tool, within it there
is \ref AddingACLTool "information" as to how to go about implemneting a new tool.
*/

class CLTool {
private:
/// The name of this command line tool
  const std::string name;
/// The list of keywords for this CLTool
  const Keywords& keywords;
/// The data read in from the command line stored in a map with the keywords
  std::map<std::string,std::string> inputData;
/// Read the arguments from the command line
  bool readCommandLineArgs( int argc, char**argv, FILE*out );
/// Read the arguments from an input file specified on the command line
  bool readInputFile( int argc, char**argv, FILE* in, FILE*out );
/// Set arguments from the default options provided to Keywords
  void setRemainingToDefault(FILE* out);
public:
/// Set the input data:
  void setInputData(const std::map<std::string,std::string>&inputData) {
    this->inputData=inputData;
  }
  const std::map<std::string,std::string>&getInputData() {
    return this->inputData;
  }
protected:
/// Get the value of one of the command line arguments
  template<class T>
  bool parse(const std::string&key,T&t);
/// Find out whether one of the command line flags is present or not
  void parseFlag(const std::string&key,bool&t);
/// Crash the command line tool with an error
  void error(const std::string& msg);
  template<class T>
  bool parseVector(const std::string&key,std::vector<T>&t);
public:
/// How is the input specified on the command line or in an input file
  enum {unset,commandline,ifile} inputdata;
/// Create the help keywords
  static void registerKeywords( Keywords& keys );
  explicit CLTool(const CLToolOptions& co );
/// Read the arguments from the command line
  bool readInput( int argc, char**argv, FILE* in, FILE*out );
/// virtual function mapping to the specific main for each tool
  virtual int main( FILE* in, FILE*out, Communicator&pc )=0;
/// virtual function returning a one-line descriptor for the tool
  virtual std::string description()const {return "(no description available)";}
/// virtual destructor to allow inheritance
  virtual ~CLTool() {}
};

template<class T>
bool CLTool::parse(const std::string&key,T&t) {
  plumed_massert(keywords.exists(key),"keyword " + key + " has not been registered");
  if(keywords.style(key,"compulsory") ) {
    if(inputData.count(key)==0) error("missing data for keyword " + key);
    bool check=Tools::convert(inputData[key],t);
    if(!check) error("data input for keyword " + key + " has wrong type");
    return true;
  }
  if( inputData.count(key)==0 ) return false;
  Tools::convert(inputData[key],t);
  return true;
}
// very limited support and check: take more from core/Action.h parseVector
template<class T>
bool CLTool::parseVector(const std::string&key,std::vector<T>&t) {

  // Check keyword has been registered
  plumed_massert(keywords.exists(key), "keyword " + key + " has not been registered");
  // initial size
  unsigned size=t.size();
  bool skipcheck=false;
  if(size==0) skipcheck=true; // if the vector in input has size zero, skip the check if size of input vector is the same of argument read

  // check if there is some value

  plumed_massert(inputData[key]!="false","compulsory keyword "+std::string(key)+"has no data");
  std::vector<std::string> words=Tools::getWords(inputData[key],"\t\n ,");
  t.resize(0);
  if(words.size()==0)return false;

  for(unsigned i=0; i<words.size(); ++i) {
    T v;
    Tools::convert(words[i],v);
    t.push_back(v);
  }
  // check the size
  if( !skipcheck && t.size()!=size ) {
    plumed_merror("vector read in for keyword  "+key+" has wrong size" );
  }
  std::string def;
  T val;
  if ( keywords.style(key,"compulsory") && t.size()==0 ) {
    if( keywords.getDefaultValue(key,def) ) {
      if( def.length()==0 || !Tools::convert(def,val) ) {
        plumed_merror("ERROR in keyword "+key+ " has weird default value" );
      } else {
        for(unsigned i=0; i<t.size(); ++i) t[i]=val;
      }
    } else {
      plumed_merror("keyword " + key + " is compulsory for this action");
    }
  }
  return true;
}



}

#endif
