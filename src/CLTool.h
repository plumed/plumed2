/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

   This file is part of plumed, version 2.0.

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

#ifndef __PLUMED_CLTool_h
#define __PLUMED_CLTool_h
#include <cstdio>
#include <vector>
#include <string>

namespace PLMD{

class PlumedCommunicator;

class CLToolOptions{
  friend class CLToolRegister;
  std::vector<std::string> line;
public:
  CLToolOptions(const std::string &name):
    line(1,name) { }
};

/**
Interface to all the command-line tools.

This class just define an interface, and does not implement anything.
Inherits from this class to create a new command-line tool.
*/
class CLTool
{
public:
/// virtual function mapping to the specific main for each tool
  virtual int main(int argc, char **argv,FILE*in,FILE*out,PlumedCommunicator&pc)=0;
/// virtual function returning a one-line descriptor for the tool
  virtual std::string description()const{return "(no description available)";};
/// virtual destructor to allow inheritance
  virtual ~CLTool(){};
};

}


#endif
