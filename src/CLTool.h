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
#include <cstdio>
#include "Tools.h"
#include "Keywords.h"

namespace PLMD{

class PlumedCommunicator;

class CLToolOptions{
  friend class CLTool;
  friend class CLToolRegister;
private:
  std::vector<std::string> line;
/// The documentation for this command line tool
  const Keywords& keys;
  static Keywords emptyKeys;
public:
  CLToolOptions(const std::string &name);
  CLToolOptions(const CLToolOptions& co, const Keywords& k);
};

/**
\ingroup INHERIT
This is the abstract base class to use for implementing new command line tool, within it there
is information as to how to go about implemneting a new tool.

To implement a command line tool you need to create a single cpp file call CLToolNAME.cpp. You can, in a command line
tool, use functionality from plumed to perform simple post-processing tasks.  For example, sum_hills uses the 
functionality inside from inside the biasing PLMD::Action, PLMD::BiasMetaD to calculate free energy surfaces.
Regardless, of what you are endeavoring to do your CLToolNAME.cpp file should be formatted in accordance with
the following template:

\verbatim
#include "CLTool.h"
#include "CLToolRegister.h"
#include "PlumedConfig.h"
#include "ActionRegister.h"

using namespace std;

namespace PLMD {
/**
//+PLUMEDOC TOOLS name
\endverbatim

Insert the documentation for your new tool here

\verbatim
\par Examples
\endverbatim

Insert some examples of how to use your tool here 

\verbatim
*/
//+ENDPLUMEDOC

/******* This is how you should declare the class for your command line tool.  main() does
         all the analsysis you require. The constructor and the registerKeywords routine 
         only do anything if you are using one of the template forms of input described below.
         However, reigsterKeywords() must always be declared.

class CLToolNAME:
public CLTool
{
public:
  static void registerKeywords( Keywords& keys );
  CLToolNAME(const CLToolOptions& co );
  int main(FILE* in, FILE*out,PlumedCommunicator& pc);
  string description()const{
    return "a description of the particular kind of task CLToolNAME performs";
  }
};

PLUMED_REGISTER_CLTOOL(CLToolNAME,"name")
\endverbatim

Insert the code for main, registerKeywords and the constructor here.

\verbatim
}  /--- don't forget to close the namespace PLMD { at the end of the file
\endverbatim

\section input Input

The input stream is passed to main so you can create tools that have an input in any form.
However, there are certain standardized forms of input that can be used for command line
tools. If your tool takes its input in one of these forms we strongly encourage you to
use the code that is already present.  If you do so a great deal of manual generation will 
be looked after automatically.

There are two command line input types that are implemented in the base class.  The first
is for tools such as driver where the input is specified using a series of command line flags
e.g.

\verbatim
plumed driver --plumed plumed.dat --ixyz trajectory.xyz --dump-forces
\endverbatim

The other are tools like simplmd that take an input file that contains one directive per line
and a corresponding value for that directive.  For both these forms of input it is possible to
read in everything that is required to run the calculation prior to the actual running of the 
calculation.  In other words these the user is not prompted to provide input data once the main 
calculation has started running (N.B. you can do tools with input of this sort though as the input stream
is passed to main).  

If you wish to use driver-like or simple-md like input then you have to specify this in the constructor.
For driver-like input you would write the following for the constructor:

\verbatim
CLToolNAME::CLToolNAME(const CLToolOptions& co ):
CLTool(co)
{
 inputdata=commandline;
}
\endverbatim

For simplemd-like input you write the following for the constructor:

\verbatim
CLToolNAME( const CLToolOptions& co ) :
CLTool(co)
{
  inputdata=ifile;
}
\endverbatim

If you are not using one of these input forms then you don't need to write a constructor although
you may choose to for reasons not connected to the input of data.

Once you have created the constructor the actual readin and manual creation is done in much the same
manner as it is done for Actions in the main plumed code (\ref usingDoxygen). You write a
registerKeywords routine as follows:

\verbatim
void CLToolNAME::registerKeywords( Keywords& keys ){
  CLTool::registerKeywords( keys );
}  
\endverbatim

Inside this routine you add descriptions of all your various command-line flags (driver-like) or
input directives (simple-md-like) and these descriptions are used to create the manual. The code 
will automatically check if a particular flag is present and read any input directive connected 
with a particular flag (i.e. the data after the space). Within main you can recover the read in 
data using CLTool::parse and CLTool:parseFlag.   

\section getplumed Re-using plumed

To re-use the functionality that is present in plumed you use the same tools that are used to 
patch the various MD codes (\ref HowToPlumedYourMD).  Alternatively, if you want to create
an instance of a particular Action you can do so by issuing the following commands:

\verbatim
PlumedMain* plumed=new PlumedMain(); std::vector<std::string> words;
Action* action=actionRegister().create(ActionOptions(plumed,words));
delete plumed; delete action;
\endverbatim

Please be aware that words should contain everything that would be required in an input
line to make a valid instance of the Action you require.
*/

class CLTool {
// This is a fried so we can create input data quickly
// when you do debug-float. 
template <typename real>
friend class CLToolDriver;
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
protected:
/// Get the value of one of the command line arguments 
  template<class T>
  bool parse(const std::string&key,T&t);
/// Find out whether one of the command line flags is present or not
  void parseFlag(const std::string&key,bool&t);  
public:
/// How is the input specified on the command line or in an input file
  enum {unset,commandline,ifile} inputdata;
/// Create the help keywords 
  static void registerKeywords( Keywords& keys );
  CLTool(const CLToolOptions& co ); 
/// Read the arguments from the command line
  bool readInput( int argc, char**argv, FILE* in, FILE*out );
/// virtual function mapping to the specific main for each tool
  virtual int main( FILE* in, FILE*out, PlumedCommunicator&pc )=0;
/// virtual function returning a one-line descriptor for the tool
  virtual std::string description()const{return "(no description available)";};
/// virtual destructor to allow inheritance
  virtual ~CLTool(){};
};

template<class T>
bool CLTool::parse(const std::string&key,T&t){
  plumed_massert(keywords.exists(key),"keyword " + key + " has not been registered");
  if(keywords.style(key,"compulsory") ){
     plumed_assert(inputData.count(key)>0);
     plumed_assert( Tools::convert(inputData[key],t) );
     return true;
  }
  if( inputData.count(key)==0 ) return false; 
  Tools::convert(inputData[key],t);
  return true;
}

}

#endif
