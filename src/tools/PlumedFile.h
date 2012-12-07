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
#ifndef __PLUMED_tools_PlumedFile_h
#define __PLUMED_tools_PlumedFile_h

#include <cstdio>
#include <vector>
#include <string>
#include <sstream>

namespace PLMD{

class PlumedCommunicator;
class PlumedMain;
class Action;
class Value;

/**
Base class for dealing with files.

This class just provides things which are common among PlumedOFile and PlumedIFile
*/

class PlumedFileBase{
/// Copy constructor is disabled (private and unimplemented)
  PlumedFileBase(const PlumedFileBase&);
/// Assignment operator is disabled (private and unimplemented)
  PlumedFileBase& operator=(const PlumedFileBase&);
protected:
/// Internal tool.
/// Base for PlumedIFile::Field and PlumedOFile::Field
  class FieldBase{
// everything is public to simplify usage
  public:
    std::string name;
    std::string value;
    bool constant;
    FieldBase(): constant(false){}
  };

/// file pointer
  FILE* fp;
/// communicator. NULL if not set
  PlumedCommunicator* comm;
/// pointer to main plumed object. NULL if not linked
  PlumedMain* plumed;
/// pointer to corresponding action. NULL if not linked
  Action* action;
/// Control closing on destructor.
/// If true, file will not be closed in destructor
  bool cloned;
/// Private constructor.
/// In this manner one cannot instantiate a PlumedFileBase object
  PlumedFileBase();
/// Set to true when end of file is encountered
  bool eof;
/// Set to true when error is encountered
  bool err;
/// path of the opened file
  std::string path;
/// Set to true if you want flush to be heavy (close/reopen)
  bool heavyFlush;
public:
/// Link to an already open filed
  PlumedFileBase& link(FILE*);
/// Link to a PlumedMain object
/// Automatically links also the corresponding PlumedCommunicator.
  PlumedFileBase& link(PlumedMain&);
/// Link to a PlumedCommunicator object
  PlumedFileBase& link(PlumedCommunicator&);
/// Link to an Action object.
/// Automatically links also the corresponding PlumedMain and PlumedCommunicator.
  PlumedFileBase& link(Action&);
/// Flushes the file to disk
  PlumedFileBase& flush();
/// Closes the file
/// Should be used only for explicitely opened files.
  void        close();
/// Virtual destructor (allows inheritance)
  virtual ~PlumedFileBase();
/// Runs a small testcase
  static void test();
/// Check for error/eof.
  operator bool () const;
/// Set heavyFlush flag
  void setHeavyFlush(){ heavyFlush=true;};
/// Opens the file (without auto-backup)
  PlumedFileBase& open(const std::string&name,const std::string& mode);
/// Check if the file exists
  bool FileExist(const std::string& path);
/// Check if a file is open
  bool isOpen();
};

/**
\ingroup TOOLBOX
Class for output files

This class provides features similar to those in the standard C "FILE*" type,
but only for sequential output. See PlumedIFile for sequential input.

See the example here for a possible use:
\verbatim
#include "PlumedFile.h"

int main(){
  PLMD::PlumedOFile pof;
  pof.open("ciao","w");
  pof.printf("%s\n","test1");
  pof.setLinePrefix("plumed: ");
  pof.printf("%s\n","test2");
  pof.setLinePrefix("");
  pof.addConstantField("x2").printField("x2",67.0);
  pof.printField("x1",10.0).printField("x3",20.12345678901234567890).printField();
  pof.printField("x1",10.0).printField("x3",-1e70*20.12345678901234567890).printField();
  pof.printField("x3",10.0).printField("x2",777.0).printField("x1",-1e70*20.12345678901234567890).printField();
  pof.printField("x3",67.0).printField("x1",18.0).printField();
  pof.close();
  return 0;
}
\endverbatim

This program is expected to produce a file "ciao" which reads
\verbatim
test1
plumed: test2
#! FIELDS x1 x3
#! SET x2                      67
                     10      20.12345678901234
                     10 -2.012345678901235e+71
#! FIELDS x1 x3
#! SET x2                     777
 -2.012345678901235e+71                     10
                     18                     67
\endverbatim

Notes
- "x2" is declared as "constant", which means that it is written using the "SET"
keyword. Thus, everytime it is modified, all the headers are repeated in the output file.
- printField() without arguments is used as a "newline".
- most methods return a reference to the PlumedOFile itself, to allow chaining many calls on the same line
(this is similar to << operator in std::ostream)

*/

class PlumedOFile:
public virtual PlumedFileBase{
/// Pointer to a linked OFile.
/// see link(PlumedOFile&)
  PlumedOFile* linked;
/// Internal buffer for printf
  char* buffer_string;
/// Internal buffer (generic use)
  char* buffer;
/// Internal buffer length
  int buflen;
/// Class identifying a single field for fielded output
  class Field:
  public FieldBase{
  };
/// Low-level write
  size_t llwrite(const char*,size_t);
/// True if fields has changed.
/// This could be due to a change in the list of fields or a reset
/// of a nominally constant field
  bool fieldChanged;
/// Format for fields writing
  std::string fieldFmt;
/// All the previously defined variable fields
  std::vector<Field> previous_fields;
/// All the defined variable fields
  std::vector<Field> fields;
/// All the defined constant fields
  std::vector<Field> const_fields;
/// Prefix for line (e.g. "PLUMED: ")
  std::string linePrefix;
/// Temporary ostringstream for << output
  std::ostringstream oss;
/// Find field index given name
  unsigned findField(const std::string&name)const;
public:
/// Constructor
  PlumedOFile();
/// Destructor
  ~PlumedOFile();
/// Allows overloading of link
  using PlumedFileBase::link;
/// Allows overloading of open
  using PlumedFileBase::open;
/// Allows linking this PlumedOFile to another one.
/// In this way, everything written to this PlumedOFile will be immediately
/// written on the linked PlumedOFile. Notice that a PlumedOFile should
/// be either opened explicitly, linked to a FILE or linked to a PlumedOFile
  PlumedOFile& link(PlumedOFile&);
/// Opens the file using automatic append/backup
  PlumedOFile& open(const std::string&name);
/// Set the prefix for output.
/// Typically "PLUMED: ". Notice that lines with a prefix cannot
/// be parsed using fields in a PlumedIFile.
  PlumedOFile& setLinePrefix(const std::string&);
/// Set the format for writing double precision fields
  PlumedOFile& fmtField(const std::string&);
/// Reset the format for writing double precision fields to its default
  PlumedOFile& fmtField();
/// Set the value of a double precision field
  PlumedOFile& printField(const std::string&,double);
/// Set the value of a int field
  PlumedOFile& printField(const std::string&,int);
/// Set the value of a string field
  PlumedOFile& printField(const std::string&,const std::string&);
///
  PlumedOFile& addConstantField(const std::string&);
/// Used to setup printing of values
  PlumedOFile& setupPrintValue( Value *val );
/// Print a value
  PlumedOFile& printField( Value* val, const double& v );
/** Close a line.
Typically used as
\verbatim
  of.printField("a",a).printField("b",b).printField();
\endverbatim
*/
  PlumedOFile& printField();
/**
Resets the list of fields.
As it is only possible to add new constant fields (addConstantField()),
this method can be used to clean the field list.
*/
  PlumedOFile& clearFields();
/// Formatted output with explicit format - a la printf
  int printf(const char*fmt,...);
/// Formatted output with << operator
  template <class T>
  friend PlumedOFile& operator<<(PlumedOFile&,const T &);
};


/**
\ingroup TOOLBOX
Class for input files

This class provides features similar to those in the standard C "FILE*" type,
but only for sequential input. See PlumedOFile for sequential output.

*/
class PlumedIFile:
/// Class identifying a single field for fielded output
public virtual PlumedFileBase{
  class Field:
  public FieldBase{
  public:
    bool read;
    Field(): read(false) {}
  };
/// Low-level read.
/// Note: in parallel, all processes read
  size_t llread(char*,size_t);
/// All the defined fields
  std::vector<Field> fields;
/// Flag set in the middle of a field reading
  bool inMiddleOfField;
/// Set to true if you want to allow fields to be ignored in the read in file
  bool ignoreFields;
/// Advance to next field (= read one line)
  PlumedIFile& advanceField();
/// Find field index by name
  unsigned findField(const std::string&name)const;
public:
/// Constructor
  PlumedIFile();
/// Destructor
  ~PlumedIFile();
/// Opens the file 
  PlumedIFile& open(const std::string&name);
/// Gets the list of all fields
  PlumedIFile& scanFieldList(std::vector<std::string>&);
/// Read a double field
  PlumedIFile& scanField(const std::string&,double&);
/// Read a int field
  PlumedIFile& scanField(const std::string&,int&);
/// Read a string field
  PlumedIFile& scanField(const std::string&,std::string&);
/**
 Ends a field-formatted line.

Typically used as
\verbatim
  if.scanField("a",a).scanField("b",b).scanField();
\endverbatim
*/
  PlumedIFile& scanField();
/// Get a full line as a string
  PlumedIFile& getline(std::string&);
/// Reset end of file                                                              
  void reset(bool);
/// Check if a field exist                                                       
  bool FieldExist(const std::string& s);
/// Read in a value
  PlumedIFile& scanField(Value* val);
/// Allow some of the fields in the input to be ignored
  void allowIgnoredFields();
};

/// Write using << syntax
template <class T>
PlumedOFile& operator<<(PlumedOFile&of,const T &t){
  of.oss<<t;
  of.printf("%s",of.oss.str().c_str());
  of.oss.str("");
  return of;
}


}

#endif
