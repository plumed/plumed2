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
#ifndef __PLUMED_tools_File_h
#define __PLUMED_tools_File_h

#include <cstdio>
#include <vector>
#include <string>
#include <sstream>

namespace PLMD{

class Communicator;
class PlumedMain;
class Action;
class Value;

/**
Base class for dealing with files.

This class just provides things which are common among OFile and IFile
*/

class FileBase{
/// Copy constructor is disabled (private and unimplemented)
  FileBase(const FileBase&);
/// Assignment operator is disabled (private and unimplemented)
  FileBase& operator=(const FileBase&);
protected:
/// Internal tool.
/// Base for IFile::Field and OFile::Field
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
  Communicator* comm;
/// pointer to main plumed object. NULL if not linked
  PlumedMain* plumed;
/// pointer to corresponding action. NULL if not linked
  Action* action;
/// Control closing on destructor.
/// If true, file will not be closed in destructor
  bool cloned;
/// Private constructor.
/// In this manner one cannot instantiate a FileBase object
  FileBase();
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
  FileBase& link(FILE*);
/// Link to a PlumedMain object
/// Automatically links also the corresponding Communicator.
  FileBase& link(PlumedMain&);
/// Link to a Communicator object
  FileBase& link(Communicator&);
/// Link to an Action object.
/// Automatically links also the corresponding PlumedMain and Communicator.
  FileBase& link(Action&);
/// Flushes the file to disk
  FileBase& flush();
/// Closes the file
/// Should be used only for explicitely opened files.
  void        close();
/// Virtual destructor (allows inheritance)
  virtual ~FileBase();
/// Runs a small testcase
  static void test();
/// Check for error/eof.
  operator bool () const;
/// Set heavyFlush flag
  void setHeavyFlush(){ heavyFlush=true;};
/// Opens the file (without auto-backup)
  FileBase& open(const std::string&name,const std::string& mode);
/// Check if the file exists
  bool FileExist(const std::string& path);
/// Check if a file is open
  bool isOpen();
};

/**
\ingroup TOOLBOX
Class for output files

This class provides features similar to those in the standard C "FILE*" type,
but only for sequential output. See IFile for sequential input.

See the example here for a possible use:
\verbatim
#include "File.h"

int main(){
  PLMD::OFile pof;
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
- most methods return a reference to the OFile itself, to allow chaining many calls on the same line
(this is similar to << operator in std::ostream)

*/

class OFile:
public virtual FileBase{
/// Pointer to a linked OFile.
/// see link(OFile&)
  OFile* linked;
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
  OFile();
/// Destructor
  ~OFile();
/// Allows overloading of link
  using FileBase::link;
/// Allows overloading of open
  using FileBase::open;
/// Allows linking this OFile to another one.
/// In this way, everything written to this OFile will be immediately
/// written on the linked OFile. Notice that a OFile should
/// be either opened explicitly, linked to a FILE or linked to a OFile
  OFile& link(OFile&);
/// Opens the file using automatic append/backup
  OFile& open(const std::string&name);
/// Set the prefix for output.
/// Typically "PLUMED: ". Notice that lines with a prefix cannot
/// be parsed using fields in a IFile.
  OFile& setLinePrefix(const std::string&);
/// Set the format for writing double precision fields
  OFile& fmtField(const std::string&);
/// Reset the format for writing double precision fields to its default
  OFile& fmtField();
/// Set the value of a double precision field
  OFile& printField(const std::string&,double);
/// Set the value of a int field
  OFile& printField(const std::string&,int);
/// Set the value of a string field
  OFile& printField(const std::string&,const std::string&);
///
  OFile& addConstantField(const std::string&);
/// Used to setup printing of values
  OFile& setupPrintValue( Value *val );
/// Print a value
  OFile& printField( Value* val, const double& v );
/** Close a line.
Typically used as
\verbatim
  of.printField("a",a).printField("b",b).printField();
\endverbatim
*/
  OFile& printField();
/**
Resets the list of fields.
As it is only possible to add new constant fields (addConstantField()),
this method can be used to clean the field list.
*/
  OFile& clearFields();
/// Formatted output with explicit format - a la printf
  int printf(const char*fmt,...);
/// Formatted output with << operator
  template <class T>
  friend OFile& operator<<(OFile&,const T &);
};


/**
\ingroup TOOLBOX
Class for input files

This class provides features similar to those in the standard C "FILE*" type,
but only for sequential input. See OFile for sequential output.

*/
class IFile:
/// Class identifying a single field for fielded output
public virtual FileBase{
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
  IFile& advanceField();
/// Find field index by name
  unsigned findField(const std::string&name)const;
public:
/// Constructor
  IFile();
/// Destructor
  ~IFile();
/// Opens the file 
  IFile& open(const std::string&name);
/// Gets the list of all fields
  IFile& scanFieldList(std::vector<std::string>&);
/// Read a double field
  IFile& scanField(const std::string&,double&);
/// Read a int field
  IFile& scanField(const std::string&,int&);
/// Read a string field
  IFile& scanField(const std::string&,std::string&);
/**
 Ends a field-formatted line.

Typically used as
\verbatim
  if.scanField("a",a).scanField("b",b).scanField();
\endverbatim
*/
  IFile& scanField();
/// Get a full line as a string
  IFile& getline(std::string&);
/// Reset end of file                                                              
  void reset(bool);
/// Check if a field exist                                                       
  bool FieldExist(const std::string& s);
/// Read in a value
  IFile& scanField(Value* val);
/// Allow some of the fields in the input to be ignored
  void allowIgnoredFields();
};

/// Write using << syntax
template <class T>
OFile& operator<<(OFile&of,const T &t){
  of.oss<<t;
  of.printf("%s",of.oss.str().c_str());
  of.oss.str("");
  return of;
}


}

#endif
