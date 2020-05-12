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
#ifndef __PLUMED_tools_OFile_h
#define __PLUMED_tools_OFile_h

#include "FileBase.h"
#include <vector>
#include <sstream>
#include <memory>

namespace PLMD {

class Value;

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
  pof.open("ciao");
  pof.printf("%s\n","test1");
  pof.setLinePrefix("plumed: ");
  pof.printf("%s\n","test2");
  pof.setLinePrefix("");
  pof.addConstantField("x2").printField("x2",67.0);
  pof.printField("x1",10.0).printField("x3",20.12345678901234567890).printField();
  pof.printField("x1",10.0).printField("x3",-1e70*20.12345678901234567890).printField();
  pof.printField("x3",10.0).printField("x2",777.0).printField("x1",-1e70*20.12345678901234567890).printField();
  pof.printField("x3",67.0).printField("x1",18.0).printField();
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

\section using-correctly-ofile Using correctly OFile in PLUMED

When a OFile object is used in PLUMED it can be convenient to link() it
to the Action object where it is defined, or to the PlumedMain object.
This will save in the OFile a pointer to the linked object and will
allow to have some extra information. E.g., if PLUMED is restarting,
files will be appended. Notice that one can enforce this behavior using
the enforceRestart() method before opening a file.

To have all files managed consistently, it is important to use OFile in the proper way.
This should allow multi-replica plumed, restart and backups to work in
the expected way. For this reason all the operations in OFile and IFile
are synchronizing all the processors of the group, so call to OFile functions
should always be performed by all processes; for this reason is also not usefull
to use Log for debugging because only master threads will actually write.
For debugging is better to use the standard stderr.

\verbatim
int main(){
// this is a growing file, containing a full history
// (frames are appended, as in traditional HILLS and COLVAR)
  OFile grw;
// this is a single-snapshopt file used e.g. for checkpointing
// (rewritten every time)
  OFile snp;

// open both files at the beginning
// (will go in \ref Action constructor)
  grw.open("growing");
  snp.open("snapshot");

// trajectory loop
  for(int i=0;i<nsteps;i++){

// files should be writen in the update() method of an \ref Action

// write on growing file
    grw<<"data at step "<<i<<\n";

// flushing
// it takes time, so do it only if data is critical
// better to leave this choice to the user with the FLUSH keyword
//    grw.flush();

// write on snapshot file
    snp.rewind();
    snp<<"snapshot at step "<<i<<"\n";
    snp.flush();
// the only difference is that snp is rewound
// notice that it should be rewound just before writing
// because rewind is going to move the file out of the way
// to have a safe copy of the file ("bck.last.filename")
// Also notice that snapshots should be flushed
// for this reason, it is better to write them only
// rarely to avoid excessive slow down

  }
}

\notice
Notice that it is not necessary to explicitely close files, since they are closed implicitely
when the object goes out of scope. In case you need to explicitly close the file before it is
destroyed, please check it the procedure is exception safe and, if necessary, add some `try/catch`
statement.

\endverbatim
*/

class OFile:
  public virtual FileBase {
/// Pointer to a linked OFile.
/// see link(OFile&)
  OFile* linked;
/// Internal buffer for printf
  std::unique_ptr<char[]> buffer_string;
/// Internal buffer (generic use)
  std::unique_ptr<char[]> buffer;
/// Internal buffer length
  int buflen;
/// This variables stores the actual buffer length
  int actual_buffer_length;
/// Class identifying a single field for fielded output
  class Field:
    public FieldBase {
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
/// The string used for backing up files
  std::string backstring;
/// Find field index given name
  unsigned findField(const std::string&name)const;
/// check if we are restarting
  bool checkRestart()const;
/// True if restart behavior should be forced
  bool enforceRestart_;
/// True if backup behavior (i.e. non restart) should be forced
  bool enforceBackup_;
public:
/// Constructor
  OFile();
/// Allows overloading of link
  using FileBase::link;
/// Allows overloading of open
  using FileBase::open;
/// Allows linking this OFile to another one.
/// In this way, everything written to this OFile will be immediately
/// written on the linked OFile. Notice that a OFile should
/// be either opened explicitly, linked to a FILE or linked to a OFile
  OFile& link(OFile&);
/// Set the string name to be used for automatic backup
  void setBackupString( const std::string& );
/// Backup a file by giving it a different name
  void backupFile( const std::string& bstring, const std::string& fname );
/// This backs up all the files that would have been created with the
/// name str.  It is used in analysis when you are not restarting.  Analysis
/// output files at different times, which are names analysis.0.<filename>,
/// analysis.1.<filename> and <filename>, are backed up to bck.0.analysis.0.<filename>,
/// bck.0.analysis.1.<filename> and bck.0.<filename>
  void backupAllFiles( const std::string& str );
/// Opens the file using automatic append/backup
  OFile& open(const std::string&name) override;
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
/// Rewind a file
  OFile&rewind();
/// Flush a file
  FileBase&flush() override;
/// Enforce restart, also if the attached plumed object is not restarting.
/// Useful for tests
  OFile&enforceRestart();
/// Enforce backup, even if the attached plumed object is restarting.
  OFile&enforceBackup();
};

/// Write using << syntax
template <class T>
OFile& operator<<(OFile&of,const T &t) {
  of.oss<<t;
  of.printf("%s",of.oss.str().c_str());
  of.oss.str("");
  return of;
}


}

#endif
