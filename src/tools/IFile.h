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
#ifndef __PLUMED_tools_IFile_h
#define __PLUMED_tools_IFile_h

#include "FileBase.h"
#include <vector>
#include <cstddef>

namespace PLMD {

class Value;

/**
\ingroup TOOLBOX
Class for input files

This class provides features similar to those in the standard C "FILE*" type,
but only for sequential input. See OFile for sequential output.

*/
class IFile:
/// Class identifying a single field for fielded output
  public virtual FileBase {
  class Field:
    public FieldBase {
  public:
    bool read;
    Field(): read(false) {}
  };
/// Low-level read.
/// Note: in parallel, all processes read
  std::size_t llread(char*,std::size_t);
/// All the defined fields
  std::vector<Field> fields;
/// Flag set in the middle of a field reading
  bool inMiddleOfField;
/// Set to true if you want to allow fields to be ignored in the read in file
  bool ignoreFields;
/// Set to true to allow files without end-of-line at the end
  bool noEOL;
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
  IFile& open(const std::string&name) override;
/// Gets the list of all fields
  IFile& scanFieldList(std::vector<std::string>&);
/// Read a double field
  IFile& scanField(const std::string&,double&);
/// Read a int field
  IFile& scanField(const std::string&,int&);
/// Read a long int field
  IFile& scanField(const std::string&,long int&);
/// Read a unsigned field
  IFile& scanField(const std::string&,unsigned&);
/// Read a long unsigned field
  IFile& scanField(const std::string&,long unsigned&);
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
/// Allow files without EOL at the end.
/// This in practice should be only used when opening
/// plumed input files
  void allowNoEOL();
};

}

#endif
