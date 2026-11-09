/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2026 The plumed team
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
#ifndef __PLUMED_tools_H5Table_h
#define __PLUMED_tools_H5Table_h

#include <memory>
#include <string>
#include <vector>

namespace PLMD {
/// (not "h5": libraries such as ITensor use a global namespace h5)
namespace h5table {

/**
\ingroup TOOLBOX
Tables of fields (the files written with OFile::printField) stored in HDF5 files.

A file whose name ends with `.h5` is written in HDF5 format (with libhdf5, see
`./configure --enable-libhdf5`). Layout, version 1:
- the root group has the attribute `plumed_format` = "plumed-table" and `plumed_format_version` = 1;
- the table is made of blocks, `block_0`, `block_1`, ...: a new block starts when the fields
  change (where a text file would print a new `#! FIELDS` header), i.e. when the names or the
  types of the columns or the values of the constant fields change;
- every column of a block is an extendable 1-D dataset with the name of the field (doubles,
  64-bit integers, or strings), in the order given by the attribute `plumed_columns` of the block
  (names separated by spaces);
- the constant fields (`#! SET` in text files) are attributes of the block, listed in order in
  the attribute `plumed_constants`.
Values are stored with full precision, whatever FMT says.

The reader (Reader) presents the table as the lines of the equivalent text file, with the numbers
written with all their digits, so that everything that reads text files (IFile) reads HDF5 files too.
*/

/// True if PLUMED was compiled with libhdf5
bool available();

/// True if the file is an HDF5 table (by its extension, .h5)
bool isH5(const std::string& path);

/// A value of a field
struct Value {
  enum Type {Double,Int,String};
  Type type=String;
  double d=0.0;
  long long i=0;
  std::string s;
  /// the value as written in a text file, with all the digits for the numbers
  std::string text() const;
  bool operator==(const Value& v) const;
};

/// A field: name and value
struct Field {
  std::string name;
  Value value;
};

/// Writes a table: rows are buffered and written in chunks (flush() writes them all)
class Writer {
  struct Impl;
  std::unique_ptr<Impl> impl;
public:
/// Create the file, or open it to append to its last block (append)
  Writer(const std::string& path,bool append);
  ~Writer();
  Writer(const Writer&)=delete;
  Writer& operator=(const Writer&)=delete;
/// Add a row with these fields and constants (a new block starts if they changed)
  void addRow(const std::vector<Field>& fields,const std::vector<Field>& constants);
/// Write the buffered rows and flush the file
  void flush();
/// Number of rows of every block in the file (after writing the buffered rows)
  std::vector<unsigned long long> rows();
};

/// Reads a table as the lines of the equivalent text file
class Reader {
  struct Impl;
  std::unique_ptr<Impl> impl;
public:
  explicit Reader(const std::string& path);
  ~Reader();
  Reader(const Reader&)=delete;
  Reader& operator=(const Reader&)=delete;
/// The next line; false at the end of the table
  bool getline(std::string& line);
};

/// Number of rows of every block of a table
std::vector<unsigned long long> rows(const std::string& path);

/// Keep the first keep[b] rows of block b (the blocks beyond keep.size() are removed).
/// The rows removed are written, with the same structure, in the table backup (if not empty).
void truncate(const std::string& path,const std::vector<unsigned long long>& keep,const std::string& backup);

}
}

#endif
