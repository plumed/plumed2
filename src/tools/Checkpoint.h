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
#ifndef __PLUMED_tools_Checkpoint_h
#define __PLUMED_tools_Checkpoint_h

#include <string>
#include <vector>
#include <map>
#include <type_traits>

namespace PLMD {

/**
\ingroup TOOLBOX
Writer for the PLUMED checkpoint file.

The checkpoint is a sequence of named, typed entries (64-bit integers, doubles,
strings, and arrays of integers or doubles with a shape) grouped in nested
sections. Actions write their own state in a section named after their label
and never see the storage format, which is a minimal little-endian binary
layout described in Checkpoint.cpp. The content is accumulated in memory and
written to disk atomically by writeFile().
*/
class CheckpointWriter {
  std::string buffer;
  std::vector<std::string> sections;
  std::vector<unsigned> entriesInSection;
  std::string currentPath() const;
  void putName(const std::string& name);
  void putShape(const std::vector<unsigned long long>& dims);
  void writeIntArray(const std::string& name,const long long* data,const std::vector<unsigned long long>& dims);
  void writeDoubleArray(const std::string& name,const double* data,const std::vector<unsigned long long>& dims);
public:
/// plumedVersion is stored in the header, for inspection only.
  explicit CheckpointWriter(const std::string& plumedVersion="");
/// Open a section. Sections can be nested.
  void beginSection(const std::string& name);
/// Close the current section.
  void endSection();
/// Number of entries written so far in the current section (sub-sections count as one).
  unsigned getNumberOfEntriesInSection() const;
  void write(const std::string& name,double value);
  void write(const std::string& name,long long value);
  void write(const std::string& name,const std::string& value);
  void write(const std::string& name,const char* value) {
    write(name,std::string(value));
  }
  template<typename T>
  typename std::enable_if<std::is_integral<T>::value && !std::is_same<T,long long>::value>::type
  write(const std::string& name,T value) {
    write(name,static_cast<long long>(value));
  }
  template<typename T>
  typename std::enable_if<std::is_floating_point<T>::value && !std::is_same<T,double>::value>::type
  write(const std::string& name,T value) {
    write(name,static_cast<double>(value));
  }
/// Write an array with an explicit shape (dims empty means one-dimensional).
  template<typename T>
  void write(const std::string& name,const std::vector<T>& value,std::vector<unsigned long long> dims= {}) {
    if(dims.empty()) {
      dims.push_back(value.size());
    }
    if(std::is_integral<T>::value) {
      std::vector<long long> tmp(value.begin(),value.end());
      writeIntArray(name,tmp.data(),dims);
    } else {
      std::vector<double> tmp(value.begin(),value.end());
      writeDoubleArray(name,tmp.data(),dims);
    }
  }
/// Return the complete file content (header, entries, trailer with checksum).
/// All sections must be closed.
  std::string finalize() const;
/// Write the content to path atomically: the data go to a temporary file,
/// which is then renamed to path; the previous file, if any, is kept as prevPath
/// (no backup if prevPath is empty).
  void writeFile(const std::string& path,const std::string& prevPath) const;
};

/**
\ingroup TOOLBOX
Reader for the PLUMED checkpoint file written by CheckpointWriter.

The whole content is parsed on construction, which detects truncated or
corrupted files. Entries are then accessed by name, relative to the current
section; unknown entries and sections are ignored, so that newer files can be
read by older code as long as the entries they need are present.
*/
class CheckpointReader {
public:
  enum class Type { Int, Double, String };
  struct Entry {
    // parseEntry assigns the type for every successfully parsed entry.
    // cppcheck-suppress uninitMemberVarNoCtor
    Type type;
    std::vector<unsigned long long> dims;
    std::vector<long long> ints;
    std::vector<double> doubles;
    std::string str;
    bool scalar() const {
      return dims.empty();
    }
  };
private:
  std::map<std::string,Entry> entries;
  std::vector<std::string> order;
  std::vector<std::string> sectionList;
  std::vector<std::string> sections;
  unsigned formatVersion=0;
  std::string plumedVersion;
  std::string currentPath() const;
  const Entry& find(const std::string& name,Type type,bool scalar) const;
public:
/// Parse the content of a checkpoint file. Throws on any inconsistency.
  explicit CheckpointReader(const std::string& content);
/// Read a file and parse it.
  static std::string readFile(const std::string& path);
/// CRC32 of a string, the checksum used in the checkpoint file.
  static unsigned checksum(const std::string& data);
  unsigned getFormatVersion() const {
    return formatVersion;
  }
  const std::string& getPlumedVersion() const {
    return plumedVersion;
  }
  void beginSection(const std::string& name);
  void endSection();
/// True if a section with this name exists inside the current section.
  bool hasSection(const std::string& name) const;
/// True if an entry with this name exists in the current section.
  bool has(const std::string& name) const;
  void read(const std::string& name,double& value) const;
  void read(const std::string& name,long long& value) const;
  void read(const std::string& name,std::string& value) const;
  template<typename T>
  typename std::enable_if<std::is_integral<T>::value && !std::is_same<T,long long>::value>::type
  read(const std::string& name,T& value) const {
    long long tmp;
    read(name,tmp);
    value=static_cast<T>(tmp);
  }
  template<typename T>
  typename std::enable_if<std::is_floating_point<T>::value && !std::is_same<T,double>::value>::type
  read(const std::string& name,T& value) const {
    double tmp;
    read(name,tmp);
    value=static_cast<T>(tmp);
  }
/// Read an array; value is resized to the stored number of elements.
  template<typename T>
  void read(const std::string& name,std::vector<T>& value) const {
    if(std::is_integral<T>::value) {
      const Entry& e=find(name,Type::Int,false);
      value.assign(e.ints.begin(),e.ints.end());
    } else {
      const Entry& e=find(name,Type::Double,false);
      value.assign(e.doubles.begin(),e.doubles.end());
    }
  }
/// Shape of a stored array.
  std::vector<unsigned long long> getShape(const std::string& name) const;
  template<typename T>
  T get(const std::string& name) const {
    T value;
    read(name,value);
    return value;
  }
/// Full names of all sections, in file order (for inspection tools).
  const std::vector<std::string>& getSections() const {
    return sectionList;
  }
/// Full names of all entries, in file order (for inspection tools).
  const std::vector<std::string>& getEntryNames() const {
    return order;
  }
  const Entry& getEntry(const std::string& fullName) const;
};

}

#endif
