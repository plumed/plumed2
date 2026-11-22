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
#ifndef __PLUMED_tools_H5Archive_h
#define __PLUMED_tools_H5Archive_h
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace PLMD {
namespace h5table {
/// Restricted, uncompressed USTAR container of immutable PLUMED HDF5 tables.
class Archive {
  struct Impl;
  std::unique_ptr<Impl> impl;
public:
  struct Member {
    std::uint64_t offset=0,size=0;
    bool recovery=false;
    std::uint64_t block=0,rows=0;
  };
  explicit Archive(const std::string& path,bool writing=false,bool append=true);
  ~Archive();
  Archive(const Archive&)=delete;
  Archive& operator=(const Archive&)=delete;
  const std::vector<Member>& members() const;
  std::vector<char> image(std::size_t index) const;
  void add(const std::vector<char>& image,std::uint64_t block,std::uint64_t rows);
  void finish();
  /// Replace only the tail, using a durable, replayable recovery record.
  void replaceTail(std::size_t first,const std::vector<char>& image,std::uint64_t block,std::uint64_t rows);
  static bool recognizes(const std::string& path);
};

/// Private staging directory; only one batch, never the accumulated archive.
class BatchFile {
  std::string directory;
public:
  std::string path;
  explicit BatchFile(const std::string& parent="");
  ~BatchFile();
  BatchFile(const BatchFile&)=delete;
  BatchFile& operator=(const BatchFile&)=delete;
  std::vector<char> image() const;
  /// Synchronize and atomically replace a destination, preserving its permissions.
  void publish(const std::string& destination);
};
}
}
#endif
