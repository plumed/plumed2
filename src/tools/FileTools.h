/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2025 The plumed team
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
#ifndef __PLUMED_tools_FileTools_h
#define __PLUMED_tools_FileTools_h
#include <memory>
#include <cstdio>
#include <filesystem>
#include <string>
namespace PLMD {

struct fileDeleter {
  void operator()(FILE* fp) {
    std::fclose(fp);
  }
};
// call fclose when this goes out of scope
using unique_FILE = std::unique_ptr<FILE,fileDeleter>;

/// \name File names in UTF-8
/// PLUMED stores file names as UTF-8 strings (as the python interface passes
/// them). On Windows the narrow functions of the C and C++ libraries use the
/// ANSI code page instead: these functions convert the names to UTF-16 there,
/// and are the standard functions on the other systems.
/// A name that is not valid UTF-8 (e.g. from a program using the ANSI code
/// page) is read in the ANSI code page, as the narrow functions would do.
///@{
/// std::fopen
std::FILE* fopenUTF8(const std::string& path,const char* mode);
/// std::rename (0 on success, errno set on failure)
int renameUTF8(const std::string& from,const std::string& to);
/// a std::filesystem::path
std::filesystem::path pathFromUTF8(const std::string& path);
/// the name of a std::filesystem::path, in UTF-8
std::string pathToUTF8(const std::filesystem::path& path);
/// std::getenv, false if the variable is not set
bool getenvUTF8(const char* name,std::string& value);
/// gzopen (returns a gzFile, nullptr if PLUMED is compiled without zlib)
void* gzopenUTF8(const std::string& path,const char* mode);
///@}
}
#endif //

