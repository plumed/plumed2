/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2020 The plumed team
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
#ifndef __PLUMED_config_Config_h
#define __PLUMED_config_Config_h

#include <string>

namespace PLMD {
namespace config {

/// Return the extension of shared libraries on this system.
/// It is either "so" or "dylib". In case shared libraries are disabled, it returns an empty string.
std::string getSoExt();

/// Return path to the root of the PLUMED package.
std::string getPlumedRoot();

/// Return path to documentation.
/// User documentation is located in `getPlumedHtmldir()+"/user-doc/html/index.html"`
/// Developer documentation is located in `getPlumedHtmldir()+"/developer-doc/html/index.html"`
std::string getPlumedHtmldir();

/// Return path to the include directory.
/// The header file for PLUMED wrapper is in `getPlumedIncludedir()+getPlumedProgramName()+"/wrapper/Plumed.h"`
std::string getPlumedIncludedir();

/// Return the name used for installing PLUMED.
/// E.g. if PLUMED has been compiled with `./configure --program-suffix _mpi`
/// it returns "plumed_mpi"
std::string getPlumedProgramName();

/**
Return a string containing a sequence of environment variables.

The returned string has the form:
\verbatim
env PLUMED_ROOT=/path env PLUMED_HTMLDIR=/path ... etc
\endverbatim
This string is used internally in PLUMED to run scripts located in plumedroot/script.
For instance, the `patch` script can be run executing the following command:
\verbatim
config::getEnvCommand()+" \""+getPlumedRoot()+"\"/scripts/patch.sh";
\endverbatim
Notice that the getPlumedRoot() output is enclosed in escaped quotes. This
allows this directory to have spaces in its name.
*/
std::string getEnvCommand();

/// Return the content of `Makefile.conf` in a single string.
/// Can be used to inspect the variables used to customize PLUMED.
/// Notice that this reflects the content of the `Makefile.conf` file.
/// Since some of the there defined variables can be overwritten at install
/// (e.g., one can use `make install prefix=/new/path` to change the installation path)
/// their values could be not faithful
std::string getMakefile();

/// Return the short PLUMED version.
/// E.g. "2.2"
std::string getVersion();

/// Return the long PLUMED version.
/// E.g. "2.2.3"
std::string getVersionLong();

/// Return the git PLUMED version
/// E.g. "c5badb091cd30"
std::string getVersionGit();

/// Return the day PLUMED was compiled.
/// E.g. "Apr 16 2018"
std::string getCompilationDate();

/// Return the time at which PLUMED was compiled.
/// E.g. "13:27:58"
std::string getCompilationTime();

bool hasMatheval();

bool hasDlopen();

bool isInstalled();

bool hasCregex();

bool hasMolfile();

bool hasExternalMolfile();

bool hasZlib();

bool hasXdrfile();
}
}

#endif
