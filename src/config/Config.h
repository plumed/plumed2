/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2015 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

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

namespace PLMD{
namespace config{

std::string getSoExt();

std::string getPlumedRoot();

std::string getPlumedHtmldir();

std::string getPlumedIncludedir();

std::string getPlumedProgramName();

std::string getEnvCommand();

std::string getMakefile();

std::string getVersion();

std::string getVersionLong();

std::string getVersionGit();

bool hasMatheval();

bool hasDlopen();

bool hasAlmost();

bool isInstalled();

bool hasCregex();

bool hasMolfile();

bool hasExternalMolfile();

bool hasZlib();

bool hasXdrfile();
}
}

#endif
