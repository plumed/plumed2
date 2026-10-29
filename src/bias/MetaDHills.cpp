/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2026 The plumed team
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
#include "MetaDHills.h"
#include "tools/Communicator.h"
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <system_error>

namespace PLMD {
namespace bias {

std::string exactString(double x) {
  char buffer[64];
  for(int precision=15; precision<17; precision++) {
    std::snprintf(buffer,sizeof(buffer),"%.*g",precision,x);
    if(std::strtod(buffer,nullptr)==x) {
      return buffer;
    }
  }
  std::snprintf(buffer,sizeof(buffer),"%.17g",x);
  return buffer;
}

void bcastString(Communicator& c, std::string& s) {
  unsigned n=s.size();
  c.Bcast(n,0);
  s.resize(n);
  if(n>0) {
    c.Bcast(s,0);
  }
}

std::string removeIncompleteLastLine(const std::string& path) {
  // compressed files cannot be fixed in place
  if(path.length()>3 && path.substr(path.length()-3)==".gz") {
    return "";
  }
  std::ifstream f(path,std::ios::binary);
  if(!f) {
    return "";
  }
  f.seekg(0,std::ios::end);
  const std::streamoff size=f.tellg();
  if(size<=0) {
    return "";
  }
  char last;
  f.seekg(size-1);
  f.get(last);
  if(last=='\n') {
    return "";
  }
  // look backward for the end of the last complete line
  std::streamoff pos=size-1;
  while(pos>0) {
    char c;
    f.seekg(pos-1);
    f.get(c);
    if(c=='\n') {
      break;
    }
    pos--;
  }
  std::string tail(size-pos,' ');
  f.seekg(pos);
  f.read(&tail[0],size-pos);
  f.close();
  std::string message="  WARNING: the last line of "+path+" is incomplete (the previous run was probably interrupted while writing it) and is removed:\n";
  message+="  "+tail+"\n";
  std::error_code ec;
  std::filesystem::resize_file(path,pos,ec);
  if(ec) {
    message+="  WARNING: the incomplete last line of "+path+" could not be removed ("+ec.message()+")\n";
  }
  return message;
}

}
}
