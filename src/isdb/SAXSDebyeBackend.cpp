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
#include "SAXSDebyeBackend.h"

#include <algorithm>

namespace PLMD {
namespace isdb {

SAXSDebyeBackendRegister& saxsDebyeBackendRegister() {
  static SAXSDebyeBackendRegister ans;
  return ans;
}

namespace {
// The register must exist before a plugin is loaded: a register created while
// a library is being registered (during LOAD) would deadlock on the registration lock.
[[maybe_unused]] const SAXSDebyeBackendRegister& registerAtLoad=saxsDebyeBackendRegister();

// the backends registered by a plugin are stored as "<image>:<name>"
std::string stripImage(const std::string& key) {
  const auto pos=key.rfind(':');
  return pos==std::string::npos ? key : key.substr(pos+1);
}
}

std::vector<std::string> SAXSDebyeBackendRegister::getBackendNames() const {
  std::vector<std::string> names;
  for(const auto & k : getKeys()) {
    const auto name=stripImage(k);
    if(std::find(names.begin(),names.end(),name)==names.end()) {
      names.push_back(name);
    }
  }
  return names;
}

std::unique_ptr<SAXSDebyeBackend> SAXSDebyeBackendRegister::create(const std::string& name, Action& action) const {
  const auto names=getBackendNames();
  std::string list;
  for(const auto & n : names) {
    list+=" "+n;
  }
  std::string selected=name;
  if(selected.empty()) {
    if(names.empty()) {
      throw ExceptionRegisterError()<<"no GPU backend for SAXS/SANS is available: it is provided by the plugin in plugins/openaccPTM of the PLUMED sources, please compile it and load it with LOAD before this action";
    }
    if(names.size()>1) {
      throw ExceptionRegisterError()<<"more than one GPU backend for SAXS/SANS is available, select one with GPU_BACKEND:"<<list;
    }
    selected=names[0];
  }
  // the last registered key with this name (a plugin loaded twice registers it twice)
  std::string key;
  for(const auto & k : getKeys()) {
    if(stripImage(k)==selected) {
      key=k;
    }
  }
  if(key.empty()) {
    throw ExceptionRegisterError()<<"GPU backend for SAXS/SANS \""<<selected<<"\" is not available, the available ones are:"<<(list.empty()?std::string(" none"):list);
  }
  return get(key).create(action);
}

}
}
