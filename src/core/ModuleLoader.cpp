/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2023 The plumed team
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
#include "ModuleLoader.h"
#include "ActionRegister.h"
#include "CLToolRegister.h"
#include "config/Config.h"
#include "tools/DLLoader.h"
#include "tools/Exception.h"
#include "tools/FileTools.h"
#include <algorithm>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <map>
#include <mutex>
#include <set>
#include <sstream>
#ifdef __PLUMED_HAS_DLOPEN
#include <dlfcn.h>
#endif

namespace PLMD {
struct ModuleLoaderContextData {
  std::vector<void*> images;
  std::vector<std::string> modules;
  explicit ModuleLoaderContextData(const DLLoader& loader):
    images(loader.getHandles()), modules(loader.getModuleNames()) {}
};

namespace {
namespace fs=std::filesystem;
// Context owns this heap snapshot and restores its parent before freeing it.
// Keep TLS trivial: a TLS destructor would pin the kernel on ELF.
thread_local const ModuleLoaderContextData* contextData=nullptr;
struct Module {
  std::string name, build, variant;
  fs::path library;
  std::vector<std::string> dependencies;
};

bool validName(const std::string& name) {
  return !name.empty() && name.find_first_not_of("abcdefghijklmnopqrstuvwxyz0123456789_-")==std::string::npos;
}

#ifdef __PLUMED_HAS_DLOPEN
std::string moduleSymbol(const std::string& name) {
  if(name.find('-')==std::string::npos) {
    return "plumed_module_"+name;
  }
  // A disjoint prefix preserves existing symbols and avoids aliasing a-b/a_b.
  std::string result="plumed_encoded_module_";
  constexpr char digits[]="0123456789abcdef";
  for(unsigned char c:name) {
    result+=digits[c>>4];
    result+=digits[c&15];
  }
  return result;
}
#endif

struct State {
  std::recursive_mutex mutex;
  std::once_flag discovery, automatic;
  std::map<std::string,Module> catalog;
  DLLoader retained;
  std::string activating;

  struct Activation {
    std::string& busy;
    explicit Activation(std::string& flag,const std::string& name):busy(flag) {
      plumed_assert(busy.empty())<<"Recursive module activation of "<<name<<" while loading "<<busy;
      busy=name;
    }
    ~Activation() {
      busy.clear();
    }
  };

  void discover() {
    std::call_once(discovery,[this] {
      std::vector<fs::path> directories;
      const auto kernel=config::getLibraryPath();
      if(!kernel.empty()) {
        directories.push_back(fs::canonical(pathFromUTF8(kernel)).parent_path()/"plumed"/"plugins");
      }
      if(const char* paths=std::getenv("PLUMED_PLUGIN_PATH")) {
        std::istringstream list(paths);
        std::string entry;
#ifdef _WIN32
        const char separator=';';
#else
        const char separator=':';
#endif
        while(std::getline(list,entry,separator)) {
          if(!entry.empty()) {
            directories.emplace_back(pathFromUTF8(entry));
          }
        }
      }
      // A troubleshooting override selects an exact library for an ambiguous
      // module. It never bypasses build compatibility or dependency checks.
      std::map<std::string,fs::path> overrides;
      if(const char* selections=std::getenv("PLUMED_MODULE_OVERRIDE")) {
        std::istringstream list(selections);
        std::string selection;
        while(std::getline(list,selection,';')) {
          auto equals=selection.find('=');
          plumed_assert(equals!=std::string::npos)<<"Expected module=/absolute/library in PLUMED_MODULE_OVERRIDE";
          overrides[selection.substr(0,equals)]=fs::canonical(pathFromUTF8(selection.substr(equals+1)));
        }
      }
      std::map<std::string,Module> found;
      std::set<fs::path> scanned;
      for(const auto& directory:directories) {
        if(!fs::exists(directory)) {
          continue;
        }
        if(!scanned.insert(fs::canonical(directory)).second) {
          continue;
        }
        std::vector<fs::path> files;
        for(const auto& entry:fs::directory_iterator(directory)) {
          if(entry.path().extension()==".plumed") {
            files.push_back(entry.path());
          }
        }
        std::sort(files.begin(),files.end());
        for(const auto& file:files) {
          Module module;
          std::ifstream input(file);
          plumed_assert(input)<<"Cannot read module registration "<<file.string();
          std::string line;
          while(std::getline(input,line)) {
            if(line.empty() || line[0]=='#') {
              continue;
            }
            const auto space=line.find(' ');
            plumed_assert(space!=std::string::npos)<<"Malformed module registration "<<file.string();
            auto key=line.substr(0,space), value=line.substr(space+1);
            if(key=="module") {
              module.name=value;
            } else if(key=="build") {
              module.build=value;
            } else if(key=="variant") {
              module.variant=value;
            } else if(key=="library") {
              module.library=file.parent_path()/pathFromUTF8(value);
            } else if(key=="requires") {
              std::istringstream deps(value);
              std::string dep;
              while(deps>>dep) {
                plumed_assert(validName(dep))<<"Invalid dependency in "<<file.string();
                module.dependencies.push_back(dep);
              }
            } else {
              plumed_error()<<"Unknown registration field "<<key<<" in "<<file.string();
            }
          }
          plumed_assert(validName(module.name) && !module.build.empty() && !module.library.empty())
              <<"Incomplete module registration "<<file.string();
          // Validate compatibility on selection, allowing explicit loading with
          // autoload disabled to avoid unrelated broken installations.
          auto override=overrides.find(module.name);
          if(override!=overrides.end() && fs::weakly_canonical(module.library)!=override->second) {
            continue;
          }
          auto previous=found.find(module.name);
          if(previous!=found.end()) {
            plumed_assert(fs::equivalent(previous->second.library,module.library)
                          && previous->second.build==module.build
                          && previous->second.variant==module.variant
                          && previous->second.dependencies==module.dependencies)
                <<"Ambiguous module "<<module.name<<"; select one installation or use PLUMED_MODULE_OVERRIDE";
          } else {
            found.emplace(module.name,std::move(module));
          }
        }
      }
      for(const auto& override:overrides) {
        plumed_assert(found.count(override.first))
            <<"Override did not select an installed module: "<<override.first;
      }
      catalog=std::move(found);
    });
  }

  bool builtin(const std::string& name) const {
    const auto names=config::getBuiltInModules();
    return std::find(names.begin(),names.end(),name)!=names.end();
  }

  void closure(const std::string& name,std::set<std::string>& visiting,
               std::set<std::string>& done,std::vector<std::string>& result) {
    if(builtin(name) || done.count(name)) {
      return;
    }
    plumed_assert(visiting.insert(name).second)<<"Cyclic module dependency at "<<name;
    auto entry=catalog.find(name);
    plumed_assert(entry!=catalog.end())<<"Module "<<name<<" not found; install plumed-ext-"<<name;
    const auto& module=entry->second;
    plumed_assert(module.build==config::getModuleBuildID())
        <<"Module "<<name<<" was built for a different PLUMED base ("<<module.build<<" versus "<<config::getModuleBuildID()<<")";
    plumed_assert(fs::is_regular_file(module.library))<<"Missing module library "<<module.library.string();
    for(const auto& dep:module.dependencies) {
      closure(dep,visiting,done,result);
    }
    visiting.erase(name);
    done.insert(name);
    result.push_back(name);
  }

  void activate(const std::string& name,DLLoader& loader) {
    if(builtin(name)) {
      return;
    }
    Activation activation(activating,name);
    std::set<std::string> visiting,done;
    std::vector<std::string> dependencies;
    closure(name,visiting,done,dependencies);
    // Validate the complete closure before any constructors run. Open the root
    // first: its native links, rather than PLUMED load order, resolve symbols.
    std::reverse(dependencies.begin(),dependencies.end());
    DLLoader staged;
    std::map<std::string,void*> images;
    for(const auto& dep:dependencies) {
      auto& module=catalog.at(dep);
      auto visible=loader.getHandles();
      visible.insert(visible.end(),staged.getHandles().begin(),staged.getHandles().end());
      visible.insert(visible.end(),retained.getHandles().begin(),retained.getHandles().end());
      if(contextData) {
        visible.insert(visible.end(),contextData->images.begin(),contextData->images.end());
      }
      auto image=staged.load(pathToUTF8(fs::canonical(module.library)));
#ifdef __PLUMED_HAS_DLOPEN
      const auto symbol=moduleSymbol(dep);
      using Identity=const char*(*)();
      auto identity=reinterpret_cast<Identity>(dlsym(image,symbol.c_str()));
      plumed_assert(identity && identity()==module.build)
          <<"Missing or incorrect module identity for "<<dep<<" in "<<module.library.string();
#ifdef __PLUMED_HAS_DLADDR
      Dl_info owner;
      plumed_assert(dladdr(reinterpret_cast<const void*>(identity),&owner) && owner.dli_fname
                    && fs::equivalent(owner.dli_fname,module.library))
          <<"Module "<<dep<<" resolved a different native library than selected: "<<module.library.string();
#endif
      using Kernel=const void*(*)();
      auto kernel=reinterpret_cast<Kernel>(dlsym(image,(symbol+"_kernel").c_str()));
      plumed_assert(kernel && kernel()==reinterpret_cast<const void*>(&config::getModuleBuildID))
          <<"Module "<<dep<<" links to a different PLUMED kernel";
#endif
      for(const auto& key:actionRegister().getKeysWithDLHandle(image)) {
        if(std::find(visible.begin(),visible.end(),image)==visible.end()) {
          plumed_assert(!actionRegister().check(visible,key))
              <<"Conflicting action "<<key<<" while loading module "<<dep;
        }
      }
      for(const auto& key:cltoolRegister().getKeysWithDLHandle(image)) {
        if(std::find(visible.begin(),visible.end(),image)==visible.end()) {
          plumed_assert(!cltoolRegister().check(visible,key))
              <<"Conflicting command "<<key<<" while loading module "<<dep;
        }
      }
      staged.markModule(dep);
      images[dep]=image;
    }
#ifdef __PLUMED_HAS_DLOPEN
    for(const auto& dep:dependencies) {
      for(const auto& requirement:catalog.at(dep).dependencies) {
        if(builtin(requirement)) {
          continue;
        }
        const auto symbol=moduleSymbol(requirement);
        plumed_assert(dlsym(images.at(dep),symbol.c_str())==dlsym(images.at(requirement),symbol.c_str()))
            <<"Module "<<dep<<" resolved a different native variant of dependency "<<requirement;
      }
    }
#endif
    loader.append(std::move(staged));
  }
};
State& state() {
  static State instance;
  return instance;
}
}

ModuleLoader::Context::Context(const DLLoader& loader):
  data(std::make_unique<ModuleLoaderContextData>(loader)), previous(contextData) {
  contextData=data.get();
  Register::setContextImages(&data->images);
}

ModuleLoader::Context::~Context() {
  contextData=previous;
  Register::setContextImages(previous?&previous->images:nullptr);
}

void ModuleLoader::inheritContext(DLLoader& loader) {
  if(contextData) {
    for(const auto& name:contextData->modules) {
      load(name,loader);
    }
  }
}

void ModuleLoader::autoload(DLLoader& loader) {
  auto& shared=state();
  std::lock_guard<std::recursive_mutex> lock(shared.mutex);
  const char* enabled=std::getenv("PLUMED_AUTOLOAD");
  if(enabled && std::string(enabled)=="0") {
    return;
  }
  plumed_assert(shared.activating.empty())<<"Recursive module autoload while loading "<<shared.activating;
  shared.discover();
  std::call_once(shared.automatic,[&shared] {
    for(const auto& entry:shared.catalog) {
      plumed_assert(!shared.builtin(entry.first))
          <<"Module "<<entry.first<<" is already built into the kernel";
      shared.activate(entry.first,shared.retained);
    }
    Register::setAutomaticImages(shared.retained.getHandles());
  });
  for(const auto& entry:shared.catalog) {
    shared.activate(entry.first,loader);
  }
}

void ModuleLoader::load(const std::string& name,DLLoader& loader) {
  plumed_assert(validName(name))<<"Invalid module name "<<name;
  auto& shared=state();
  std::lock_guard<std::recursive_mutex> lock(shared.mutex);
  if(shared.builtin(name) || available(name,loader)) {
    return;
  }
  shared.discover();
  shared.activate(name,loader);
}

bool ModuleLoader::available(const std::string& name,const DLLoader& loader) {
  const auto names=modules(loader);
  return std::find(names.begin(),names.end(),name)!=names.end();
}

std::vector<std::string> ModuleLoader::modules(const DLLoader& loader) {
  auto& shared=state();
  std::lock_guard<std::recursive_mutex> lock(shared.mutex);
  auto names=config::getBuiltInModules();
  const auto& extra=loader.getModuleNames();
  names.insert(names.end(),extra.begin(),extra.end());
  std::sort(names.begin(),names.end());
  names.erase(std::unique(names.begin(),names.end()),names.end());
  return names;
}

std::vector<std::string> ModuleLoader::features(const DLLoader& loader, const std::string& module) {
  std::set<std::string> result;
  const auto add=[&result](const std::string& words) {
    std::istringstream stream(words);
    std::string feature;
    while(stream>>feature) {
      result.insert(feature);
    }
  };
  const auto builtin=config::getBuiltInModules();
  const bool inKernel=std::find(builtin.begin(),builtin.end(),module)!=builtin.end();
  if(module.empty() || inKernel) {
    std::istringstream configuration(config::getMakefile());
    std::string line;
    const std::string prefix="-D__PLUMED_" "HAS_";
    while(std::getline(configuration,line)) {
      if(line.rfind("CPPFLAGS=",0)!=0) {
        continue;
      }
      std::istringstream flags(line.substr(9));
      std::string flag;
      while(flags>>flag) {
        if(flag.rfind(prefix,0)!=0 || flag.size()<=prefix.size()) {
          continue;
        }
        auto feature=flag.substr(prefix.size());
        const auto equals=feature.find('=');
        if(equals!=std::string::npos) {
          if(feature.substr(equals+1)=="0") {
            continue;
          }
          feature.erase(equals);
        }
        std::transform(feature.begin(),feature.end(),feature.begin(),[](unsigned char c) {
          return c>='A' && c<='Z'?c-'A'+'a':c;
        });
        add(feature=="mpi"?"mpi real_mpi":feature);
        if(feature=="thread_mpi") {
          add("mpi");
        }
      }
    }
  }
  if(inKernel) {
    return {result.begin(),result.end()};
  }
#ifdef __PLUMED_HAS_DLOPEN
  for(const auto& name:loader.getModuleNames()) {
    if(!module.empty() && name!=module) {
      continue;
    }
    const auto symbol=moduleSymbol(name)+"_features";
    using Features=const char*(*)();
    for(auto image:loader.getHandles()) {
      auto features=reinterpret_cast<Features>(dlsym(image,symbol.c_str()));
      if(features) {
        add(features());
        break;
      }
    }
  }
#endif
  return {result.begin(),result.end()};
}

}
