/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2023 The plumed team
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
#include "CLToolMain.h"
#include "config/Config.h"
#include "tools/Exception.h"
#include "tools/Communicator.h"
#include "CLTool.h"
#include "CLToolRegister.h"
#include "tools/Tools.h"
#include "tools/DLLoader.h"
#include "ModuleLoader.h"
#include "tools/Subprocess.h"
#include "../mpiwrap/mpiwrap.h"
#include <string>
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <climits>
#include <cerrno>
#include <iostream>
#include <algorithm>
#include <memory>

namespace PLMD {

namespace {

#if defined(__PLUMED_MPIWRAP) && defined(__PLUMED_HAS_THREAD_MPI) && !defined(__PLUMED_HAS_MPI)
struct ThreadMPICommand {
  std::vector<std::string> arguments;
  FILE *input;
  FILE *output;
};

int runThreadMPICommand(int, const void *native_world, std::size_t native_size,
                        void *opaque) {
  if(!native_world || native_size==0) {
    throw ExceptionError("internal launcher did not supply the native world communicator");
  }
  const ThreadMPICommand &context=*static_cast<ThreadMPICommand *>(opaque);
  // The outer CLToolMain is only the dispatcher. Construct a fresh rank-local
  // command object after tMPI initialization so it gets the real rank communicator.
  CLToolMain command;
  command.cmd("setMPIComm",TypesafePtr::unchecked(native_world));
  int argc=static_cast<int>(context.arguments.size());
  command.cmd("setArgc",&argc);
  command.cmd("setIn",context.input);
  command.cmd("setOut",context.output);
  std::size_t bytes=0;
  for(const auto &arg:context.arguments) {
    bytes+=arg.size()+1;
  }
  std::vector<char> storage(bytes);
  std::vector<char *> argv(argc);
  char *cursor=storage.data();
  for(int i=0; i<argc; ++i) {
    argv[i]=cursor;
    const std::string &arg=context.arguments[static_cast<std::size_t>(i)];
    std::copy(arg.begin(),arg.end(),cursor);
    cursor+=arg.size();
    *cursor++='\0';
  }
  command.cmd("setArgv",argv.data());
  int result=0;
  command.cmd("run",&result);
  return result;
}
#endif

bool parseThreadMPIRanks(const char *text,int &ranks) {
  if(!text || !text[0]) {
    return false;
  }
  errno=0;
  char *end=nullptr;
  const long value=std::strtol(text,&end,10);
  if(errno || end==text || *end!='\0' || value<1 || value>INT_MAX) {
    return false;
  }
  ranks=static_cast<int>(value);
  return true;
}

} // namespace

CLToolMain::CLToolMain():
  argc(0),
  in(stdin),
  out(stdout) {
}

CLToolMain::~CLToolMain() {
// empty destructor to delete unique_ptr
}

#define CHECK_NULL(val,word) plumed_assert(val) << "NULL pointer received in cmd(\"CLTool " << word <<  "\"";

void CLToolMain::cmd(std::string_view word,const TypesafePtr & val) {

// Enumerate all possible commands:
  enum {
#include "CLToolMainEnum.inc"
  };

// Static object (initialized once) containing the map of commands:
  const static Tools::FastStringUnorderedMap<int> word_map = {
#include "CLToolMainMap.inc"
  };

  gch::small_vector<std::string_view> words;
  Tools::getWordsSimple(words,word);
  unsigned nw=words.size();
  if(nw==0) {
    // do nothing
  } else {
    int iword=-1;
    const char*const*v;
    const char*vv;
    const auto it=word_map.find(words[0]);
    if(it!=word_map.end()) {
      iword=it->second;
    }
    switch(iword) {
    case cmd_setArgc:
      CHECK_NULL(val,word);
      argc=val.get<int>();
      break;
    case cmd_setArgv:
      CHECK_NULL(val,word);
      v=val.get<const char*const*>({argc});
      for(int i=0; i<argc; ++i) {
        argv.push_back(std::string(v[i]));
      }
      break;
    case cmd_setArgvLine:
      CHECK_NULL(val,word);
      vv=val.getCString();
      argv=Tools::getWords(vv);
      break;
    case cmd_setIn:
      CHECK_NULL(val,word);
      in=val.get<FILE*>();
      break;
    case cmd_setOut:
      CHECK_NULL(val,word);
      out=val.get<FILE*>();
      break;
    case cmd_setMPIComm:
      comm.Set_comm(val);
      break;
    case cmd_setMPIFComm:
      comm.Set_fcomm(val);
      break;
    case cmd_run:
      CHECK_NULL(val,word);
      argc=argv.size();
      {
        int n=0;
        for(int i=0; i<argc; ++i) {
          n+=argv[i].length()+1;
        }
        std::vector<char> args(n);
        std::vector<char*> vvv(argc);
        char* ptr=args.data();
        for(int i=0; i<argc; ++i) {
          vvv[i]=ptr;
          for(unsigned c=0; c<argv[i].length(); ++c) {
            *ptr=argv[i][c];
            ptr++;
          }
          *ptr=0;
          ptr++;
        }
        val.set(int(run(argc,vvv.data(),in,out,comm)));
      }
      break;
    default:
      plumed_error() << "cannot interpret cmd(\"CLTool " << word << "\"). check plumed developers manual to see the available commands.";
      break;
    }
  }
}

/**
This is the entry point to the command line tools
included in the plumed library.
*/

int CLToolMain::run(int argc, char **argv,FILE*in,FILE*out,Communicator& pc) {
  int i;
  bool printhelp=false;

  int thread_mpi_ranks=0;
  bool thread_mpi_requested=false;
  bool invalid_thread_mpi_option=false;
  bool explicit_mpi=false;
  bool explicit_no_mpi=false;
  std::vector<std::string> filtered_arguments;
  if(argc>0) {
    filtered_arguments.emplace_back(argv[0]);
  }
  i=1;
  while(i<argc) {
    if(std::strcmp(argv[i],"--tmpi-threads")==0) {
      int parsed=0;
      if(thread_mpi_requested || i+1>=argc ||
          !parseThreadMPIRanks(i+1<argc ? argv[i+1] : nullptr,parsed)) {
        invalid_thread_mpi_option=true;
      } else {
        thread_mpi_requested=true;
        thread_mpi_ranks=parsed;
      }
      if(i+1<argc) {
        ++i;
      }
    } else if(std::strncmp(argv[i],"--tmpi-threads=",15)==0) {
      int parsed=0;
      if(thread_mpi_requested || !parseThreadMPIRanks(argv[i]+15,parsed)) {
        invalid_thread_mpi_option=true;
      } else {
        thread_mpi_requested=true;
        thread_mpi_ranks=parsed;
      }
    } else {
      explicit_mpi=explicit_mpi || std::strcmp(argv[i],"--mpi")==0;
      explicit_no_mpi=explicit_no_mpi || std::strcmp(argv[i],"--no-mpi")==0;
      filtered_arguments.emplace_back(argv[i]);
    }
    ++i;
  }
  if(invalid_thread_mpi_option) {
    std::fprintf(stderr,"ERROR: --tmpi-threads requires one positive rank count and may be used once\n");
    return 1;
  }
  if(thread_mpi_requested) {
    if(explicit_mpi || explicit_no_mpi) {
      std::fprintf(stderr,"ERROR: --tmpi-threads cannot be combined with --mpi or --no-mpi\n");
      return 1;
    }
#if defined(__PLUMED_MPIWRAP) && defined(__PLUMED_HAS_THREAD_MPI) && !defined(__PLUMED_HAS_MPI)
    ThreadMPICommand context {filtered_arguments,in,out};
    int result=0;
    mpiwrap::launchRanks(thread_mpi_ranks,
                         &runThreadMPICommand,&context,&result,
                         mpiwrap::MPILifecyclePolicy::internal_thread_mpi);
    return result;
#else
    std::fprintf(stderr,"ERROR: --tmpi-threads requires an embedded thread-MPI build (requested %d ranks)\n",thread_mpi_ranks);
    return 1;
#endif
  }

  DLLoader dlloader;
  ModuleLoader::autoload(dlloader);

  std::string root=config::getPlumedRoot();

  bool standalone_executable=false;

// Start parsing options
  std::string prefix("");
  std::string a("");
  for(i=1; i<argc; i++) {
    a=prefix+argv[i];
    if(a.length()==0) {
      continue;
    }
#if defined(__PLUMED_MPIWRAP) && !defined(__PLUMED_HAS_MPI)
    // Explicit CLI inquiries honor injection or an environment-selected provider
    // without initializing MPI. Embedded thread-MPI needs an explicit launch.
    if(a=="--has-mpi" || a=="--has-real-mpi" ||
        a=="--has-thread-mpi" || a=="--mpi-ranks-are-threads") {
      mpiwrap::ensure_provider_selected();
    }
#endif
    if(a=="help" || a=="-h" || a=="--help") {
      printhelp=true;
      break;
    } else if(a=="--has-mpi") {
      return Communicator::plumedHasMPI()?0:1;
    } else if(a=="--has-real-mpi") {
      return mpiwrap::provider_installed() && !mpiwrap::provider_ranks_are_threads()?0:1;
    } else if(a=="--has-thread-mpi") {
#if defined(__PLUMED_HAS_THREAD_MPI)
      return 0;
#else
      return mpiwrap::provider_ranks_are_threads()?0:1;
#endif
    } else if(a=="--mpi-ranks-are-threads") {
      return mpiwrap::provider_ranks_are_threads()?0:1;
    } else if(a=="--has-cregex") {
      return (config::hasCregex()?0:1);
    } else if(a=="--has-dlopen") {
      return (config::hasDlopen()?0:1);
    } else if(a=="--has-molfile") {
      return (config::hasMolfile()?0:1);
    } else if(a=="--has-external-molfile") {
      return (config::hasExternalMolfile()?0:1);
    } else if(a=="--has-zlib") {
      return (config::hasZlib()?0:1);
    } else if(a=="--has-xdrfile") {
      return 0; // always ok
    } else if(a=="--is-installed") {
      return (config::isInstalled()?0:1);
    } else if(a=="--no-mpi") {
// this is ignored, as it is parsed in main
      continue;
    } else if(a=="--mpi") {
// this is ignored, as it is parsed in main
      continue;
    } else if(a=="--standalone-executable") {
      standalone_executable=true;
    } else if(Tools::startWith(a,"--load-module=")) {
      ModuleLoader::load(a.substr(a.find('=')+1),dlloader);
      prefix="";
    } else if(a=="--load-module") {
      prefix="--load-module=";
    } else if(a=="--module-build-id") {
      std::fprintf(out,"%s\n",config::getModuleBuildID().c_str());
      return 0;
    } else if(a=="--list-modules") {
      for(const auto& name:ModuleLoader::modules(dlloader)) {
        std::fprintf(out,"%s\n",name.c_str());
      }
      return 0;
    } else if(Tools::startWith(a,"--load=")) {
      a.erase(0,a.find("=")+1);
      prefix="";
      dlloader.load(a);
    } else if(a=="--load") {
      prefix="--load=";
    } else if(a[0]=='-') {
      std::string msg="ERROR: Unknown option " +a;
      std::fprintf(stderr,"%s\n",msg.c_str());
      return 1;
    } else {
      break;
    }
  }

  if(!prefix.empty()) {
    std::fprintf(stderr,"ERROR: missing argument for %s\n",prefix.c_str());
    return 1;
  }

  ModuleLoader::Context moduleContext(dlloader);

// Check if plumedRoot/patches/ directory exists (as a further check)
  if(!standalone_executable) {
    std::vector<std::string> files=Tools::ls(root);
    if(find(files.begin(),files.end(),"patches")==files.end()) {
      std::string msg=
        "WARNING: I cannot find "+root+"/patches/ directory. Set PLUMED_ROOT or reinstall PLUMED\n\n";
      std::fprintf(stderr,"%s",msg.c_str());
    }
  }

// Build list of available C++ tools:
  std::vector<std::string> availableCxx=cltoolRegister().list();
  for(auto image:dlloader.getHandles()) {
    const auto extra=cltoolRegister().getKeysWithDLHandle(image);
    availableCxx.insert(availableCxx.end(),extra.begin(),extra.end());
  }
  std::sort(availableCxx.begin(),availableCxx.end());
  availableCxx.erase(std::unique(availableCxx.begin(),availableCxx.end()),availableCxx.end());
// Build list of available shell tools:
  std::vector<std::string> availableShell;
  if(!standalone_executable) {
    std::vector<std::string> tmp;
    tmp=Tools::ls(std::string(root+"/scripts"));
    for(unsigned j=0; j<tmp.size(); ++j) {
      size_t ff=tmp[j].find(".sh");
      if(ff==std::string::npos) {
        tmp[j].erase();
      } else {
        tmp[j].erase(ff);
      }
    }
    for(unsigned j=0; j<tmp.size(); ++j)
      if(tmp[j].length()>0) {
        availableShell.push_back(tmp[j]);
      }
  }

  if(printhelp) {
    std::string msg=
      "Usage: plumed [options] [command] [command options]\n"
      "  plumed [command] -h|--help: to print help for a specific command\n"
      "Options:\n"
      "  [help|-h|--help]          : to print this help\n"
      "  [--is-installed]          : fails if plumed is not installed\n"
      "  [--has-mpi]               : fails if neither real nor internal thread MPI is available\n"
      "  [--has-real-mpi]          : fails if no real MPI provider is available\n"
      "  [--has-thread-mpi]        : fails if no thread MPI capability is available\n"
      "  [--mpi-ranks-are-threads] : succeeds if the selected MPI provider uses threads\n"
      "  [--has-dlopen]            : fails if plumed is compiled without dlopen\n"
      "  [--load-module MODULE]   : loads an installed module and its hard dependencies\n"
      "  [--load LIB]              : loads a shared object (typically a plugin library)\n"
      "  [--standalone-executable] : tells plumed not to look for commands implemented as scripts\n"
      "Commands:\n";
    std::fprintf(out,"%s",msg.c_str());
    for(unsigned j=0; j<availableCxx.size(); ++j) {
      auto cl=cltoolRegister().create(dlloader.getHandles(),CLToolOptions(availableCxx[j]));
      plumed_assert(cl);
      std::string manual=availableCxx[j]+" : "+cl->description();
      std::fprintf(out,"  plumed %s\n", manual.c_str());
    }
    for(unsigned j=0; j<availableShell.size(); ++j) {
      std::string manual;
      if(Subprocess::available()) {
        auto cmd=config::getEnvCommand()+" \""+root+"/scripts/"+availableShell[j]+".sh\" --description";
        auto proc=Subprocess(cmd);
        // we only need the first line, so this is fine:
        proc.getline(manual);
        // process will be killed immediately after
      } else {
        manual="(doc not avail)";
      }
      manual= availableShell[j]+" : "+manual;
      std::fprintf(out,"  plumed %s\n", manual.c_str());
    }
    return 0;
  }
  if(i==argc) {
    std::fprintf(out,"%s","Nothing to do. Use 'plumed help' for help\n");
    return 0;
  }

// this is the command to be executed:
  std::string command(argv[i]);
  // Module availability is runtime state; other config inquiries retain their
  // build-provenance semantics in config.sh.
  if(command=="config") {
    bool query=false,quiet=false;
    std::vector<std::string> names;
    for(int j=i+1; j<argc; ++j) {
      const std::string option=argv[j];
      if(option=="module") {
        query=true;
      } else if(option=="-q" || option=="--quiet") {
        quiet=true;
      } else {
        names.push_back(option);
      }
    }
    if(query) {
      if(names.empty()) {
        return 1;
      }
      int result=0;
      for(const auto& name:names) {
        const bool present=ModuleLoader::available(name,dlloader);
        if(!quiet) {
          std::fprintf(out,"%s %s\n",name.c_str(),present?"on":"off");
        }
        if(!present) {
          result=1;
        }
      }
      return result;
    }
  }

  if(find(availableCxx.begin(),availableCxx.end(),command)!=availableCxx.end()) {
    auto cl=cltoolRegister().create(dlloader.getHandles(),CLToolOptions(command));
    plumed_assert(cl);
    cl->setDLHandles(dlloader.getHandles());
    int ret=0;
    // NOTE: catch PLUMED exceptions from CLTool readInput() and main() so that user-input
    // errors terminate cleanly with a message and nonzero exit status
    // instead of an uncaught-exception abort (0xc0000409 on Windows).
    try {
      // Read the command line options (returns false if we are just printing help)
      if( !cl->readInput( argc-i,&argv[i],in,out ) ) {
        return 0;
      }
      ret=cl->main(in,out,pc);
    } catch(const PLMD::Exception& e) {
      std::fprintf(stderr,"ERROR : in input for command line tool %s : %s\n",command.c_str(),e.what());
      std::fflush(stderr);
      return 1;
    } catch(const std::exception& e) {
      std::fprintf(stderr,"ERROR : in command line tool %s : %s\n",command.c_str(),e.what());
      std::fflush(stderr);
      return 1;
    }
    return ret;
  }

  if(find(availableShell.begin(),availableShell.end(),command)!=availableShell.end()) {
    plumed_massert(in==stdin,"shell tools can only work on stdin");
    plumed_massert(out==stdout,"shell tools can only work on stdin");
    std::string cmd=config::getEnvCommand()+" \""+root+"/scripts/"+command+".sh\"";
    for(int j=i+1; j<argc; j++) {
      // mklib accepts literal paths, including spaces and shell metacharacters.
      // Preserve the legacy argument handling of the other shell tools here.
      cmd+=" "+(command=="mklib" ? Tools::quoteForShell(argv[j]) : std::string(argv[j]));
    }
    int r=std::system(cmd.c_str());
// this is necessary since system seems to return numbers which are multiple
// of 256. this would make the interpretation by the shell wrong
// I just return 1 in case of failure and 0 in case of success
    if(r!=0) {
      return 1;
    } else {
      return 0;
    }
  }

  std::string msg="ERROR: unknown command " + command + ". Use 'plumed help' for help";
  std::fprintf(stderr,"%s\n",msg.c_str());
  return 1;

}
}
