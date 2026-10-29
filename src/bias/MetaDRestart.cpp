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
#include "MetaDRestart.h"
#include "MetaDHills.h"
#include "core/Action.h"
#include "tools/Communicator.h"
#include "tools/File.h"
#include "tools/Keywords.h"
#include "tools/Log.h"
#include <filesystem>
#include <numeric>
#include <vector>

namespace PLMD {
namespace bias {

void MetaDWalkers::registerKeywords(Keywords& keys) {
  keys.add("optional","WALKERS_ID", "walker id");
  keys.add("optional","WALKERS_N", "number of walkers");
  keys.add("optional","WALKERS_DIR", "shared directory with the hills files from all the walkers");
  keys.add("optional","WALKERS_RSTRIDE","stride for reading hills files");
  keys.addFlag("WALKERS_MPI",false,"Switch on MPI version of multiple walkers - not compatible with WALKERS_* options other than WALKERS_DIR");
}

void MetaDWalkers::parse(Action& action) {
  action.parse("WALKERS_N",n);
  action.parse("WALKERS_ID",id);
  if(n<1) {
    action.error("the number of walkers should be positive");
  }
  if(id<0 || n<=id) {
    action.error("walker ID should be a numerical value between 0 and the total number of walkers minus one");
  }
  action.parse("WALKERS_DIR",dir);
  action.parse("WALKERS_RSTRIDE",rstride);
  if(rstride<=0) {
    action.error("WALKERS_RSTRIDE should be positive");
  }
  action.parseFlag("WALKERS_MPI",mpi);
  // if PLUMED is not compiled with MPI the user is informed and we exit gracefully
  if(mpi) {
    if(!Communicator::plumedHasMPI()) {
      action.error("Invalid walkers configuration: WALKERS_MPI flag requires MPI compilation");
    }
    if(!Communicator::initialized()) {
      action.error("Invalid walkers configuration: WALKERS_MPI needs the communicator correctly initialized.");
    }
  }
}

void MetaDWalkers::setup(Action& action) {
  Log& log=action.log;
  if(n>1) {
    if(mpi) {
      action.error("MPI version of multiple walkers is not compatible with filesystem version of multiple walkers");
    }
    log.printf("  %d multiple walkers active\n",n);
    log.printf("  walker id %d\n",id);
    log.printf("  reading stride %d\n",rstride);
    if(dir!="") {
      log.printf("  directory with hills files %s\n",dir.c_str());
    }
  } else if(mpi) {
    log.printf("  Multiple walkers active using MPI communication\n");
    if(dir!="") {
      log.printf("  directory with hills files %s\n",dir.c_str());
    }
    mpi_nw=action.getNumberOfReplicas();
    mpi_id=action.getReplicaIndex();
  }
}

std::string MetaDWalkers::getFileName(const std::string& name, int walker) const {
  return MetaDHillsFiles::getWalkerFileName(name,dir,n,walker,mpi);
}

std::string MetaDWalkers::firstReplicaDirectory(Action& action) const {
  std::string cwd=std::filesystem::current_path().string();
  // the paths of the replicas can have different lengths
  if(action.comm.Get_rank()==0) {
    bcastString(action.multi_sim_comm,cwd);
  }
  bcastString(action.comm,cwd);
  return cwd;
}

void MetaDWalkers::checkSameChoice(Action& action, bool choice, const std::string& what) const {
  if(!action.getRestart() || !mpi) {
    return;
  }
  std::vector<int> all(mpi_nw,0);
  if(action.comm.Get_rank()==0) {
    action.multi_sim_comm.Allgather(int(choice),all);
  }
  action.comm.Bcast(all,0);
  const unsigned nchoice=std::accumulate(all.begin(),all.end(),0);
  if(nchoice!=0 && nchoice!=mpi_nw) {
    action.error("in this WALKERS_MPI run some replicas have restarted from "+what+" while others have not");
  }
}

void MetaDWalkers::openOutput(Action& action, OFile& file, const std::string& name) const {
  file.link(action);
  if(shareFiles()) {
    file.enforceSuffix("");
  }
  file.open(writesFiles() ? name : "/dev/null");
}

}
}
