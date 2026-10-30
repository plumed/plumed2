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
#include "MetaDGrid.h"
#include "MetaDHills.h"
#include "core/Action.h"
#include "tools/Communicator.h"
#include "tools/File.h"
#include "tools/Keywords.h"
#include "tools/Log.h"
#include <exception>
#include <filesystem>
#include <limits>
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

std::vector<long long> readBiasGrid(Action& action, const MetaDWalkers& walkers, MetaDBiasGrid& grid, const std::string& filename,
                                    const std::string& funcl, const std::vector<Value*>& args, const std::vector<std::string>& gmin,
                                    const std::vector<std::string>& gmax, const std::vector<unsigned>& gbin, bool sparse, bool spline,
                                    bool nobounds) {
  std::string path=filename;
  if(walkers.mpi && path.at(0)!='/') {
    path=walkers.firstReplicaDirectory(action)+"/"+path;
  }
  IFile gridfile;
  gridfile.link(action);
  if(walkers.shareFiles()) {
    // multiple walkers write their files without the suffix of the replica: with WALKERS_MPI all the walkers
    // read the same file
    gridfile.enforceSuffix("");
  }
  if(gridfile.FileExist(path)) {
    gridfile.open(path);
  } else {
    action.error("The GRID file you want to read: " + path + ", cannot be found!");
  }
  grid.read(gridfile,path,funcl,args,gmin,gmax,gbin,sparse,spline,nobounds);
  std::vector<long long> walkerhills;
  if(walkers.n>1) {
    // the hills of the other walkers in the grid
    walkerhills=MetaDHillsFiles::readGridComment(action.comm,gridfile.getPath(),walkers.n);
  }
  action.log.printf("  Restarting from %s\n",path.c_str());
  return walkerhills;
}

void readHillsAtStart(Action& action, MetaDWalkers& walkers, const std::vector<MetaDHillsInput>& inputs, bool restartedFromGrid,
                      MetaDHillSchedule& schedule) {
  Log& log=action.log;
  Communicator& comm=action.comm;
  const bool restart=action.getRestart();
  // with WALKERS_MPI the walkers must all read a grid, or none
  walkers.checkSameChoice(action,restartedFromGrid,"GRID");

  if(walkers.mpi && walkers.dir=="") {
    // if possible the first replica shares its current directory, so that all the walkers read the same files
    bool relative=true;
    for(const auto & input : inputs) {
      relative=relative && input.name.at(0)!='/';
    }
    if(relative) {
      walkers.dir=walkers.firstReplicaDirectory(action);
    }
  }

  // open all the files
  for(const auto & input : inputs) {
    for(int j=0; j<walkers.n; ++j) {
      // with multiple walkers the files are written without the suffix of the replica
      const unsigned k=input.files->addInput(walkers.getFileName(input.name,j),walkers.shareFiles());
      plumed_assert(k==static_cast<unsigned>(j));
    }
  }

  // When the bias is read from grids, the own hills files are read first: the time of their last hill tells
  // which hills of the other walkers (WALKERS_N) are already in the grids (for grids without their numbers)
  std::vector<int> order;
  if(restart && restartedFromGrid) {
    order.push_back(walkers.id);
  }
  for(int j=0; j<walkers.n; ++j) {
    if(order.empty() || order[0]!=j) {
      order.push_back(j);
    }
  }
  bool restartedFromHills=false;
  // the time of the last hill in the files of this walker (the hills of all the files are written at the same time:
  // the last one is used)
  bool lastHillTimeKnown=false;
  double lastHillTime=0.0;
  // (lasttime is set by readHills, and is valid only if at least one hill was read)
  auto setLastHill=[&](unsigned nhills, double lasttime) {
    if(nhills>0 && !(lastHillTimeKnown && lasttime<=lastHillTime)) {
      lastHillTimeKnown=true;
      lastHillTime=lasttime;
    }
  };
  // the actual paths of the hills files of this walker, if they were read
  std::vector<std::string> ownpath(inputs.size());
  for(const int j : order) {
    for(unsigned s=0; s<inputs.size(); ++s) {
      MetaDHillsFiles& files=*inputs[s].files;
      const unsigned k=j;
      if(!files.isOpen(k)) {
        // in case a file does not exist and we are restarting, complain that the file was not found
        if(restart && !restartedFromGrid) {
          action.error("restart file "+files.getInputName(k)+" not found");
        }
        continue;
      }
      if(restart && !restartedFromGrid) {
        log.printf("  Restarting from %s:",files.getInputName(k).c_str());
        if(j==walkers.id) {
          double lasttime=0.0;
          const unsigned nhills=files.readHills(k,&lasttime);
          setLastHill(nhills,lasttime);
          ownpath[s]=files.getPath(k);
        } else {
          files.readHills(k);
        }
        restartedFromHills=true;
      }
      if(restart && restartedFromGrid && j==walkers.id) {
        // new hills will be appended to this file also when the bias is read from a grid
        ownpath[s]=files.getPath(k);
        // the bias is read from the grid, but the hills are still needed for the time of the last hill (to decide
        // if a hill is added at the first step, and which hills of the other walkers are in the grid) and for the
        // statistics of ADAPTIVE=DIFF
        log.printf("  Reading hills from %s:",files.getInputName(k).c_str());
        double lasttime=0.0;
        try {
          const unsigned nhills=files.readHills(k,&lasttime,false);
          setLastHill(nhills,lasttime);
        } catch(const std::exception& e) {
          // the bias comes from the grid: if the hills file cannot be used, its information is not used
          log.printf("\n  WARNING: the hills file cannot be used for the time of the last hill and the statistics of ADAPTIVE=DIFF"
                     " (new hills are appended to it):\n%s\n",e.what());
          if(files.getBias().diffrestart) {
            files.getBias().diffrestart->clear();
          }
        }
      }
      if(restart && restartedFromGrid && j!=walkers.id) {
        // the hills of the other walkers that are not in the grid are added
        const std::vector<long long>& ingridcount=inputs[s].gridWalkerHills;
        const double ingrid=(lastHillTimeKnown ? lastHillTime : std::numeric_limits<double>::infinity());
        if(ingridcount.size()==static_cast<unsigned>(walkers.n) && ingridcount[j]>=0) {
          // the grid file contains the number of hills read from each walker
          log.printf("  Reading hills from %s (the first %lld are in the grid):",files.getInputName(k).c_str(),ingridcount[j]);
          files.readHills(k,nullptr,true,-std::numeric_limits<double>::infinity(),ingridcount[j]);
        } else if(lastHillTimeKnown) {
          // grid written by an older version: those up to the last hill of this walker are assumed to be in the grid
          log.printf("  Reading hills from %s (those up to time %f are assumed to be in the grid):",files.getInputName(k).c_str(),ingrid);
          files.readHills(k,nullptr,true,ingrid);
        } else {
          log.printf("  Reading hills from %s (they are assumed to be in the grid):",files.getInputName(k).c_str());
          files.readHills(k,nullptr,true,ingrid);
        }
      }
      files.reset(k);
      // close only the walker own hills file for later writing
      if(j==walkers.id) {
        files.close(k);
      }
    }
  }

  // with WALKERS_MPI the walkers must all read the hills files, or none
  walkers.checkSameChoice(action,restartedFromHills,"FILE");

  if(walkers.mpi) {
    // the decision to add a hill at the first step must be the same for all the walkers (they share the hills)
    std::vector<double> last {lastHillTimeKnown ? 1.0 : 0.0, lastHillTime};
    if(comm.Get_rank()==0) {
      action.multi_sim_comm.Bcast(last,0);
    }
    comm.Bcast(last,0);
    lastHillTimeKnown=(last[0]!=0.0);
    lastHillTime=last[1];
  }
  if(lastHillTimeKnown) {
    schedule.setLastHillTime(lastHillTime);
  }

  comm.Barrier();
  // this barrier is needed when using walkers_mpi
  // to be sure that all files have been read before
  // backing them up
  // it should not be used when walkers_mpi is false otherwise
  // it would introduce troubles when using replicas without METAD
  // (e.g. in bias exchange with a neutral replica)
  // see issue #168 on github
  if(comm.Get_rank()==0 && walkers.mpi) {
    action.multi_sim_comm.Barrier();
  }

  // all walkers have read the hills: the one that writes the files can now remove an incomplete
  // last line (left by a run that was interrupted while writing), before new hills are appended to them
  for(const auto & path : ownpath) {
    std::string message;
    if(path.length()>0 && comm.Get_rank()==0 && walkers.writesFiles()) {
      message=removeIncompleteLastLine(path);
    }
    // the log must be written by all the processes
    if(comm.Get_size()>1) {
      bcastString(comm,message);
    }
    if(message.length()>0) {
      log<<message;
    }
  }
}

void readOtherWalkersHills(Action& action, const MetaDWalkers& walkers, const std::vector<MetaDHillsFiles*>& files) {
  for(int j=0; j<walkers.n; ++j) {
    // don't read your own Gaussians
    if(j==walkers.id) {
      continue;
    }
    for(auto f : files) {
      if(!f->isOpen(j)) {
        // check if it exists now and open it!
        if(f->openIfExists(j)) {
          f->reset(j);
        }
      } else {
        // otherwise read the new Gaussians
        action.log.printf("  Reading hills from %s:",f->getInputName(j).c_str());
        f->readHills(j);
        f->reset(j);
      }
    }
  }
}

}
}
