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
#include "Checkpointing.h"
#include "PlumedMain.h"
#include "ActionSet.h"
#include "Action.h"
#include "tools/Checkpoint.h"
#include "tools/Communicator.h"
#include "tools/Exception.h"
#include "tools/OFile.h"
#include "tools/Tools.h"
#include "config/Config.h"
#include <algorithm>
#include <filesystem>
#include <fstream>
#include <limits>

namespace PLMD {

namespace {
// identity check of a trimmed file: checksum of its last bytes
const unsigned long long tailBytes=4096;
}

Checkpointing::Checkpointing(PlumedMain& p):
  plumed(p) {
}

Checkpointing::~Checkpointing() = default;

std::string Checkpointing::suffixedPath(const std::string& fileName) const {
  return FileBase::appendSuffix(fileName,plumed.getSuffix());
}

void Checkpointing::setup(const std::string& fileName,long long int checkpointStride) {
  plumed_massert(!enabled(),"only one CHECKPOINT action is allowed");
  plumed_massert(!fileName.empty(),"CHECKPOINT: FILE cannot be empty");
  plumed_massert(outputFiles.empty(),"CHECKPOINT must come before any action that writes a file");
  plumed_massert(Tools::extension(fileName)!="gz","CHECKPOINT: the checkpoint file cannot be compressed");
  file=fileName;
  stride=checkpointStride;
  path=suffixedPath(file);
  // plumed.cpt -> plumed_prev.cpt (the suffix, if any, stays at the end)
  std::string prev=file;
  const size_t slash=prev.find_last_of("/\\");
  const size_t dot=prev.find_last_of('.');
  if(dot!=std::string::npos && (slash==std::string::npos || dot>slash)) {
    prev.insert(dot,"_prev");
  } else {
    prev+="_prev";
  }
  prevPath=suffixedPath(prev);
  // a new run does not overwrite the checkpoints of a previous one
  if(!pending) {
    OFile f;
    f.link(plumed);
    f.backupFile("bck",path);
    f.backupFile("bck",prevPath);
  }
}

void Checkpointing::registerOutputFile(OFile& f) {
  const std::string fileName=f.getPath();
  if(f.getMode()=="a") {
    const auto it=lastFields.find(fileName);
    if(it!=lastFields.end()) {
      f.setLastFields(it->second);
    }
  }
  if(std::find(outputFiles.begin(),outputFiles.end(),fileName)!=outputFiles.end()) {
    return;
  }
  // devices and pipes (e.g. /dev/null, written by the processes or replicas that do not write a shared file)
  // have no size to record
  std::error_code ec;
  if(std::filesystem::exists(fileName,ec) && !std::filesystem::is_regular_file(fileName,ec)) {
    return;
  }
  if(enabled() && Tools::extension(fileName)=="gz") {
    plumed_merror("file "+fileName+": compressed output files cannot be trimmed at restart, remove the .gz extension or the CHECKPOINT action");
  }
  outputFiles.push_back(fileName);
}

void Checkpointing::markSnapshotFile(const std::string& fileName) {
  snapshotFiles.insert(fileName);
}

std::string Checkpointing::crcOfTail(const std::string& fname,unsigned long long size) {
  std::ifstream in(fname,std::ios::binary);
  if(!in) {
    plumed_merror("checkpoint: cannot open "+fname+" for reading");
  }
  const unsigned long long n=std::min(size,tailBytes);
  std::string tail(n,'\0');
  in.seekg(size-n);
  in.read(&tail[0],n);
  if(static_cast<unsigned long long>(in.gcount())!=n) {
    plumed_merror("checkpoint: cannot read the last "+std::to_string(n)+" bytes of "+fname);
  }
  return std::to_string(CheckpointReader::checksum(tail));
}

void Checkpointing::write() {
  plumed_assert(enabled());
  auto& log=plumed.getLog();
  auto& comm=plumed.comm;
  // everything written so far by step-1 must be on disk before the sizes are measured
  plumed.fflush();
  log.flush();
  for(const auto & p : plumed.getActionSet()) {
    p->fflush();
  }
  CheckpointWriter w(config::getVersionLong());
  w.write("step",plumed.getStep());
  w.write("replica",plumed.getReplicaIndex());
  w.write("nreplicas",plumed.getNumberOfReplicas());

  w.beginSection("files");
  w.write("n",outputFiles.size());
  if(comm.Get_rank()==0) {
    for(unsigned i=0; i<outputFiles.size(); ++i) {
      const std::string& f=outputFiles[i];
      const std::string num=std::to_string(i);
      std::error_code ec;
      unsigned long long size=0;
      if(std::filesystem::exists(f,ec)) {
        size=std::filesystem::file_size(f,ec);
        if(ec) {
          plumed_merror("checkpoint: cannot get the size of "+f+": "+ec.message());
        }
      }
      w.write("path"+num,f);
      w.write("size"+num,size);
      w.write("crc"+num,size>0 ? crcOfTail(f,size) : std::string("0"));
      w.write("snapshot"+num,snapshotFiles.count(f)>0);
      // header of the last line written, so that the restart does not print it again
      for(const auto & fb : plumed.getFiles()) {
        OFile* of=dynamic_cast<OFile*>(fb);
        if(of && of->isOpen() && of->getPath()==f) {
          w.write("fields"+num,of->getLastFields());
          break;
        }
      }
    }
  }
  w.endSection();

  w.beginSection("actions");
  std::vector<std::string> unported;
  for(const auto & p : plumed.getActionSet()) {
    if(!p->hasCheckpointState()) {
      continue;
    }
    w.beginSection(p->getLabel());
    p->writeCheckpoint(w);
    if(w.getNumberOfEntriesInSection()==0) {
      unported.push_back(p->getLabel());
    }
    w.endSection();
  }
  w.endSection();

  if(!unported.empty() && !warnedAboutUnportedState) {
    warnedAboutUnportedState=true;
    log<<"WARNING: the following actions have an internal state that is not stored in the checkpoint, an exact restart is not guaranteed:";
    for(const auto & l : unported) {
      log<<" "<<l;
    }
    log<<"\n";
  }

  if(comm.Get_rank()==0) {
    w.writeFile(path,prevPath);
  }
  comm.Barrier();
  lastWrittenStep=plumed.getStep();
  log<<"Checkpoint written to "<<path<<" at step "<<plumed.getStep()<<"\n";
}

void Checkpointing::restartFrom(const std::string& fileName) {
  auto& log=plumed.getLog();
  auto& comm=plumed.comm;
  plumed_massert(!pending,"only one CHECKPOINT action is allowed");
  plumed_massert(outputFiles.empty(),"CHECKPOINT must come before any action that writes a file");
  const std::string fname=suffixedPath(fileName);
  // rank 0 reads, the content is broadcast so that every rank can parse it
  std::string content;
  if(comm.Get_rank()==0) {
    content=CheckpointReader::readFile(fname);
  }
  if(comm.Get_size()>1) {
    unsigned long long n=content.size();
    comm.Bcast(n,0);
    plumed_massert(n<static_cast<unsigned long long>(std::numeric_limits<int>::max()),"checkpoint file too large");
    content.resize(n);
    comm.Bcast(&content[0],static_cast<int>(n),0);
  }
  pending=std::make_unique<CheckpointReader>(content);
  restartStep=pending->get<long long int>("step");
  const int nrep=pending->get<int>("nreplicas");
  if(nrep!=plumed.getNumberOfReplicas()) {
    plumed_merror("checkpoint "+fname+" was written by a run with "+std::to_string(nrep)+" replicas, this run has "+std::to_string(plumed.getNumberOfReplicas()));
  }
  // all the replicas must restart from the same step
  if(plumed.getNumberOfReplicas()>1) {
    // no MPI type for signed long long in Communicator
    long long unsigned step0=static_cast<long long unsigned>(restartStep);
    if(comm.Get_rank()==0) {
      plumed.multi_sim_comm.Bcast(step0,0);
    }
    comm.Bcast(step0,0);
    if(static_cast<long long int>(step0)!=restartStep) {
      plumed_merror("checkpoint "+fname+" was written at step "+std::to_string(restartStep)+" but the checkpoint of replica 0 at step "+std::to_string(step0)+": all the replicas must restart from the same step");
    }
  }
  log<<"Restarting from checkpoint "<<fname<<" written at step "<<restartStep<<" by PLUMED "<<pending->getPlumedVersion()<<"\n";
  trimFiles(*pending);
  plumed.setRestart(true);
}

void Checkpointing::trimFiles(CheckpointReader& r) {
  auto& log=plumed.getLog();
  auto& comm=plumed.comm;
  r.beginSection("files");
  const unsigned n=r.get<unsigned>("n");
  if(comm.Get_rank()==0) {
    for(unsigned i=0; i<n; ++i) {
      const std::string num=std::to_string(i);
      const std::string f=r.get<std::string>("path"+num);
      const unsigned long long size=r.get<unsigned long long>("size"+num);
      const std::string crc=r.get<std::string>("crc"+num);
      if(r.get<bool>("snapshot"+num)) {
        // rewritten from scratch at every write: the restarted run replaces it
        log<<"  file "<<f<<" is rewritten at every write and is not truncated\n";
        continue;
      }
      std::error_code ec;
      if(!std::filesystem::exists(f,ec)) {
        if(size==0) {
          continue;
        }
        plumed_merror("checkpoint: file "+f+" is missing, but it had "+std::to_string(size)+" bytes when the checkpoint was written");
      }
      const unsigned long long current=std::filesystem::file_size(f,ec);
      if(ec) {
        plumed_merror("checkpoint: cannot get the size of "+f+": "+ec.message());
      }
      if(current<size) {
        plumed_merror("checkpoint: file "+f+" has "+std::to_string(current)+" bytes, fewer than the "+std::to_string(size)+" it had when the checkpoint was written");
      }
      if(current>size) {
        std::filesystem::resize_file(f,size,ec);
        if(ec) {
          plumed_merror("checkpoint: cannot truncate "+f+" to "+std::to_string(size)+" bytes: "+ec.message());
        }
        log<<"  file "<<f<<" truncated from "<<current<<" to "<<size<<" bytes\n";
      }
      if(size>0 && crcOfTail(f,size)!=crc) {
        plumed_merror("checkpoint: the content of "+f+" differs from the one recorded in the checkpoint");
      }
    }
  }
  // on every rank, since OFile keeps the fields on every rank
  for(unsigned i=0; i<n; ++i) {
    const std::string num=std::to_string(i);
    if(r.has("fields"+num)) {
      lastFields[r.get<std::string>("path"+num)]=r.get<std::string>("fields"+num);
    }
  }
  r.endSection();
  comm.Barrier();
}

void Checkpointing::restore() {
  auto& log=plumed.getLog();
  plumed_assert(pending);
  if(plumed.getStep()!=restartStep) {
    plumed_merror("the checkpoint was written at step "+std::to_string(restartStep)+" but the MD engine restarted at step "+std::to_string(plumed.getStep())+": the PLUMED checkpoint and the restart file of the MD engine do not refer to the same step");
  }
  if(!plumed.getRestart()) {
    plumed_merror("restarting from a checkpoint, but restart has been switched off (RESTART NO after CHECKPOINT?)");
  }
  pending->beginSection("actions");
  std::vector<std::string> labels;
  for(const auto & p : plumed.getActionSet()) {
    labels.push_back(p->getLabel());
    if(!p->hasCheckpointState()) {
      continue;
    }
    // per-action RESTART=NO: the action starts from scratch at the restart step
    if(!p->getRestart()) {
      log<<"  action "<<p->getLabel()<<": RESTART=NO, state not restored from the checkpoint, files backed up\n";
      continue;
    }
    if(!pending->hasSection(p->getLabel())) {
      plumed_merror("action "+p->getLabel()+" has no section in the checkpoint: was the checkpoint written with a different input? Use RESTART=NO on the action to start it from scratch");
    }
    pending->beginSection(p->getLabel());
    p->readCheckpoint(*pending);
    pending->endSection();
  }
  pending->endSection();
  for(const auto & s : pending->getSections()) {
    // sections are "actions/label/"
    if(s.compare(0,8,"actions/")!=0 || s.size()<=8) {
      continue;
    }
    const std::string label=s.substr(8,s.size()-9);
    if(label.find('/')!=std::string::npos) {
      continue;
    }
    if(std::find(labels.begin(),labels.end(),label)==labels.end()) {
      log<<"WARNING: the checkpoint has a section for action "<<label<<", which is not in the input\n";
    }
  }
  pending.reset();
  log<<"Checkpoint state restored at step "<<restartStep<<"\n";
}

void Checkpointing::atStepStart(bool signal) {
  if(pending) {
    restore();
  }
  if(!enabled()) {
    return;
  }
  const long long int step=plumed.getStep();
  if(step==lastWrittenStep || step==restartStep) {
    return;
  }
  if(signal || (stride>0 && step%stride==0)) {
    write();
  }
}

}
