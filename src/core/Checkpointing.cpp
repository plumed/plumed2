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
#include "ActionWithValue.h"
#include "tools/Checkpoint.h"
#include "tools/FileTools.h"
#include "tools/Communicator.h"
#include "tools/Exception.h"
#include "tools/OFile.h"
#include "tools/FileBase.h"
#include "tools/Tools.h"
#include "tools/H5Table.h"
#include "config/Config.h"
#include <algorithm>
#include <filesystem>
#include <fstream>
#include <limits>

#ifdef __PLUMED_HAS_ZLIB
#include <zlib.h>
#endif

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

std::string Checkpointing::taggedName(const std::string& fileName,const std::string& tag) {
  // plumed.cpt -> plumed_prev.cpt, plumed_final.cpt (the replica suffix, if any, goes after them)
  std::string name=fileName;
  const size_t slash=name.find_last_of("/\\");
  const size_t dot=name.find_last_of('.');
  if(dot!=std::string::npos && (slash==std::string::npos || dot>slash)) {
    name.insert(dot,tag);
  } else {
    name+=tag;
  }
  return name;
}

std::string Checkpointing::finalName(const std::string& fileName) {
  return taggedName(fileName,"_final");
}

void Checkpointing::setup(const std::string& fileName,long long int checkpointStride) {
  plumed_massert(!enabled(),"only one CHECKPOINT action is allowed");
  plumed_massert(!fileName.empty(),"CHECKPOINT: FILE cannot be empty");
  plumed_massert(outputFiles.empty(),"CHECKPOINT must come before any action that writes a file");
  plumed_massert(Tools::extension(fileName)!="gz","CHECKPOINT: the checkpoint file cannot be compressed");
  file=fileName;
  stride=checkpointStride;
  path=suffixedPath(file);
  prevPath=suffixedPath(taggedName(file,"_prev"));
  finalPath=suffixedPath(finalName(file));
  // a new run does not overwrite the checkpoints of a previous one
  if(!pending) {
    OFile f;
    f.link(plumed);
    f.backupFile("bck",path);
    f.backupFile("bck",prevPath);
    f.backupFile("bck",finalPath);
  }
}

void Checkpointing::restoreLastFields(OFile& f) const {
  if(f.getMode()=="a") {
    const auto it=lastFields.find(f.getPath());
    if(it!=lastFields.end()) {
      f.setLastFields(it->second);
    }
  }
}

void Checkpointing::registerOutputFile(OFile& f) {
  const std::string fileName=f.getPath();
  restoreLastFields(f);
  if(std::find(outputFiles.begin(),outputFiles.end(),fileName)!=outputFiles.end()) {
    return;
  }
  // devices and pipes (e.g. /dev/null, written by the processes or replicas that do not write a shared file)
  // have no size to record
  std::error_code ec;
  if(std::filesystem::exists(pathFromUTF8(fileName),ec) && !std::filesystem::is_regular_file(pathFromUTF8(fileName),ec)) {
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

void Checkpointing::setEngineRestartStep(long long int step) {
  if(step<0) {
    plumed_merror("the restart step should not be negative");
  }
  engineRestartStep=step;
}

namespace {

/// The lines of a plain or gzipped file (gzipped files written in append mode are made
/// of several gzip members, which are read as a single stream)
class TrimLineReader {
  std::ifstream in;
  std::unique_ptr<h5table::Reader> h5reader;
#ifdef __PLUMED_HAS_ZLIB
  gzFile gz=nullptr;
#endif
public:
  TrimLineReader(const std::string& fname,bool compressed) {
    if(h5table::isH5(fname)) {
      h5reader=std::make_unique<h5table::Reader>(fname);
    } else if(compressed) {
#ifdef __PLUMED_HAS_ZLIB
      gz=(gzFile)gzopenUTF8(fname,"rb");
      if(!gz) {
        plumed_merror("cannot open "+fname+" to trim it at the restart step");
      }
#endif
    } else {
      in.open(pathFromUTF8(fname),std::ios::binary);
      if(!in) {
        plumed_merror("cannot open "+fname+" to trim it at the restart step");
      }
    }
  }
  ~TrimLineReader() {
#ifdef __PLUMED_HAS_ZLIB
    if(gz) {
      gzclose(gz);
    }
#endif
  }
  TrimLineReader(const TrimLineReader&)=delete;
  TrimLineReader& operator=(const TrimLineReader&)=delete;
  /// the next line without its newline; newline tells if it had one (the last line may not)
  bool getline(std::string& line,bool& newline) {
    if(h5reader) {
      newline=true;
      return h5reader->getline(line);
    }
#ifdef __PLUMED_HAS_ZLIB
    if(gz) {
      line.clear();
      char buf[4096];
      while(gzgets(gz,buf,sizeof(buf))) {
        line+=buf;
        if(!line.empty() && line.back()=='\n') {
          line.pop_back();
          newline=true;
          return true;
        }
      }
      newline=false;
      return !line.empty();
    }
#endif
    if(!std::getline(in,line)) {
      return false;
    }
    newline=!in.eof();
    return true;
  }
};

/// A free name for the lines removed from a file: bck.trimmed.N.name, in the directory of the file
std::string trimmedBackupName(const std::string& fname) {
  const std::size_t found=fname.find_last_of("/\\");
  const std::string directory=fname.substr(0,found+1);
  const std::string file=fname.substr(found+1);
  for(int i=0;; i++) {
    const std::string b=directory+"bck.trimmed."+std::to_string(i)+"."+file;
    std::error_code ec;
    if(!std::filesystem::exists(pathFromUTF8(b),ec)) {
      return b;
    }
  }
}

#ifdef __PLUMED_HAS_ZLIB
/// copy n uncompressed bytes (all the remaining ones if n<0) from in to out
void gzCopy(gzFile in,gzFile out,long long n,const std::string& fname) {
  char buf[65536];
  while(n!=0) {
    const unsigned want=(n<0 || n>static_cast<long long>(sizeof(buf))) ? sizeof(buf) : static_cast<unsigned>(n);
    const int r=gzread(in,buf,want);
    if(r<0) {
      plumed_merror("cannot read "+fname+" to trim it at the restart step");
    }
    if(r==0) {
      break;
    }
    if(gzwrite(out,buf,r)!=r) {
      plumed_merror("cannot write while trimming "+fname+" at the restart step");
    }
    if(n>0) {
      n-=r;
    }
  }
}
#endif

}

void Checkpointing::trimAtRestart(const std::string& fname,bool keepRestartStep) {
  if(!trimsAtRestart() || fname=="/dev/null" || !trimmedFiles.insert(fname).second) {
    return;
  }
  auto& log=plumed.getLog();
  std::error_code ec;
  if(!std::filesystem::exists(pathFromUTF8(fname),ec)) {
    return;
  }
  const bool compressed=Tools::extension(fname)=="gz";
#ifndef __PLUMED_HAS_ZLIB
  if(compressed) {
    log<<"  WARNING: file "<<fname<<" is compressed and PLUMED has no zlib: it is not trimmed at the restart step\n";
    return;
  }
#endif
  const auto ts=plumed.getActionSet().selectWithLabel<ActionWithValue*>("timestep");
  const double dt=ts ? ts->copyOutput(0)->get() : 0.0;
  if(!(dt>0.0)) {
    log<<"  WARNING: the timestep is not known, file "<<fname<<" is not trimmed at the restart step\n";
    return;
  }
  // the lines of the restart step are removed too, since the step is computed again
  // (unless keepRestartStep); times are printed with a few digits, hence the half timestep
  const double tcut=(engineRestartStep+(keepRestartStep?0.5:-0.5))*dt;
  // keepEnd: end of the last line with a time before tcut (the tail starts there,
  // including the comment lines that follow it); offsets in uncompressed bytes
  long long offset=0,keepEnd=0;
  long long ntail=0;
  int timeIndex=-1;
  bool hasTime=false,inTail=false,increasing=true;
  double lastTime=0.0;
  // the field header (FIELDS and SET lines) of the current section, and the one of the last line kept,
  // in the format of OFile::getLastFields: the restarted run does not print it again if it is the same
  std::string header,keptHeader;
  // for HDF5 tables: the block (a FIELDS line starts one) and the row of the last line kept
  long long block=-1,row=0,keepBlock=-1,keepRow=0;
  const bool table=h5table::isH5(fname);
  {
    TrimLineReader reader(fname,compressed);
    std::string line;
    bool newline=false;
    std::vector<std::string> words;
    while(reader.getline(line,newline)) {
      const long long end=offset+line.size()+(newline?1:0);
      if(line.compare(0,9,"#! FIELDS")==0) {
        block++;
        row=0;
        words=Tools::getWords(line.substr(9));
        timeIndex=-1;
        header.clear();
        for(unsigned i=0; i<words.size(); ++i) {
          if(words[i]=="time" && timeIndex<0) {
            timeIndex=i;
            hasTime=true;
          }
          header+=words[i]+"\t0\t\n";
        }
      } else if(line.compare(0,7,"#! SET ")==0) {
        // "#! SET name value": the value is kept as it is printed (it may have leading spaces, from FMT)
        const std::size_t sp=line.find(' ',7);
        if(sp!=std::string::npos) {
          header+=line.substr(7,sp-7)+"\t1\t"+line.substr(sp+1)+"\n";
        }
      } else if(!line.empty() && line[0]!='#' && timeIndex>=0) {
        row++;
        words=Tools::getWords(line);
        double t;
        if(words.size()>static_cast<unsigned>(timeIndex) && Tools::convertNoexcept(words[timeIndex],t)) {
          if(t<tcut) {
            keepEnd=end;
            keptHeader=header;
            keepBlock=block;
            keepRow=row;
            inTail=false;
            ntail=0;
            increasing=true;
          } else {
            if(inTail && t<lastTime) {
              increasing=false;
            }
            inTail=true;
            ntail++;
            lastTime=t;
          }
        } else if(inTail) {
          // an incomplete line (an interrupted run) at the end of the tail
          ntail++;
        }
      }
      offset=end;
    }
  }
  if(!hasTime) {
    log<<"  WARNING: file "<<fname<<" has no time field and is not trimmed at the restart step\n";
    return;
  }
  if(inTail && !increasing) {
    log<<"  WARNING: in file "<<fname<<" the time of the lines after the restart step is not increasing: the file is not trimmed\n";
    return;
  }
  // the file ends with the header of its last line: it is not printed again if the fields are the same
  // (an HDF5 table continues its last block if the fields are the same)
  if(keptHeader.length()>0 && !table) {
    lastFields[fname]=keptHeader;
  }
  if(!inTail) {
    return;
  }
  // the removed lines are kept in a separate file
  const std::string backup=trimmedBackupName(fname);
  if(table) {
    // the blocks before the one of the last line kept, and its rows up to that line
    const auto rows=h5table::rows(fname);
    std::vector<unsigned long long> keep;
    for(long long b=0; b<keepBlock; b++) {
      keep.push_back(rows[b]);
    }
    if(keepBlock>=0) {
      keep.push_back(keepRow);
    }
    h5table::truncate(fname,keep,backup);
  } else if(!compressed) {
    {
      std::ifstream in(pathFromUTF8(fname),std::ios::binary);
      std::ofstream out(pathFromUTF8(backup),std::ios::binary);
      in.seekg(keepEnd);
      out<<in.rdbuf();
      if(!out) {
        plumed_merror("cannot write "+backup+" while trimming "+fname+" at the restart step");
      }
    }
    std::filesystem::resize_file(pathFromUTF8(fname),keepEnd,ec);
    if(ec) {
      plumed_merror("cannot trim "+fname+" at the restart step: "+ec.message());
    }
  } else {
#ifdef __PLUMED_HAS_ZLIB
    // a gzipped file cannot be truncated: the kept part is compressed again in a temporary file,
    // which is then copied on the file (in place, so that the file is the same for whoever has it open)
    const std::string tmp=fname+".trim.tmp";
    {
      gzFile in=(gzFile)gzopenUTF8(fname,"rb");
      gzFile kept=(gzFile)gzopenUTF8(tmp,"wb9");
      gzFile removed=(gzFile)gzopenUTF8(backup,"wb9");
      if(!in || !kept || !removed) {
        plumed_merror("cannot open the files to trim "+fname+" at the restart step");
      }
      gzCopy(in,kept,keepEnd,fname);
      gzCopy(in,removed,-1,fname);
      gzclose(in);
      if(gzclose(kept)!=Z_OK || gzclose(removed)!=Z_OK) {
        plumed_merror("cannot write while trimming "+fname+" at the restart step");
      }
    }
    {
      std::ifstream in(pathFromUTF8(tmp),std::ios::binary);
      std::ofstream out(pathFromUTF8(fname),std::ios::binary|std::ios::trunc);
      out<<in.rdbuf();
      if(!out) {
        plumed_merror("cannot write "+fname+" while trimming it at the restart step");
      }
    }
    std::filesystem::remove(pathFromUTF8(tmp),ec);
#endif
  }
  log<<"  file "<<fname<<": "<<ntail<<" lines written "<<(keepRestartStep?"after":"at or after")<<" the restart step "<<engineRestartStep<<" moved to "<<backup<<"\n";
}

void Checkpointing::trimOutputFileAtRestart(Action& action,const std::string& name,bool nosuffix) {
  if(action.getNoTrim() || !trimsAtRestart()) {
    return;
  }
  const std::string fname=nosuffix ? name : FileBase::appendSuffix(name,plumed.getSuffix());
  if(action.comm.Get_rank()==0) {
    trimAtRestart(fname);
  }
  action.comm.Barrier();
}

std::string Checkpointing::crcOfTail(const std::string& fname,unsigned long long size) {
  std::ifstream in(pathFromUTF8(fname),std::ios::binary);
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

void Checkpointing::write(long long int step,const std::string& outPath,const std::string& outPrevPath) {
  plumed_assert(enabled());
  auto& log=plumed.getLog();
  auto& comm=plumed.comm;
  // the files still to be trimmed at the restart step are trimmed before their sizes are measured
  for(const auto & fb : plumed.getFiles()) {
    if(auto of=dynamic_cast<OFile*>(fb)) {
      of->applyPendingTrim();
    }
  }
  // everything written so far by step-1 must be on disk before the sizes are measured
  plumed.fflush();
  log.flush();
  for(const auto & p : plumed.getActionSet()) {
    p->fflush();
  }
  CheckpointWriter w(config::getVersionLong());
  w.write("step",step);
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
      if(std::filesystem::exists(pathFromUTF8(f),ec)) {
        size=std::filesystem::file_size(pathFromUTF8(f),ec);
        if(ec) {
          plumed_merror("checkpoint: cannot get the size of "+f+": "+ec.message());
        }
      }
      w.write("path"+num,f);
      if(h5table::isH5(f)) {
        // an HDF5 table: the rows of its blocks, its size in bytes is meaningless
        std::vector<long long> rows;
        for(const auto & fb : plumed.getFiles()) {
          OFile* of=dynamic_cast<OFile*>(fb);
          if(of && of->isOpen() && of->getPath()==f) {
            for(const auto r : of->getH5Rows()) {
              rows.push_back(r);
            }
            break;
          }
        }
        w.write("h5rows"+num,rows);
        size=0;
      }
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
    w.writeFile(outPath,outPrevPath);
  }
  comm.Barrier();
  lastWrittenStep=step;
  log<<"Checkpoint written to "<<outPath<<" at step "<<step<<"\n";
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
    // The byte broadcast below is limited to int counts. Use that same type
    // for its size, also on hosts whose MPI subset lacks unsigned long long.
    int n=content.size()<static_cast<std::size_t>(std::numeric_limits<int>::max())
          ? static_cast<int>(content.size()) : -1;
    comm.Bcast(n,0);
    plumed_massert(n>=0,"checkpoint file too large");
    content.resize(n);
    comm.Bcast(content.data(),n,0);
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
      if(r.has("h5rows"+num)) {
        // an HDF5 table: the blocks are cut to the rows they had
        std::vector<long long> rows;
        r.read("h5rows"+num,rows);
        if(!std::filesystem::exists(pathFromUTF8(f),ec)) {
          if(rows.empty()) {
            continue;
          }
          plumed_merror("checkpoint: file "+f+" is missing, but it had rows when the checkpoint was written");
        }
        const auto current=h5table::rows(f);
        if(current.size()<rows.size()) {
          plumed_merror("checkpoint: file "+f+" has fewer blocks than when the checkpoint was written");
        }
        bool longer=current.size()>rows.size();
        std::vector<unsigned long long> keep;
        for(unsigned b=0; b<rows.size(); b++) {
          if(current[b]<static_cast<unsigned long long>(rows[b])) {
            plumed_merror("checkpoint: file "+f+" has fewer rows than when the checkpoint was written");
          }
          longer=longer || current[b]>static_cast<unsigned long long>(rows[b]);
          keep.push_back(rows[b]);
        }
        if(longer) {
          h5table::truncate(f,keep,"");
          log<<"  file "<<f<<" truncated to the rows it had at the checkpoint\n";
        }
        continue;
      }
      if(!std::filesystem::exists(pathFromUTF8(f),ec)) {
        if(size==0) {
          continue;
        }
        plumed_merror("checkpoint: file "+f+" is missing, but it had "+std::to_string(size)+" bytes when the checkpoint was written");
      }
      const unsigned long long current=std::filesystem::file_size(pathFromUTF8(f),ec);
      if(ec) {
        plumed_merror("checkpoint: cannot get the size of "+f+": "+ec.message());
      }
      if(current<size) {
        plumed_merror("checkpoint: file "+f+" has "+std::to_string(current)+" bytes, fewer than the "+std::to_string(size)+" it had when the checkpoint was written");
      }
      if(current>size) {
        std::filesystem::resize_file(pathFromUTF8(f),size,ec);
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
  // the final checkpoint of the previous run is out of date once the run continues
  // (RESTART_FROM=AUTO reads it when it exists)
  if(enabled()) {
    OFile f;
    f.link(plumed);
    f.backupFile("bck",finalPath);
  }
}

void Checkpointing::atStepStart(bool signal) {
  if(pending) {
    restore();
  }
  const long long int step=plumed.getStep();
  if(!engineRestartStepChecked && trimsAtRestart()) {
    engineRestartStepChecked=true;
    if(step!=engineRestartStep) {
      plumed_merror("the MD engine said that it restarts from step "+std::to_string(engineRestartStep)+", and the files were trimmed at that step, but the first step is "+std::to_string(step));
    }
  }
  if(!enabled()) {
    return;
  }
  signalReceived=signalReceived || signal;
  stepComplete=false;
  if(step!=lastStep) {
    if(lastStep>=0 && step>lastStep) {
      stepIncrement=step-lastStep;
    }
    lastStep=step;
  }
  if(step==lastWrittenStep || step==restartStep) {
    return;
  }
  if(signal || (stride>0 && step%stride==0)) {
    write(step,path,prevPath);
  }
}

void Checkpointing::atRunEnd(bool fromDestructor) {
  if(runEndHandled) {
    return;
  }
  runEndHandled=true;
  // nothing can be written (nor logged) after MPI_Finalize; finalizing PLUMED after MPI is
  // not allowed anyway (the communicators are freed)
  if(fromDestructor && Communicator::initialized() && Communicator::finalized()) {
    return;
  }
  // an MD engine that sends the checkpointing signal writes its own restart files
  // at the steps of the signals (GROMACS also at the last step): a checkpoint
  // after the last step would not match them
  if(!enabled() || signalReceived || lastStep<0) {
    return;
  }
  auto& log=plumed.getLog();
  if(!stepComplete) {
    log<<"Final checkpoint not written: the last step ("<<lastStep<<") was not completed\n";
    return;
  }
  // the state after the update of the last step is the one at the start of the next step
  write(lastStep+stepIncrement,finalPath,"");
  plumed.getLog()<<"  (final checkpoint, read at restart with RESTART_FROM=AUTO: the MD engine must restart from step "<<lastStep+stepIncrement<<")\n";
}

}
