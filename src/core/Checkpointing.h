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
#ifndef __PLUMED_core_Checkpointing_h
#define __PLUMED_core_Checkpointing_h

#include <string>
#include <vector>
#include <memory>
#include <map>
#include <set>

namespace PLMD {

class PlumedMain;
class CheckpointReader;
class OFile;

/**
Exact-restart checkpoint of a PLUMED run, configured by the CHECKPOINT action.

The checkpoint for step S captures the state at the start of S (after the
update of step S-1 and before anything of S is computed) together with the
sizes of all the output files after a flush at that moment. On restart the
files are truncated to those sizes and the MD engine presents step S again,
which is replayed in full. No action needs to know that a step is repeated.

The checkpoint is written when the MD engine sends the checkpointing signal
(cmd "doCheckPoint") or, with STRIDE, every STRIDE steps. It is a separate
file, one per replica (usual suffix), written atomically; the previous one is
kept with a _prev suffix.
*/
class Checkpointing {
  PlumedMain& plumed;
/// Base name given by the user (empty: checkpointing is off)
  std::string file;
/// Full path (with the replica suffix) and path of the previous checkpoint
  std::string path;
  std::string prevPath;
  long long int stride=0;
/// Paths of all the files opened for writing in this run
  std::vector<std::string> outputFiles;
/// Files rewritten from scratch at every write (OFile::rewind): they are not
/// trimmed at restart, since the next write replaces them anyway
  std::set<std::string> snapshotFiles;
/// Field headers recorded in the checkpoint being restored, by path
  std::map<std::string,std::string> lastFields;
/// Checkpoint being restored, consumed at the first step
  std::unique_ptr<CheckpointReader> pending;
/// Step stored in the checkpoint being restored (-1: no restart)
  long long int restartStep=-1;
  long long int lastWrittenStep=-1;
  bool warnedAboutUnportedState=false;
  void write();
  void restore();
  void trimFiles(CheckpointReader& reader);
  static std::string crcOfTail(const std::string& fname,unsigned long long size);
public:
  explicit Checkpointing(PlumedMain& plumed);
  ~Checkpointing();
  bool enabled() const {
    return !file.empty();
  }
/// Configure checkpoint writing (called by the CHECKPOINT action).
  void setup(const std::string& fileName,long long int checkpointStride);
/// Read a checkpoint, trim the output files and switch restart on
/// (called by the CHECKPOINT action, before any other action opens a file).
  void restartFrom(const std::string& fileName);
/// Record a file opened for writing (called by OFile::open) and, after a
/// restart, restore the header of its last line so that it is not printed again.
  void registerOutputFile(OFile& file);
/// Mark a file as rewritten from scratch at every write (called by OFile::rewind).
  void markSnapshotFile(const std::string& fileName);
/// To be called before any computation of the current step, possibly more
/// than once per step. Restores the pending checkpoint at the first step
/// after a restart, and writes a checkpoint when requested by the signal or
/// by the stride.
  void atStepStart(bool signal);
/// True after restartFrom() has been called (a checkpoint restart is in progress or done).
  bool isRestartingFromCheckpoint() const {
    return restartStep>=0;
  }
  const std::string& getPath() const {
    return path;
  }
/// Full path of the checkpoint file for this replica, given its base name.
  std::string suffixedPath(const std::string& fileName) const;
};

}

#endif
