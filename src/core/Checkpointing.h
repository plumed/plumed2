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
class Action;

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
kept with a _prev suffix. If the MD engine never sends the signal, a checkpoint
of the state after the last step (i.e. for the next step) is written at the end
of the run (cmd "runFinalJobs", or the destructor of PlumedMain), with a _final suffix.
*/
class Checkpointing {
  PlumedMain& plumed;
/// Base name given by the user (empty: checkpointing is off)
  std::string file;
/// Full path (with the replica suffix) and path of the previous checkpoint
  std::string path;
  std::string prevPath;
/// Path of the checkpoint written at the end of the run (atRunEnd)
  std::string finalPath;
  long long int stride=0;
/// A checkpointing signal was received from the MD engine in this run
  bool signalReceived=false;
/// The end of the run has been handled (atRunEnd called)
  bool runEndHandled=false;
/// The last step presented has been completed (update done)
  bool stepComplete=false;
/// Last step presented and difference between the last two steps (to know the next one)
  long long int lastStep=-1;
  long long int stepIncrement=1;
/// Paths of all the files opened for writing in this run
  std::vector<std::string> outputFiles;
/// Files rewritten from scratch at every write (OFile::rewind): they are not
/// trimmed at restart, since the next write replaces them anyway
  std::set<std::string> snapshotFiles;
/// Field headers recorded in the checkpoint being restored, or found when
/// trimming a file at the restart step, by path
  std::map<std::string,std::string> lastFields;
/// Checkpoint being restored, consumed at the first step
  std::unique_ptr<CheckpointReader> pending;
/// Step stored in the checkpoint being restored (-1: no restart)
  long long int restartStep=-1;
  long long int lastWrittenStep=-1;
  bool warnedAboutUnportedState=false;
/// Step from which the MD engine restarts (cmd "setRestartStep"), -1 if not given
  long long int engineRestartStep=-1;
/// The first step presented after a restart has been compared with engineRestartStep
  bool engineRestartStepChecked=false;
/// False with RESTART NO_TRIM
  bool trimEnabled=true;
/// Files already trimmed at the restart step (or left as they are)
  std::set<std::string> trimmedFiles;
  void write(long long int step,const std::string& outPath,const std::string& outPrevPath);
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
/// After a restart, restore the field header of the last line of a file appended
/// to (from the checkpoint, or found when trimming it), so that it is not printed again.
  void restoreLastFields(OFile& file) const;
/// Mark a file as rewritten from scratch at every write (called by OFile::rewind).
  void markSnapshotFile(const std::string& fileName);
/// To be called before any computation of the current step, possibly more
/// than once per step. Restores the pending checkpoint at the first step
/// after a restart, and writes a checkpoint when requested by the signal or
/// by the stride.
  void atStepStart(bool signal);
/// To be called at the end of the update of a step.
  void atStepEnd() {
    stepComplete=true;
  }
/// At the end of the run (cmd "runFinalJobs", before the final jobs, or, if it was not
/// called, when PLUMED is finalized, fromDestructor): if the MD engine never sent the
/// checkpointing signal, write the state after the last step, i.e. at the start of the
/// next one, to FILE with the suffix _final (plumed_final.cpt). Only the first call does
/// something; nothing is written if the last step was not completed or, from the
/// destructor, if MPI is already finalized.
  void atRunEnd(bool fromDestructor=false);
/// True after restartFrom() has been called (a checkpoint restart is in progress or done).
  bool isRestartingFromCheckpoint() const {
    return restartStep>=0;
  }
  const std::string& getPath() const {
    return path;
  }
/// plumed.cpt -> plumed_prev.cpt for the tag _prev (without the replica suffix)
  static std::string taggedName(const std::string& fileName,const std::string& tag);
/// Name of the checkpoint written at the end of the run (plumed_final.cpt), without the replica suffix
  static std::string finalName(const std::string& fileName);
/// Full path of the checkpoint file for this replica, given its base name.
  std::string suffixedPath(const std::string& fileName) const;
/// Step from which the MD engine restarts (cmd "setRestartStep", before init).
  void setEngineRestartStep(long long int step);
  long long int getEngineRestartStep() const {
    return engineRestartStep;
  }
/// RESTART NO_TRIM: the files are appended as they are.
  void disableTrimming() {
    trimEnabled=false;
  }
/// True when the files appended at restart are trimmed at the restart step:
/// the MD engine gave the step, trimming was not switched off, and the run does
/// not restart from a checkpoint (which trims the files to the recorded sizes).
  bool trimsAtRestart() const {
    return engineRestartStep>=0 && trimEnabled && !isRestartingFromCheckpoint();
  }
/// Remove from a file the lines written at or after the restart step, i.e. the
/// last lines whose time is not smaller than the time of the restart step (the
/// earlier lines, also those of earlier runs whose time restarted from zero, are
/// kept). Files without a time field, compressed files, and files where the
/// time of these lines is not increasing are left as they are, with a line in
/// the log. Each file is processed only once. To be called by one process only,
/// before the file is read or written by this run (see OFile).
/// keepRestartStep: the lines of the restart step itself are kept (e.g. the hill of
/// that step, when the bias is read from a grid written at the checkpoint).
  void trimAtRestart(const std::string& fileName,bool keepRestartStep=false);
/// For an action that reads one of its output files at restart before appending
/// to it (e.g. the hills of METAD): trim the file at the restart step first, if
/// it is trimmed (see trimAtRestart and NO_TRIM). The path is the one the action
/// writes: name with the replica suffix, unless nosuffix. To be called by all the
/// processes of the action.
  void trimOutputFileAtRestart(Action& action,const std::string& name,bool nosuffix=false);
};

}

#endif
