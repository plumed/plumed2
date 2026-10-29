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
#ifndef __PLUMED_bias_MetaDHills_h
#define __PLUMED_bias_MetaDHills_h

#include "MetaDKernel.h"
#include "tools/Communicator.h"
#include <deque>
#include <functional>
#include <limits>
#include <memory>
#include <string>
#include <vector>

namespace PLMD {

class ActionWithValue;
class IFile;
class OFile;
class Value;

namespace bias {

class MetaDBiasGrid;
class MetaDDiffusionRestart;

/// shortest representation of x that is converted back to the same double
std::string exactString(double x);

/// broadcast a string from process 0: the length is sent first, as on the other processes the string
/// can have a different length (Communicator::Bcast of a string uses the local length)
void bcastString(Communicator& c, std::string& s);

/// remove the last line of a hills file if it is incomplete (left by a run that was interrupted while
/// writing it); returns the message for the log (empty if nothing was done)
std::string removeIncompleteLastLine(const std::string& path);

/// The hills files of METAD and PBMETAD for a set of arguments: writing of the hills, and reading of the files
/// of a restart or of the other walkers.
///
/// With several MPI processes the files are opened and read by the first process only, which also decides if
/// they exist, and the hills are broadcast to the others: the file of another walker can grow while it is read,
/// and if each process read it they could read a different number of hills, and their biases would differ.
/// The file operations do not synchronize the processes (the files are linked to a communicator of the process
/// alone). Everything written on the log is written by all the processes.
class MetaDHillsFiles {
public:
  /// called for each hill that is read: time, multivariate, height, center, sigma (the widths, or the upper
  /// triangle of the inverse covariance matrix for multivariate hills)
  using Callback=std::function<void(double,bool,double,const std::vector<double>&,const std::vector<double>&)>;
  /// The bias to which the hills that are read are added (see readHills)
  struct Bias {
    /// add a hill to the bias
    std::function<void(const MetaDGaussian&)> add;
    /// the grid of the bias, if any (it completes the line of the log when it is extended)
    MetaDBiasGrid* grid=nullptr;
    /// ADAPTIVE=DIFF: when it is pending, the hills of the own file are kept here for the statistics of the fluctuations
    MetaDDiffusionRestart* diffrestart=nullptr;
    /// the height of a hill is the one in the file times heightscale (well-tempered metadynamics)
    double heightscale=1.0;
    /// times closer than this are the same time (times are written with a finite precision): half a PACE
    double tolerance=0.0;
  };
private:
  ActionWithValue& action_;
  Communicator& comm_;
  /// communicator of the process alone, linked to the input files
  Communicator selfcomm_;
  std::vector<Value*> args_;
  std::vector<std::unique_ptr<IFile>> ifiles_;
  std::vector<std::string> ifilesnames_;
  /// number of hills read from each file
  std::vector<unsigned long long> nread_;
  Bias bias_;
  bool noStretchWarningDone_=false;
  /// a hill with a Gaussian (not stretched) kernel was read
  bool gaussianKernelsFound_=false;
  void noStretchWarning();
  /// an error in a hills file
  [[noreturn]] void inputError(const std::string& msg) const;
  /// the numbers of hills of the other walkers in a grid file (on one process)
  static std::vector<long long> readGridComment(const std::string& path, unsigned nwalkers);
  bool scanOneHill(IFile* ifile, std::vector<Value>& v, std::vector<double>& center, std::vector<double>& sigma, double& height,
                   bool& multivariate, double& time);
public:
  MetaDHillsFiles(ActionWithValue& action, Communicator& comm, const std::vector<Value*>& args);
  /// the name of the hills file of a walker: name.walker with WALKERS_N, in the directory dir if given
  /// (with WALKERS_N or WALKERS_MPI)
  static std::string getWalkerFileName(const std::string& name, const std::string& dir, int nwalkers, int walker, bool mpi);
  ~MetaDHillsFiles();
  /// add an input file, which is opened if it exists; returns its index. With nosuffix, the name is used without
  /// the suffix of the replica (with multiple walkers, as for the files that are written)
  unsigned addInput(const std::string& name, bool nosuffix=false);
  const std::string& getInputName(unsigned i) const {
    return ifilesnames_[i];
  }
  /// true if the file is open (the same on all the processes)
  bool isOpen(unsigned i);
  /// path of an open file (on the first process only, empty on the others)
  std::string getPath(unsigned i) const;
  /// open the file if it exists now; returns true if it has been opened
  bool openIfExists(unsigned i);
  /// read the hills that are in the file (an open file) up to now and call f for each one; returns their number.
  /// An error while reading (e.g. a wrong field) is thrown by all the processes.
  unsigned read(unsigned i, const Callback& f);
  /// With WALKERS_N, the grid written by a walker contains the hills it has read from the files of the other walkers:
  /// their numbers are written in a comment line after the grid, so that when restarting from the grid only the
  /// hills that are not in it are added. The comment line (own is the index of this walker)
  std::string getGridComment(unsigned own) const;
  /// the numbers of hills of the other walkers in the grid file path (written by getGridComment), -1 for this walker;
  /// empty if the file does not contain them
  /// (read by the first process of comm and broadcast)
  static std::vector<long long> readGridComment(Communicator& comm, const std::string& path, unsigned nwalkers);
  void setBias(const Bias& bias) {
    bias_=bias;
  }
  const Bias& getBias() const {
    return bias_;
  }
  /// Read the hills of file i that are in the file up to now (an open file) and add them to the bias, writing their
  /// number on the log (the caller has started the line). Returns their number.
  /// lasttime: if given, the time of the last hill is returned, and the hills are kept for ADAPTIVE=DIFF (own file);
  /// addtobias: false to only scan the file; the hills up to time skipupto, and the first skipcount, are not added
  /// (they are already in the grid read from file)
  unsigned readHills(unsigned i, double* lasttime=nullptr, bool addtobias=true,
                     double skipupto=-std::numeric_limits<double>::infinity(), unsigned long long skipcount=0);
  /// prepare an open file for the next read (the hills read up to now are not read again)
  void reset(unsigned i);
  void close(unsigned i);
  /// write a hill; height is the one written in the file
  void write(OFile& file, double time, const MetaDGaussian& hill, double height, double biasf, bool clock) const;
};

/// ADAPTIVE=DIFF: when restarting, the statistics of the fluctuations of the variables (running average and
/// covariance) are initialized from the last hills of this walker in its hills file. The hills are kept here while the
/// file is read, until the first update, when the current time is known.
class MetaDDiffusionRestart {
  std::deque<double> time_;
  std::deque<std::vector<double>> center_;
  /// as stored in the hills: upper triangle of the inverse covariance matrix
  std::deque<std::vector<double>> sigma_;
  std::string file_;
public:
  /// the statistics are still to be initialized
  bool pending=false;
  /// maximum number of hills to keep (the older ones have a negligible weight)
  unsigned maxhills=0;
  void clear() {
    time_.clear();
    center_.clear();
    sigma_.clear();
  }
  /// a hill read from file (hills that are not multivariate cannot be used: the previous ones are discarded)
  void add(const std::string& file, double time, bool multivariate, const std::vector<double>& center, const std::vector<double>& sigma);
  const std::string& getFile() const {
    return file_;
  }
  /// The statistics from the hills of this walker. The hills with times within tol form blocks (one deposition):
  /// with blocksize 0 all the hills are used, otherwise the hill at position offset of the blocks with blocksize hills.
  /// Only the hills up to now, and the last sequence with increasing times, are used. The covariance is that of the
  /// last hill; the average is the average of the centers with the exponential weights of the running average
  /// (decay per update, updatetime is the time between updates). difference(i,a,b) and bringBackInPbc(i,x) are those
  /// of variable i. Returns the number of hills used (0: no statistics), and the time of the last one.
  unsigned compute(double tol, double now, unsigned blocksize, unsigned offset, double decay, double updatetime,
                   const std::function<double(unsigned,double,double)>& difference,
                   const std::function<double(unsigned,double)>& bringBackInPbc,
                   std::vector<double>& average, std::vector<double>& variance, double& lasttime) const;
};

/// The steps when a hill is added: every stride steps, except the first step of a run. At step 0 no hill is
/// added (also when the bias of a previous simulation is read with RESTART). When continuing a simulation, the
/// hill is not added only if the hills file already contains it, i.e. if the previous run already performed this
/// step (some MD codes compute again the last step of the previous run when restarting, others do not).
class MetaDHillSchedule {
  bool runStarted_=false;
  long long int runStartStep_=0;
  /// the time of the last hill in the hills file read when restarting
  bool lastHillTimeKnown_=false;
  double lastHillTime_=0.0;
public:
  void setLastHillTime(double time) {
    lastHillTimeKnown_=true;
    lastHillTime_=time;
  }
  /// to be called at every calculate (the first call starts the run)
  void start(long long int step) {
    if(!runStarted_) {
      runStarted_=true;
      runStartStep_=step;
    }
  }
  /// true before the first calculate and at the first step of the run
  bool isRunFirstStep(long long int step) const {
    return !runStarted_ || step==runStartStep_;
  }
  /// tolerance: times are written with a finite precision (FMT), so that a hill within tolerance (half a PACE)
  /// is considered to be the hill of this step; a later hill means that the previous run went beyond this step
  bool isHillStep(long long int step, long long int stride, double time, double tolerance) const {
    if(step%stride!=0) {
      return false;
    }
    if(isRunFirstStep(step)) {
      if(step==0 || !lastHillTimeKnown_) {
        return false;
      }
      return lastHillTime_<time-tolerance;
    }
    return true;
  }
};

}
}

#endif
