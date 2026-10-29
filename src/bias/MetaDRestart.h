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
#ifndef __PLUMED_bias_MetaDRestart_h
#define __PLUMED_bias_MetaDRestart_h

#include <string>

namespace PLMD {

class Action;
class Keywords;
class OFile;

namespace bias {

/// The multiple walkers of METAD and PBMETAD: they share their hills through the hills files (WALKERS_N), each
/// walker reading the files of the others every WALKERS_RSTRIDE steps, or through MPI (WALKERS_MPI), the walkers
/// being the replicas, which all add the hills of all of them and write them in the same files.
struct MetaDWalkers {
  /// number of walkers sharing hills through files (WALKERS_N)
  int n=1;
  /// directory with the hills files (WALKERS_DIR; with WALKERS_MPI the directory of the first replica, if not given)
  std::string dir;
  /// id of this walker (WALKERS_ID)
  int id=0;
  /// stride for reading the hills files of the other walkers (WALKERS_RSTRIDE)
  int rstride=1;
  /// MPI version (WALKERS_MPI)
  bool mpi=false;
  /// number of MPI walkers, and index of this one (the replica index)
  unsigned mpi_nw=0;
  unsigned mpi_id=0;
  static void registerKeywords(Keywords& keys);
  /// parse the WALKERS_* keywords
  void parse(Action& action);
  /// check the keywords, write them on the log, and set the MPI walkers
  void setup(Action& action);
  /// the walkers share files (hills and grids), which are written and read without the suffix of the replica
  bool shareFiles() const {
    return n>1 || mpi;
  }
  /// with WALKERS_MPI the shared files are written by the first walker only
  bool writesFiles() const {
    return !(mpi && mpi_id>0);
  }
  /// with WALKERS_N, the steps when the hills files of the other walkers are read (and the own one flushed)
  bool isReadStep(long long int step) const {
    return n>1 && step%rstride==0;
  }
  /// the name of the file of a walker: name.walker with WALKERS_N, in dir if given (WALKERS_N or WALKERS_MPI)
  std::string getFileName(const std::string& name, int walker) const;
  /// the current directory of the first replica, so that all the walkers use the same files (the replicas
  /// can run in different directories); to be called by all the processes of all the replicas
  std::string firstReplicaDirectory(Action& action) const;
  /// with WALKERS_MPI and RESTART, an error if the walkers did not all make the same choice (e.g. reading a grid);
  /// to be called by all the processes of all the replicas
  void checkSameChoice(Action& action, bool choice, const std::string& what) const;
  /// open a file written by this walker (hills or grid): shared files are opened without the suffix of the replica,
  /// and the walkers that do not write them write on /dev/null
  void openOutput(Action& action, OFile& file, const std::string& name) const;
};

}
}

#endif
