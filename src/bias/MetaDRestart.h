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
#include <vector>

namespace PLMD {

class Action;
class Keywords;
class OFile;
class Value;

namespace bias {

class MetaDBiasGrid;
class MetaDHillSchedule;
class MetaDHillsFiles;

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

/// Read the bias from a grid file (METAD: GRID_RFILE, PBMETAD: one of GRID_RFILES), on the lattice of the grid given
/// in input (see MetaDBiasGrid::read). With WALKERS_MPI, a relative path is in the directory of the first replica, so
/// that all the walkers read the same file. Returns the numbers of hills of the other walkers in the grid (WALKERS_N,
/// see MetaDHillsFiles::getGridComment), or an empty vector.
std::vector<long long> readBiasGrid(Action& action, const MetaDWalkers& walkers, MetaDBiasGrid& grid, const std::string& filename,
                                    const std::string& funcl, const std::vector<Value*>& args, const std::vector<std::string>& gmin,
                                    const std::vector<std::string>& gmax, const std::vector<unsigned>& gbin, bool sparse, bool spline,
                                    bool nobounds);

/// The hills files of a bias when the action starts (METAD: one, PBMETAD: one for each family)
struct MetaDHillsInput {
  /// the hills files, with the bias to which their hills are added (see MetaDHillsFiles::setBias)
  MetaDHillsFiles* files=nullptr;
  /// the name of the hills file (FILE)
  std::string name;
  /// the numbers of hills of the other walkers in the grid read from file (see readBiasGrid)
  std::vector<long long> gridWalkerHills;
};

/// Open the hills files of all the walkers when the action starts, and read them:
/// - with RESTART, the hills of all the files are added to the bias (a missing file is an error);
/// - with RESTART and the bias read from grid files, the hills files of this walker are scanned, for the time of the
///   last hill and for the statistics of ADAPTIVE=DIFF, and the hills of the other walkers that are not in the grids
///   are added (WALKERS_N: their numbers are in the grid files; for older grid files, those up to the last hill of
///   this walker are assumed to be in the grids).
/// The time of the last hill of the files of this walker is given to schedule (with WALKERS_MPI, the same for all the
/// walkers). Then the incomplete last line of these files (an interrupted run) is removed, and they are closed, to be
/// written. With WALKERS_MPI and no WALKERS_DIR, the files with a relative path are in the directory of the first
/// replica (walkers.dir is set). To be called by all the processes of all the replicas.
void readHillsAtStart(Action& action, MetaDWalkers& walkers, const std::vector<MetaDHillsInput>& inputs, bool restartedFromGrid,
                      MetaDHillSchedule& schedule);

/// WALKERS_N: add the hills that the other walkers have written in their files since the last read
/// (the files that did not exist yet are opened, if they exist now)
void readOtherWalkersHills(Action& action, const MetaDWalkers& walkers, const std::vector<MetaDHillsFiles*>& files);

}
}

#endif
