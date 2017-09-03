/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    Copyright (c) 2017 of Haochuan Chen

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published
    by the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#ifdef __PLUMED_HAS_BOOST_SERIALIZATION
#include "cltools/CLTool.h"
#include "cltools/CLToolRegister.h"
#include "config/Config.h"
#include "core/ActionRegister.h"
#include "DRR.h"
#include "tools/Tools.h"
#include "tools/Units.h"
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

using namespace std;
using namespace PLMD::DRR;

namespace PLMD {
namespace cltools {

//+PLUMEDOC TOOLS drrtool
/*
 - Extract .grad and .count files from the binary output .drrstate
 - Merge windows

\par Examples

The following command will extract .grad and .count files.
\verbatim
plumed drrtool --extract eabf.drrstate
\endverbatim

The following command will merge windows of two .drrstate file, and output the
.grad and .count files.
\verbatim
plumed drrtool --merge win1.drrstate,win2.drrstate
\endverbatim

After getting the .grad and .count file, you can do numerical integration by
using abf_integrate tool from
https://github.com/Colvars/colvars/tree/master/colvartools
\verbatim
abf_integrate eabf.czar.grad
\endverbatim

*/
//+ENDPLUMEDOC

class drrtool : public CLTool {
public:
  static void registerKeywords(Keywords &keys);
  explicit drrtool(const CLToolOptions &co);
  int main(FILE *in, FILE *out, Communicator &pc);
  void extractdrr(const std::vector<std::string> &filename);
  void mergewindows(const std::vector<std::string> &filename);
  string description() const { return "Extract or merge the drrstate files."; }

private:
  bool verbosity;
};

PLUMED_REGISTER_CLTOOL(drrtool, "drrtool")

void drrtool::registerKeywords(Keywords &keys) {
  CLTool::registerKeywords(keys);
  keys.add("optional", "--extract", "Extract drrstate file(s)");
  keys.add("optional", "--merge", "Merge eABF windows");
  keys.addFlag("-v", false, "Verbose output");
}

drrtool::drrtool(const CLToolOptions &co) : CLTool(co) {
  inputdata = commandline;
}

int drrtool::main(FILE *in, FILE *out, Communicator &pc) {
  parseFlag("-v", verbosity);
  vector<string> stateFilesToExtract;
  bool doextract = parseVector("--extract", stateFilesToExtract);
  if (doextract) {
    extractdrr(stateFilesToExtract);
  }
  vector<string> stateFilesToMerge;
  bool domerge = parseVector("--merge", stateFilesToMerge);
  if (domerge) {
    mergewindows(stateFilesToMerge);
  }
  return 0;
}

void drrtool::extractdrr(const std::vector<std::string> &filename) {
#pragma omp parallel for
  for (size_t i = 0; i < filename.size(); ++i) {
    std::ifstream in;
    in.open(filename[i]);
    boost::archive::binary_iarchive ia(in);
    long long int step;
    std::vector<double> fict;
    std::vector<double> vfict;
    std::vector<double> vfict_laststep;
    std::vector<double> ffict;
    ABF abfgrid;
    CZAR czarestimator;
    ia >> step >> fict >> vfict >> vfict_laststep >> ffict >> abfgrid >>
        czarestimator;
    in.close();
    if (verbosity) {
      std::cout << "Dumping information of extended variables..." << '\n';
      std::cout << "Step: " << step << '\n';
      for (size_t i = 0; i < fict.size(); ++i) {
        std::cout << "Dimension[" << i + 1 << "]:\n"
                  << "  Coordinate: " << fict[i] << '\n'
                  << "  Velocity: " << vfict[i] << '\n'
                  << "  Velocity(laststep): " << vfict_laststep[i] << '\n'
                  << "  Force: " << ffict[i] << '\n';
      }
      std::cout << "Dumping counts and gradients from grids..." << '\n';
    }
    std::string outputname(filename[i]);
    const std::string suffix(".drrstate");
    outputname = outputname.substr(0, outputname.length() - suffix.length());
    if (verbosity)
      std::cout << "Writing ABF(naive) estimator files..." << '\n';
    abfgrid.writeAll(outputname);
    if (verbosity)
      std::cout << "Writing CZAR estimator files..." << '\n';
    czarestimator.writeAll(outputname);
  }
}

void drrtool::mergewindows(const std::vector<std::string> &filename) {
  if (filename.size() < 2) {
    std::cerr << "ERROR!" << std::endl;
    std::abort();
  }
  // Read grid into abfs and czars;
  std::vector<ABF> abfs;
  std::vector<CZAR> czars;
  for (size_t i = 0; i < filename.size(); ++i) {
    std::ifstream in;
    in.open(filename[i]);
    boost::archive::binary_iarchive ia(in);
    long long int step;
    std::vector<double> fict;
    std::vector<double> vfict;
    std::vector<double> vfict_laststep;
    std::vector<double> ffict;
    ABF abfgrid;
    CZAR czarestimator;
    ia >> step >> fict >> vfict >> vfict_laststep >> ffict >> abfgrid >>
        czarestimator;
    abfs.push_back(abfgrid);
    czars.push_back(czarestimator);
    in.close();
  }
  CZAR cmerged = CZAR::mergewindow(czars[0], czars[1]);
  ABF amerged = ABF::mergewindow(abfs[0], abfs[1]);
  for (size_t i = 2; i < czars.size(); ++i) {
    cmerged = CZAR::mergewindow(cmerged, czars[i]);
    amerged = ABF::mergewindow(amerged, abfs[i]);
  }
  std::string coutputname("czar.");
  std::string aoutputname("abf.");
  cmerged.writeAll(coutputname);
  amerged.writeAll(aoutputname);
}

} // End of namespace
}

#endif
