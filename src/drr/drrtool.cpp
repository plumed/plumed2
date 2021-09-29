/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    Copyright (c) 2017 of Haochuan Chen (excluding colvar_UIestimator.h)
    Copyright (c) 2017 of Haohao Fu (colvar_UIestimator.h)

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

using namespace PLMD;
using namespace cltools;

namespace PLMD {
namespace drr {

//+PLUMEDOC EABFMOD_TOOLS drr_tool
/*
 - Extract .grad and .count files from the binary output .drrstate
 - Merge windows

\par Examples

The following command will extract .grad and .count files.
\verbatim
plumed drr_tool --extract eabf.drrstate
\endverbatim

The following command will merge windows of two .drrstate file, and output the
.grad and .count files.
\verbatim
plumed drr_tool --merge win1.drrstate,win2.drrstate
\endverbatim

After getting the .grad and .count file, you can do numerical integration by
using abf_integrate tool from
https://github.com/Colvars/colvars/tree/master/colvartools
\verbatim
abf_integrate eabf.czar.grad
\endverbatim
\note
The abf_integrate in colvartools is in kcal/mol, so it may be better to use --units kcal/mol when running drr_tool

*/
//+ENDPLUMEDOC

using std::vector;
using std::string;

class drrtool : public CLTool {
public:
  static void registerKeywords(Keywords &keys);
  explicit drrtool(const CLToolOptions &co);
  int main(FILE *in, FILE *out, Communicator &pc);
  void extractdrr(const vector<string> &filename);
  void mergewindows(const vector<string> &filename, string outputname);
  void calcDivergence(const vector<string> &filename);
  string description() const { return "Extract or merge the drrstate files."; }

private:
  bool verbosity;
  Units units;
  const string suffix{".drrstate"};
};

PLUMED_REGISTER_CLTOOL(drrtool, "drr_tool")

void drrtool::registerKeywords(Keywords &keys) {
  CLTool::registerKeywords(keys);
  keys.add("optional", "--extract", "Extract drrstate file(s)");
  keys.add("optional", "--merge", "Merge eABF windows");
  keys.add("optional", "--merge_output", "The output filename of the merged result");
  keys.add("optional", "--divergence", "Calculate divergence of gradient field (experimental)");
  keys.add("compulsory","--units","kj/mol","the units of energy can be kj/mol, kcal/mol, j/mol, eV or the conversion factor from kj/mol");
  keys.addFlag("-v", false, "Verbose output");
}

drrtool::drrtool(const CLToolOptions &co) : CLTool(co) {
  inputdata = commandline;
  verbosity = false;
}

int drrtool::main(FILE *in, FILE *out, Communicator &pc) {
  parseFlag("-v", verbosity);
  vector<string> stateFilesToExtract;
  string unitname;
  parse("--units",unitname);
  units.setEnergy( unitname );
  bool doextract = parseVector("--extract", stateFilesToExtract);
  if (doextract) {
    extractdrr(stateFilesToExtract);
  }
  vector<string> stateFilesToMerge;
  bool domerge = parseVector("--merge", stateFilesToMerge);
  if (domerge) {
    string merge_outputname;
    parse("--merge_output", merge_outputname);
    mergewindows(stateFilesToMerge, merge_outputname);
  }
  vector<string> stateFilesToDivergence;
  bool dodivergence = parseVector("--divergence", stateFilesToDivergence);
  if (dodivergence) {
    calcDivergence(stateFilesToDivergence);
  }
  return 0;
}

void drrtool::extractdrr(const vector<string> &filename) {
  #pragma omp parallel for
  for (size_t j = 0; j < filename.size(); ++j) {
    std::ifstream in;
    in.open(filename[j]);
    boost::archive::binary_iarchive ia(in);
    long long int step;
    vector<double> fict;
    vector<double> vfict;
    vector<double> vfict_laststep;
    vector<double> ffict;
    ABF abfgrid;
    CZAR czarestimator;
    ia >> step >> fict >> vfict >> vfict_laststep >> ffict >> abfgrid >>
       czarestimator;
    in.close();
    abfgrid.setOutputUnit(units.getEnergy());
    czarestimator.setOutputUnit(units.getEnergy());
    if (verbosity) {
      std::cout << "Output units factor: " << units.getEnergy() << '\n';
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
    string outputname(filename[j]);
    outputname = outputname.substr(0, outputname.length() - suffix.length());
    if (verbosity)
      std::cout << "Writing ABF(naive) estimator files..." << '\n';
    abfgrid.writeAll(outputname);
    if (verbosity)
      std::cout << "Writing CZAR estimator files..." << '\n';
    czarestimator.writeAll(outputname);
    czarestimator.writeZCount(outputname);
  }
}

void drrtool::mergewindows(const vector<string> &filename, string outputname) {
  if (filename.size() < 2) {
    std::cerr << "ERROR! You need at least two .drrstate file to merge windows!" << std::endl;
    std::abort();
  }
  // Read grid into abfs and czars;
  vector<ABF> abfs;
  vector<CZAR> czars;
  for (auto it_fn = filename.begin(); it_fn != filename.end(); ++it_fn) {
    std::ifstream in;
    in.open((*it_fn));
    boost::archive::binary_iarchive ia(in);
    long long int step;
    vector<double> fict;
    vector<double> vfict;
    vector<double> vfict_laststep;
    vector<double> ffict;
    ABF abfgrid;
    CZAR czarestimator;
    ia >> step >> fict >> vfict >> vfict_laststep >> ffict >> abfgrid >>
       czarestimator;
    abfgrid.setOutputUnit(units.getEnergy());
    czarestimator.setOutputUnit(units.getEnergy());
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
  if (outputname.empty()) {
    // Generate new file name for merged grad and count
    vector<string> tmp_name = filename;
    std::transform(std::begin(tmp_name), std::end(tmp_name), std::begin(tmp_name),
    [&](string s) {return s.substr(0, s.find(suffix));});
    outputname = std::accumulate(std::begin(tmp_name), std::end(tmp_name), string(""),
    [](const string & a, const string & b) {return a + b + "+";});
    outputname = outputname.substr(0, outputname.size() - 1);
    std::cerr << "You have not specified an output filename for the merged"
              << " result, so the default name \"" + outputname
              << "\" is used here, which may yield unexpected behavior.\n";
  }
  cmerged.writeAll(outputname);
  cmerged.writeZCount(outputname);
  amerged.writeAll(outputname);
}

void drrtool::calcDivergence(const vector<string> &filename) {
  #pragma omp parallel for
  for (size_t j = 0; j < filename.size(); ++j) {
    std::ifstream in;
    in.open(filename[j]);
    boost::archive::binary_iarchive ia(in);
    long long int step;
    vector<double> fict;
    vector<double> vfict;
    vector<double> vfict_laststep;
    vector<double> ffict;
    ABF abfgrid;
    CZAR czarestimator;
    ia >> step >> fict >> vfict >> vfict_laststep >> ffict >> abfgrid >>
       czarestimator;
    in.close();
    abfgrid.setOutputUnit(units.getEnergy());
    czarestimator.setOutputUnit(units.getEnergy());
    if (verbosity) {
      std::cout << "Output units factor: " << units.getEnergy() << '\n';
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
    string outputname(filename[j]);
    outputname = outputname.substr(0, outputname.length() - suffix.length());
    abfgrid.writeDivergence(outputname);
    czarestimator.writeDivergence(outputname);
  }
}

} // End of namespace
}

#endif
