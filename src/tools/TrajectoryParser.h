/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2025 The plumed team
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

#ifndef __PLUMED_tools_TrajectoryParser_h
#define __PLUMED_tools_TrajectoryParser_h
#include "Keywords.h"
#include "Tools.h"

#include <memory>
#include <optional>
#include <string_view>
#include <utility>
#include <vector>



namespace PLMD {
class fileParser;
struct RawFrame;

/// A frame of a trajectory, as read by TrajectoryParser.
/// Lengths are in nm and times in ps (the readers convert from the units of the file);
/// masses and charges are in the units of the file.
/// The caller sets the inputs and provides the buffers, sized for
/// TrajectoryParser::nOfAtoms() atoms (known after readHeader()).
template <typename real>
struct TrajectoryFrame {
  /// input: the trajectory stride of the caller.
  /// With 0, the step stored in xtc and trr frames replaces `step`
  int stride=1;
  /// input: false when the box is given in another way (the driver's --box).
  /// The dlp4 reader then expects no cell vectors in the file
  bool readBox=true;
  /// input: only the positions of the atoms from `first` to `first+count-1` are copied
  /// (the driver uses it to test particle decomposition); with count<0 all the atoms are copied
  int first=0;
  /// input: see `first`
  int count=-1;
  /// output buffer: 3*natoms values, x y z of each atom
  real* positions=nullptr;
  /// output buffer: 9 values, the cell vectors a, b and c, one per row
  real* box=nullptr;
  /// output buffer: natoms values, written only if the file stores the masses (see hasMassesAndCharges)
  real* masses=nullptr;
  /// output buffer: natoms values, written only if the file stores the charges (see hasMassesAndCharges)
  real* charges=nullptr;
  /// input and output: the step counter of the caller.
  /// The step stored in the file replaces it for dlp4 (always, in readHeader())
  /// and for xtc, trr and the formats read by chemfiles that store it (only with stride 0, in readAtoms())
  long long step=0;
  /// output of readHeader(): the time step of the simulation, for the formats that store it (dlp4)
  std::optional<double> timestep;
  /// output of readAtoms(): the time of the frame, for the formats that store it (xtc, trr, some of the formats read by chemfiles)
  std::optional<double> time;
  /// output of readAtoms(): true if masses and charges were read from the file (dlp4, and TPR and LAMMPS Data with chemfiles)
  bool hasMassesAndCharges=false;
  /// output of readAtoms(): set together with "EOF" when the last frame of the file is incomplete
  /// (a truncated xtc, trr or dcd file, as written by a simulation that stopped while writing)
  std::optional<std::string> warning;
};

/// Reads trajectories for the driver and the benchmark.
/// The functions that read return std::nullopt on success, "EOF" at the end of the file,
/// and an error message otherwise.
class TrajectoryParser {
public:
  enum class trajfmt {
    xdr_xtc,
    xdr_trr,
    xyz,
    gro,
    dlp4,
    pdb,
    dcd,
    /// any format read by chemfiles (--itraj)
    chemfiles,
    error
  };
  static trajfmt FMTfromString(std::string_view fmt);
  static std::string toString(trajfmt fmt);
private:
  std::unique_ptr<fileParser> parser;
  /// the frame as read by the parser, before the copy to the buffers of the caller
  std::unique_ptr<RawFrame> raw;
public:
  TrajectoryParser();
  ~TrajectoryParser();
  static void registerKeywords(Keywords& keys);
  /// The formats, each read with the option --i<format> ("traj" is any format read by chemfiles)
  static std::vector<std::string> trajectoryOptions();
  /// True if PLUMED was compiled with chemfiles (--itraj)
  static bool hasChemfiles();
  /// The options of the molfile plugins, removed in v2.11, that are still accepted for one release:
  /// pairs {option, format}, e.g. {"--mf_xtc","xtc"} (the same as --ixtc)
  static const std::vector<std::pair<std::string,std::string>> & deprecatedOptions();
  /// The options of the molfile plugins, removed in v2.11, whose formats are not read anymore
  /// (e.g. "--mf_crd"). registerKeywords() registers them as hidden options: the driver and the benchmark
  /// stop with removedOptionMessage() instead of reporting an unknown option
  static const std::vector<std::string> & removedOptions();
  static std::string removedOptionMessage(std::string_view option);

  /// Open the trajectory file fname in the format fmt (e.g. "xtc", see trajectoryOptions()).
  /// With fmt "traj" (chemfiles), format is the format of the file as named by chemfiles,
  /// or empty to guess it from the extension of the file
  std::optional<std::string> init(std::string_view fmt,
                                  std::string_view fname,
                                  std::string_view format="");
  //Driver with filename ="-"
  std::optional<std::string> init(std::string_view fmt,FILE* fileHandle);

  /// Read the header of the next frame: after it nOfAtoms() is the number of atoms of the frame
  template <typename real>
  std::optional<std::string> readHeader(TrajectoryFrame<real>& frame);
  /// Read the atoms of the frame whose header was read last
  template <typename real>
  std::optional<std::string> readAtoms(TrajectoryFrame<real>& frame);
  /// readHeader() followed by readAtoms(), when the buffers are already sized
  template <typename real>
  std::optional<std::string> readFrame(TrajectoryFrame<real>& frame);
  /// Return the number of atoms
  int nOfAtoms() const;
  /// Return the file pointer to the initial position (use at your own risk)
  std::optional<std::string> rewind();
};
} //namespace PLMD
#endif //__PLUMED_tools_TrajectoryParser_h
