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



namespace PLMD {
class fileParser;

class TrajectoryParser {
public:
  enum class trajfmt {
    molfile,
    xdr_xtc,
    xdr_trr,
    xyz,
    gro,
    dlp4,
    error
  };
  static trajfmt FMTfromString(std::string_view fmt);
  static std::string toString(trajfmt fmt);
private:
  std::unique_ptr<fileParser> parser;
public:
  TrajectoryParser();
  ~TrajectoryParser();
  static void registerKeywords(Keywords& keys);
  static std::vector<std::string> trajectoryOptions();
  static std::vector<std::string> getMolfilePluginsnames();

  std::optional<std::string> init(std::string_view fmt,
                                  std::string_view fname,
                                  bool useMolfile=false,
                                  int command_line_natoms=-1);
  //Driver with filename ="-"
  std::optional<std::string> init(std::string_view fmt,FILE* fileHandle);

  std::optional<std::string> readHeader(long long int &step,
                                        double &timeStep);
  std::optional<std::string> readHeader(long long int &step,
                                        float &timeStep);
  std::optional<std::string> readAtoms(int stride,
                                       bool dont_read_pbc,
                                       bool debug_pd,
                                       int pd_start,
                                       int pd_nlocal,
                                       long long int &step,
                                       double* masses,
                                       double* charges,
                                       double* coordinates,
                                       double* cell );
  std::optional<std::string> readAtoms(int stride,
                                       bool dont_read_pbc,
                                       bool debug_pd,
                                       int pd_start,
                                       int pd_nlocal,
                                       long long int &step,
                                       float* masses,
                                       float* charges,
                                       float* coordinates,
                                       float* cell );
  std::optional<std::string> readFrame(int stride,
                                       bool dont_read_pbc,
                                       bool debug_pd,
                                       int pd_start,
                                       int pd_nlocal,
                                       long long int &step,
                                       double&timeStep,
                                       double* masses,
                                       double* charges,
                                       double* coordinates,
                                       double* cell );
  std::optional<std::string> readFrame(int stride,
                                       bool dont_read_pbc,
                                       bool debug_pd,
                                       int pd_start,
                                       int pd_nlocal,
                                       long long int &step,
                                       float&timeStep,
                                       float* masses,
                                       float* charges,
                                       float* coordinates,
                                       float* cell );
  /// Return the number of atoms
  int nOfAtoms() const;
  /// Return the file pointer to the initial position (use at your own risk)
  std::optional<std::string> rewind();
};
} //namespace PLMD
#endif //__PLUMED_tools_TrajectoryParser_h
