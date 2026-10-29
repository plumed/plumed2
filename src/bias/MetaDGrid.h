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
#ifndef __PLUMED_bias_MetaDGrid_h
#define __PLUMED_bias_MetaDGrid_h

#include "MetaDKernel.h"
#include "tools/Grid.h"
#include <memory>
#include <string>
#include <vector>

namespace PLMD {

class Action;
class IFile;
class Value;

namespace bias {

/// The grid on which METAD and PBMETAD accumulate their hills.
///
/// The grid is extended along the variables that are not periodic whenever a hill does not fit in it, so that it
/// always contains all the points where the hills are not zero (plus one point on each side): the bias is zero
/// outside the grid. The spacing does not change and the new points are on the same lattice.
/// With INTERVAL (one dimensional only) the bias beyond the limits is constant: the grid contains the interval and
/// at least one point beyond each limit, and the bias beyond the grid is the value at its ends.
class MetaDBiasGrid {
  /// the action that owns the grid (errors, log)
  Action& action_;
  std::unique_ptr<GridBase> grid_;
  bool sparse_=false;
  /// maximum number of points (GRID_SIZE_LIMIT)
  double sizelimit_=1.0e8;
  /// no range was given and there are no hills yet: the grid is placed around the first one
  bool pristine_=false;
  /// a line of the log is not complete (an extension starts a new line)
  bool midline_=false;
  /// INTERVAL (one dimensional)
  bool interval_=false;
  double lower_=0.0;
  double upper_=0.0;
  /// with INTERVAL, true if the bias for cv is the constant value at an end of the grid
  bool isIntervalBeyondGrid(double cv) const;
  /// with INTERVAL, the bias for a CV beyond the ends of the grid
  double getIntervalBiasOutsideGrid(double cv) const;
  /// add the hill to the points within its support (no INTERVAL)
  void addHillToPoints(const MetaDGaussian& hill, const MetaDKernel& kernel);
public:
  /// The number of bins from GRID_MIN, GRID_MAX, GRID_BIN and GRID_SPACING (empty if not given): with neither GRID_BIN
  /// nor GRID_SPACING the spacing is 1/5 of SIGMA (of SIGMA_MIN with ADAPTIVE). With GRID_SPACING alone the grid is
  /// placed around the first hill: gmin, gmax and gbin are set to a first range and true is returned.
  static bool setupBins(Action& action, const std::vector<Value*>& args, std::vector<std::string>& gmin, std::vector<std::string>& gmax,
                        std::vector<unsigned>& gbin, std::vector<double>& gspacing, bool adaptive, const std::vector<double>& sigma0,
                        const std::vector<double>& sigma0min);
  /// error if GRID_SIZE_LIMIT is not valid
  static void checkSizeLimit(Action& action, double limit);
  explicit MetaDBiasGrid(Action& action);
  /// a new grid; nobounds: no initial range was given (GRID_SPACING only), the grid is placed around the first hill
  void create(const std::string& funcl, const std::vector<Value*>& args, const std::vector<std::string>& gmin,
              const std::vector<std::string>& gmax, const std::vector<unsigned>& gbin, bool sparse, bool spline, bool nobounds);
  /// read a grid from a file: it must be on the lattice of the grid given in input (gmin, gmax, gbin), with any range
  void read(IFile& file, const std::string& filename, const std::string& funcl, const std::vector<Value*>& args,
            const std::vector<std::string>& gmin, const std::vector<std::string>& gmax, const std::vector<unsigned>& gbin,
            bool sparse, bool spline, bool nobounds);
  /// maximum number of points (an initial grid that is larger raises the limit to twice its size)
  void setSizeLimit(double limit);
  /// INTERVAL: the grid is extended to contain the interval and at least one point beyond each limit
  void setInterval(double lower, double upper, double time);
  /// extend the grid, if needed, so that it contains all the points within halfwidth[i] from x along each
  /// variable i that is not periodic, plus one point on each side
  void extend(const std::vector<double>& x, const std::vector<double>& halfwidth, double time);
  /// add a hill (the grid is extended if needed)
  void addHill(const MetaDGaussian& hill, const MetaDKernel& kernel, double time);
  /// the bias (zero outside the grid, or the constant value beyond the limits with INTERVAL)
  double getBias(const std::vector<double>& cv) const;
  /// the bias and its derivatives
  double getBiasAndDerivatives(const std::vector<double>& cv, std::vector<double>& der) const;
  /// a line of the log has been started and is not complete (e.g. while reading hills)
  void setMidline(bool midline) {
    midline_=midline;
  }
  bool isSparse() const {
    return sparse_;
  }
  GridBase& getGrid() {
    return *grid_;
  }
  const GridBase& getGrid() const {
    return *grid_;
  }
};

}
}

#endif
