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
#ifndef __PLUMED_tools_GridExtension_h
#define __PLUMED_tools_GridExtension_h

#include <cstddef>
#include <string>
#include <vector>

namespace PLMD {

/// Plans how a grid on a lattice grows so that it contains a range of points: the rules shared by the grid of
/// METAD/PBMETAD (bias/MetaDGrid) and that of KDE (gridtools/KDEGridTools).
///
/// Along each direction the caller gives the lattice indices of the first and the last point that the grid must
/// contain (relative to the current first point, so the first can be negative and the last beyond the current last
/// point). A direction can be:
/// - fixed: not changed (periodic, not extendable, ...);
/// - placed: the grid becomes exactly the range from the first to the last point (also fewer points than now);
/// - extended: the points that are missing are added, and the grid grows by at least a quarter of its current number
///   of points at that end, so that it is not extended at every step. If the grown grid has more points than the
///   limit (GRID_SIZE_LIMIT), only the missing points are added; if even that is too many, the plan is too large.
///
/// The planner only computes numbers of points: the caller applies them (GridBase::resize,
/// GridCoordinatesObject::extend) and writes its own errors and log lines.
class GridExtension {
public:
  enum class Mode {fixed, place, extend};
  struct Plan {
    /// some points are added (or removed, when a direction is placed)
    bool changed=false;
    /// points to add at the lower and upper end of each direction (negative: points removed, only when placed)
    std::vector<long long> lower, upper;
    /// number of points of the planned grid
    double npoints=0.0;
    /// the grid has more points than the limit even without growth (lower and upper are the minimal ones)
    bool toolarge=false;
    /// the first direction whose range is not finite or too far from the grid (npos if none): nothing is planned
    std::size_t nonfinite=std::string::npos;
  };
  /// an error message if limit is not a valid GRID_SIZE_LIMIT (grid points are indexed with unsigned integers),
  /// an empty string otherwise
  static std::string checkSizeLimit(double limit);
  /// set the maximum number of points; npoints are the points of the initial grid along each direction: an initial
  /// grid larger than the limit raises it to twice its size (but not beyond the largest unsigned).
  /// Returns true if the limit has been raised.
  bool setSizeLimit(double limit, const std::vector<std::size_t>& npoints);
  double getSizeLimit() const {
    return sizelimit_;
  }
  /// npoints: the current number of points along each direction (bins+1 for a direction that is not periodic);
  /// first, last: lattice indices of the first and last point needed (ignored for fixed directions)
  Plan plan(const std::vector<std::size_t>& npoints, const std::vector<Mode>& mode,
            const std::vector<double>& first, const std::vector<double>& last) const;
private:
  double sizelimit_=1.0e8;
};

}

#endif
