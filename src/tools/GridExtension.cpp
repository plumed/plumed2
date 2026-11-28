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
#include "GridExtension.h"
#include "Exception.h"
#include "Tools.h"
#include <algorithm>
#include <cmath>
#include <limits>

namespace PLMD {

std::string GridExtension::checkSizeLimit(double limit) {
  if(!(limit>=1.0)) {
    return "GRID_SIZE_LIMIT should be at least 1";
  }
  if(limit>static_cast<double>(std::numeric_limits<unsigned>::max())) {
    return "GRID_SIZE_LIMIT cannot be larger than "+Tools::boundString(std::numeric_limits<unsigned>::max());
  }
  return "";
}

bool GridExtension::setSizeLimit(double limit, const std::vector<std::size_t>& npoints) {
  sizelimit_=limit;
  double size=1.0;
  for(const auto n : npoints) {
    size*=n;
  }
  if(size>sizelimit_) {
    sizelimit_=std::min(2.0*size,static_cast<double>(std::numeric_limits<unsigned>::max()));
    return true;
  }
  return false;
}

GridExtension::Plan GridExtension::plan(const std::vector<std::size_t>& npoints, const std::vector<Mode>& mode,
                                        const std::vector<double>& first, const std::vector<double>& last) const {
  const std::size_t ndim=npoints.size();
  plumed_assert(mode.size()==ndim && first.size()==ndim && last.size()==ndim);
  Plan p;
  // the minimal change (lower, upper) and the one with growth (lowergrow, uppergrow)
  std::vector<long long> lower(ndim,0), upper(ndim,0), lowergrow(ndim,0), uppergrow(ndim,0);
  for(std::size_t i=0; i<ndim; ++i) {
    if(mode[i]==Mode::fixed) {
      continue;
    }
    const double lo=first[i], hi=last[i];
    if(!std::isfinite(lo) || !std::isfinite(hi) || std::fabs(lo)>1e15 || std::fabs(hi)>1e15) {
      p.nonfinite=i;
      p.lower.assign(ndim,0);
      p.upper.assign(ndim,0);
      return p;
    }
    // index of the current last point
    const double end=static_cast<double>(npoints[i])-1.0;
    if(mode[i]==Mode::place) {
      lower[i]=lowergrow[i]=static_cast<long long>(-lo);
      upper[i]=uppergrow[i]=static_cast<long long>(hi-end);
      p.changed=true;
    } else {
      // geometric growth: at least a quarter of the current number of points
      const long long grow=static_cast<long long>(std::ceil(0.25*npoints[i]));
      if(lo<0.0) {
        lower[i]=static_cast<long long>(-lo);
        lowergrow[i]=std::max(lower[i],grow);
        p.changed=true;
      }
      if(hi>end) {
        upper[i]=static_cast<long long>(hi-end);
        uppergrow[i]=std::max(upper[i],grow);
        p.changed=true;
      }
    }
  }
  auto size=[&](const std::vector<long long>& lo,const std::vector<long long>& up) {
    double n=1.0;
    for(std::size_t i=0; i<ndim; ++i) {
      n*=static_cast<double>(npoints[i])+lo[i]+up[i];
    }
    return n;
  };
  if(p.changed && size(lowergrow,uppergrow)>sizelimit_) {
    lowergrow=lower;
    uppergrow=upper;
    p.toolarge=size(lower,upper)>sizelimit_;
  }
  p.npoints=size(lowergrow,uppergrow);
  p.lower.swap(lowergrow);
  p.upper.swap(uppergrow);
  return p;
}

}
