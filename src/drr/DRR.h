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
#ifndef __PLUMED_drr_DRR_h
#define __PLUMED_drr_DRR_h
// Build requirement: boost, c++11 compatible compiler.
#ifdef __PLUMED_HAS_BOOST_SERIALIZATION

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <numeric>
#include <sstream>

// boost headers for serialization
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/vector.hpp>

namespace PLMD {
namespace drr {

using std::vector;
using std::string;
using std::begin;
using std::end;

/// This class can store the minimum, maximum and bins of a dimension(axis).
class DRRAxis {
public:
  /// Default empty constructor
  DRRAxis() {
    min = max = 0.0;
    nbins = 0;
    periodic = false;
    domainMax = domainMin = 0.0;
    binWidth = 0.0;
  }
  /// Constructor using maximum value, minimum value and the number of bins(No
  /// pbc)
  DRRAxis(double l, double h, size_t n)
    : min(l), max(h), nbins(n), periodic(false), domainMax(0), domainMin(0),
      binWidth((max - min) / double(nbins)) {}
  /// PBC-aware constructor
  DRRAxis(double l, double h, size_t n, bool pbc, double dMax, double dMin)
    : min(l), max(h), nbins(n), periodic(pbc), domainMax(dMax),
      domainMin(dMin), binWidth((max - min) / double(nbins)) {}
  /// Set values
  void set(double l, double h, size_t n, bool pbc = false, double dmin = 0,
           double dmax = 0) {
    min = l;
    max = h;
    nbins = n;
    periodic = pbc;
    domainMax = dmax;
    domainMin = dmin;
    binWidth = (max - min) / nbins;
  }
  /// Set PBC data
  void setPeriodicity(double dmin, double dmax) {
    domainMax = dmax;
    domainMin = dmin;
    periodic = true;
  }
  /// Getters
  double getMin() const { return this->min; }
  double getMax() const { return this->max; }
  double getWidth() const { return binWidth; }
  double getDomainMax() const { return this->domainMax; }
  double getDomainMin() const { return this->domainMin; }
  size_t getBins() const { return this->nbins; }

  /// Check periodicity
  bool isPeriodic() const { return this->periodic; }
  /// Check real periodicity, i.e. the maximum == the domain maximum
  bool isRealPeriodic() const;

  /// Check whether x is in this axis
  bool isInBoundary(double x) const;
  /// Get an array of middle points of each bins
  vector<double> getMiddlePoints();

  /// Combine two axes if they share the same bin width.
  static DRRAxis merge(const DRRAxis &d1, const DRRAxis &d2);

  friend class DRRForceGrid;

protected:
  double min;       // Minimum value of the axis
  double max;       // Maximum value of the axis
  size_t nbins;     // Number of bins
  bool periodic;    // Periodicity
  double domainMax; // Maximum value of the CV domain
  double domainMin; // Minimum value of the CV domain
  friend class boost::serialization::access;
  /// Use boost serialization
  template <typename Archive>
  void save(Archive &ar, const unsigned int version) const {
    ar &min;
    ar &max;
    ar &nbins;
    ar &periodic;
    ar &domainMax;
    ar &domainMin;
  }
  /// Split save and load. The bin width is calculated after initialization.
  template <typename Archive>
  void load(Archive &ar, const unsigned int version) {
    ar &min;
    ar &max;
    ar &nbins;
    ar &periodic;
    ar &domainMax;
    ar &domainMin;
    binWidth = (max - min) / double(nbins);
  }
  template <typename Archive>
  void serialize(Archive &ar, const unsigned int version) {
    boost::serialization::split_member(ar, *this, version);
  }

private:
  double binWidth; // bin width
};

/// A class for collecting instantaneous forces, calculating average forces and
/// build CV histogram.
class DRRForceGrid {
public:
  /// Empty constructor
  DRRForceGrid();
  /// "Real" constructor
  /// The 2D table vector is mainly used for print grid points in grad and count
  /// file.
  /// So when use binary output we can set initializeTable to false to save
  /// memory.
  explicit DRRForceGrid(const vector<DRRAxis> &p_dimensions,
                        const string &p_suffix,
                        bool initializeTable = true);
  /// Check whether a point is in this grid
  bool isInBoundary(const vector<double> &pos) const;
  //  /// Get internal indices of a point
  //  vector<size_t> index(const vector<double> &pos) const;
  /// Get internal counts address of a point
  size_t sampleAddress(const vector<double> &pos) const;
  /// Store instantaneous forces of a point
  /// nsamples > 1 is useful for merging windows
  bool store(const vector<double> &pos, const vector<double> &f,
             unsigned long int nsamples = 1);
  /// Get accumulated forces of a point
  vector<double>
  getAccumulatedForces(const vector<double> &pos) const;
  /// Get counts of a point
  unsigned long int getCount(const vector<double> &pos,
                             bool SkipCheck = false) const;
  /// Virtual function! get gradients of a point
  /// CZAR and naive(ABF) have different gradient formulae
  virtual vector<double> getGradient(const vector<double> &pos,
                                     bool SkipCheck = false) const;
  /// Calculate divergence of the mean force field (experimental)
  double getDivergence(const vector<double> &pos) const;
  /// Calculate dln(ρ)/dz, useful for CZAR
  /// This function may be moved to CZAR class in the future
  vector<double>
  getCountsLogDerivative(const vector<double> &pos) const;
  /// Write grad file
//   void writeGrad(string filename) const;
  /// Write 1D pmf file on one dimensional occasion
  void write1DPMF(string filename) const;
  /// Write count file
//   void writeCount(string filename) const;
  /// Write necessary output file in one function (.grad and .count)
  void writeAll(const string &filename) const;
  /// Output divergence (.div) (experimental)
  void writeDivergence(const string &filename) const;
  /// merge windows
  static vector<DRRAxis> merge(const vector<DRRAxis> &dA,
                               const vector<DRRAxis> &dB);
  /// Get suffix
  string getSuffix() const { return suffix; }
  /// Set unit for .grad output
  void setOutputUnit(double unit) { outputunit = unit; }
  /// Destructor
  virtual ~DRRForceGrid() {}

protected:
  /// The output suffix appended before .grad(.czar.grad) and
  /// .count(.czar.count)
  string suffix;
  /// Number of dimensions
  size_t ndims;
  /// Store each axes
  vector<DRRAxis> dimensions;
  /// Size of samples
  size_t sampleSize;
  /// The header lines of .grad and .count files
  string headers;
  /// A table stores the middle points of all dimensions.
  /// For output in .grad and .count files
  vector<vector<double>> table;
  /// Store the average force of each bins
  vector<double> forces;
  /// Store counts of each bins
  vector<unsigned long int> samples;
  /// Only for 1D pmf output
  vector<double> endpoints;
  /// For faster indexing
  /// shifts[0] = 1, shifts[n+1] = shifts[n] * dimensions[n].nbins
  vector<size_t> shifts;
  /// For set different output units
  double outputunit;

  /// Miscellaneous helper functions
  static size_t index1D(const DRRAxis &c, double x);
  void fillTable(const vector<vector<double>> &in);

  /// Boost serialization functions
  friend class boost::serialization::access;
  template <class Archive>
  void save(Archive &ar, const unsigned int version) const {
    // Don't save all members.
    ar << suffix;
    ar << dimensions;
    ar << forces;
    ar << samples;
  }
  template <class Archive> void load(Archive &ar, const unsigned int version) {
    ar >> suffix;
    ar >> dimensions;
    ar >> forces;
    ar >> samples;
    // Restore other members.
    ndims = dimensions.size();
    sampleSize = samples.size();
    std::stringstream ss;
    ss << "# " << ndims << '\n';
    vector<vector<double>> mp(ndims);
    shifts.resize(ndims, 0);
    shifts[0] = 1;
    for (size_t i = 0; i < ndims; ++i) {
      mp[i] = dimensions[i].getMiddlePoints();
      if (i > 0) {
        shifts[i] = shifts[i - 1] * dimensions[i - 1].nbins;
      }
      ss.precision(std::numeric_limits<double>::max_digits10);
      ss << std::fixed << "# " << dimensions[i].min << ' '
         << dimensions[i].binWidth << ' ' << dimensions[i].nbins;
      if (dimensions[i].isPeriodic())
        ss << " 1" << '\n';
      else
        ss << " 0" << '\n';
    }
    fillTable(mp);
    headers = ss.str();
    outputunit = 1.0;
    // For 1D pmf
    if (ndims == 1) {
      endpoints.resize(dimensions[0].nbins + 1, 0);
      double ep = dimensions[0].min;
      double stride = dimensions[0].binWidth;
      for (auto it = begin(endpoints); it != end(endpoints); ++it) {
        (*it) = ep;
        ep += stride;
      }
    }
  }
  template <typename Archive>
  void serialize(Archive &ar, const unsigned int version) {
    boost::serialization::split_member(ar, *this, version);
  }
};

class ABF : public DRRForceGrid {
public:
  ABF() {}
  ABF(const vector<DRRAxis> &p_dimensions, const string &p_suffix,
      double fullSamples = 500.0, double maxFactor = 1.0,
      bool initializeTable = true)
    : DRRForceGrid(p_dimensions, p_suffix, initializeTable),
      mFullSamples(fullSamples), mMaxFactors(p_dimensions.size(), maxFactor) {}
  ABF(const vector<DRRAxis> &p_dimensions, const string &p_suffix,
      double fullSamples, const vector<double>& maxFactors,
      bool initializeTable = true)
    : DRRForceGrid(p_dimensions, p_suffix, initializeTable),
      mFullSamples(fullSamples), mMaxFactors(maxFactors) {}
  // Provide a setter for ABF parametres (fullsamples, maxfactor)
  void setParameters(double fullSamples, const vector<double>& maxFactors) {
    mFullSamples = fullSamples;
    mMaxFactors = maxFactors;
  }
  // Store the "instantaneous" spring force of a point and get ABF bias forces.
  bool store_getbias(const vector<double> &pos,
                     const vector<double> &f,
                     vector<double> &fbias);
  static ABF mergewindow(const ABF &aWA, const ABF &aWB);
  ~ABF() {}

private:
  // Parametres for calculate bias force
  double mFullSamples;
  vector<double> mMaxFactors;
  // Boost serialization
  friend class boost::serialization::access;
  template <typename Archive>
  void serialize(Archive &ar, const unsigned int version) {
    ar &boost::serialization::base_object<DRRForceGrid>(*this);
  }
};

class CZAR : public DRRForceGrid {
public:
  CZAR() : kbt(0) {}
  CZAR(const vector<DRRAxis> &p_dimensions, const string &p_suffix,
       double p_kbt, bool initializeTable = true)
    : DRRForceGrid(p_dimensions, p_suffix, initializeTable), kbt(p_kbt) {}
  vector<double> getGradient(const vector<double> &pos,
                             bool SkipCheck = false) const;
  double getkbt() const { return kbt; }
  void setkbt(double p_kbt) { kbt = p_kbt; }
  static CZAR mergewindow(const CZAR &cWA, const CZAR &cWB);
  void writeZCount(const string &filename) const;
  ~CZAR() {}

private:
  double kbt;
  friend class boost::serialization::access;
  template <typename Archive>
  void serialize(Archive &ar, const unsigned int version) {
    ar &boost::serialization::base_object<DRRForceGrid>(*this);
    ar &kbt;
  }
};
}
}

#endif
#endif
