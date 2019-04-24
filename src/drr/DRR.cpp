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
#include "DRR.h"

namespace PLMD {
namespace drr {

using std::vector;
using std::string;
using std::begin;
using std::end;

bool DRRAxis::isInBoundary(double x) const {
  if (x < min || x > max)
    return false;
  else
    return true;
}

bool DRRAxis::isRealPeriodic() const {
  if (periodic == true) {
    if (std::abs(domainMax - max) < binWidth &&
        std::abs(domainMin - min) < binWidth) {
      return true;
    } else {
      return false;
    }
  } else {
    return false;
  }
}

DRRAxis DRRAxis::merge(const DRRAxis &d1, const DRRAxis &d2) {
  const double newmin = std::min(d1.min, d2.min);
  const double newmax = std::max(d1.max, d2.max);
  const double newWidth = d1.binWidth;
  const size_t newbins = size_t(std::nearbyint((newmax - newmin) / newWidth));
  const bool newpbc = d1.periodic;
  const double newdmin = std::min(d1.domainMin, d2.domainMin);
  const double newdmax = std::max(d1.domainMax, d2.domainMax);
  DRRAxis result(newmin, newmax, newbins, newpbc, newdmin, newdmax);
  return result;
}

vector<double> DRRAxis::getMiddlePoints() {
  vector<double> result(nbins, 0);
  const double width = binWidth;
  double temp = min - width / 2;
  std::generate(begin(result), end(result), [&temp, &width]() {
    temp += width;
    return temp;
  });
  return result;
}

size_t DRRForceGrid::index1D(const DRRAxis &c, double x) {
  size_t idx = size_t(std::floor((x - c.min) / c.binWidth));
  idx = (idx == c.nbins) ? (c.nbins - 1) : idx;
  return idx;
}

void DRRForceGrid::fillTable(const vector<vector<double>> &in) {
  table.resize(ndims, vector<double>(sampleSize, 0));
  for (size_t i = 0; i < ndims; ++i) {
    size_t repeatAll = 1, repeatOne = 1;
    for (size_t j = i + 1; j < ndims; ++j)
      repeatOne *= in[j].size();
    for (size_t j = 0; j < i; ++j)
      repeatAll *= in[j].size();
    size_t in_i_sz = in[i].size();
    for (size_t l = 0; l < in_i_sz; ++l)
      std::fill_n(begin(table[i]) + l * repeatOne, repeatOne, in[i][l]);
    for (size_t k = 0; k < repeatAll - 1; ++k)
      std::copy_n(begin(table[i]), repeatOne * in_i_sz,
                  begin(table[i]) + repeatOne * in_i_sz * (k + 1));
  }
}

DRRForceGrid::DRRForceGrid()
  : suffix(""), ndims(0), dimensions(0), sampleSize(0),
    headers(""), table(0), forces(0), samples(0), endpoints(0), shifts(0),
    outputunit(1.0) {}

DRRForceGrid::DRRForceGrid(const vector<DRRAxis> &p_dimensions,
                           const string &p_suffix, bool initializeTable)
  : suffix(p_suffix), ndims(p_dimensions.size()), dimensions(p_dimensions) {
  sampleSize = 1;
  vector<vector<double>> mp(ndims);
  std::stringstream ss;
  ss << "# " << ndims << '\n';
  shifts.resize(ndims, 0);
  shifts[0] = 1;
  for (size_t i = 0; i < ndims; ++i) {
    sampleSize = dimensions[i].nbins * sampleSize;
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
  headers = ss.str();
  if (initializeTable)
    fillTable(mp);
  forces.resize(sampleSize * ndims, 0.0);
  samples.resize(sampleSize, 0);
  outputunit = 1.0;
  // For 1D pmf
  if (ndims == 1) {
    endpoints.resize(dimensions[0].nbins + 1, 0);
    double ep = dimensions[0].min;
    double stride = dimensions[0].binWidth;
    for (auto &i : endpoints) {
      i = ep;
      ep += stride;
    }
  }
}

bool DRRForceGrid::isInBoundary(const vector<double> &pos) const {
  bool result = true;
  for (size_t i = 0; i < ndims; ++i) {
    if (pos[i] < dimensions[i].min || pos[i] > dimensions[i].max)
      return false;
  }
  return result;
}

size_t DRRForceGrid::sampleAddress(const vector<double> &pos) const {
  size_t saddr = 0;
  for (size_t i = 0; i < ndims; ++i) {
    saddr += shifts[i] * index1D(dimensions[i], pos[i]);
  }
  return saddr;
}

bool DRRForceGrid::store(const vector<double> &pos, const vector<double> &f,
                         unsigned long int nsamples) {
  if (isInBoundary(pos)) {
    if (nsamples == 0)
      return true;
    const size_t baseaddr = sampleAddress(pos) * ndims;
    samples[baseaddr / ndims] += nsamples;
    auto it_fa = begin(forces) + baseaddr;
    std::transform(begin(f), end(f), it_fa, it_fa, std::plus<double>());
    return true;
  } else {
    return false;
  }
}

vector<DRRAxis> DRRForceGrid::merge(const vector<DRRAxis> &dA,
                                    const vector<DRRAxis> &dB) {
  vector<DRRAxis> dR(dA.size());
  std::transform(begin(dA), end(dA), begin(dB), begin(dR), DRRAxis::merge);
  return dR;
}

vector<double>
DRRForceGrid::getAccumulatedForces(const vector<double> &pos) const {
  vector<double> result(ndims, 0);
  if (!isInBoundary(pos))
    return result;
  const size_t baseaddr = sampleAddress(pos) * ndims;
  std::copy(begin(forces) + baseaddr, begin(forces) + baseaddr + ndims,
            begin(result));
  return result;
}

unsigned long int DRRForceGrid::getCount(const vector<double> &pos,
    bool SkipCheck) const {
  if (!SkipCheck) {
    if (!isInBoundary(pos)) {
      return 0;
    }
  }
  return samples[sampleAddress(pos)];
}

vector<double> DRRForceGrid::getGradient(const vector<double> &pos,
    bool SkipCheck) const {
  vector<double> result(ndims, 0);
  if (!SkipCheck) {
    if (!isInBoundary(pos)) {
      return result;
    }
  }
  const size_t baseaddr = sampleAddress(pos);
  const unsigned long int &count = samples[baseaddr];
  if (count == 0)
    return result;
  auto it_fa = begin(forces) + baseaddr * ndims;
  std::transform(it_fa, it_fa + ndims, begin(result),
  [&count](double fa) { return (-1.0) * fa / count; });
  return result;
}

double DRRForceGrid::getDivergence(const vector<double> &pos) const {
  double div = 0.0;
  vector<double> grad_deriv(ndims, 0.0);
  if (!isInBoundary(pos)) {
    return div;
  }
  const size_t force_addr = sampleAddress(pos) * ndims;
  vector<double> grad = getGradient(pos);
  for (size_t i = 0; i < ndims; ++i) {
    const double binWidth = dimensions[i].binWidth;
    vector<double> first = pos;
    first[i] = dimensions[i].min + binWidth * 0.5;
    vector<double> last = pos;
    last[i] = dimensions[i].max - binWidth * 0.5;
    const size_t force_addr_first = sampleAddress(first) * ndims;
    const size_t force_addr_last = sampleAddress(last) * ndims;
    if (force_addr == force_addr_first) {
      if (dimensions[i].isRealPeriodic() == true) {
        vector<double> next = pos;
        next[i] += binWidth;
        const vector<double> grad_next = getGradient(next);
        const vector<double> grad_prev = getGradient(last);
        grad_deriv[i] = (grad_next[i] - grad_prev[i]) / (2.0 * binWidth);
      } else {
        vector<double> next = pos;
        next[i] += binWidth;
        vector<double> next_2 = next;
        next_2[i] += binWidth;
        const vector<double> grad_next = getGradient(next);
        const vector<double> grad_next_2 = getGradient(next_2);
        grad_deriv[i] =
          (grad_next_2[i] * -1.0 + grad_next[i] * 4.0 - grad[i] * 3.0) /
          (2.0 * binWidth);
      }
    } else if (force_addr == force_addr_last) {
      if (dimensions[i].isRealPeriodic() == true) {
        vector<double> prev = pos;
        prev[i] -= binWidth;
        const vector<double> grad_next = getGradient(first);
        const vector<double> grad_prev = getGradient(prev);
        grad_deriv[i] = (grad_next[i] - grad_prev[i]) / (2.0 * binWidth);
      } else {
        vector<double> prev = pos;
        prev[i] -= binWidth;
        vector<double> prev_2 = prev;
        prev_2[i] -= binWidth;
        const vector<double> grad_prev = getGradient(prev);
        const vector<double> grad_prev_2 = getGradient(prev_2);
        grad_deriv[i] =
          (grad[i] * 3.0 - grad_prev[i] * 4.0 + grad_prev_2[i] * 1.0) /
          (2.0 * binWidth);
      }
    } else {
      vector<double> prev = pos;
      prev[i] -= binWidth;
      vector<double> next = pos;
      next[i] += binWidth;
      const vector<double> grad_next = getGradient(next);
      const vector<double> grad_prev = getGradient(prev);
      grad_deriv[i] = (grad_next[i] - grad_prev[i]) / (2.0 * binWidth);
    }
  }
  div = std::accumulate(begin(grad_deriv), end(grad_deriv), 0.0);
  return div;
}

vector<double>
DRRForceGrid::getCountsLogDerivative(const vector<double> &pos) const {
  const size_t addr = sampleAddress(pos);
  const unsigned long int count_this = samples[addr];
  vector<double> result(ndims, 0);
  for (size_t i = 0; i < ndims; ++i) {
    const double binWidth = dimensions[i].binWidth;
    const size_t addr_first =
      addr - shifts[i] * index1D(dimensions[i], pos[i]) + 0;
    const size_t addr_last = addr_first + shifts[i] * (dimensions[i].nbins - 1);
    if (addr == addr_first) {
      if (dimensions[i].isRealPeriodic() == true) {
        const unsigned long int &count_next = samples[addr + shifts[i]];
        const unsigned long int &count_prev = samples[addr_last];
        if (count_next != 0 && count_prev != 0)
          result[i] =
            (std::log(count_next) - std::log(count_prev)) / (2 * binWidth);
      } else {
        const unsigned long int &count_next = samples[addr + shifts[i]];
        const unsigned long int &count_next2 = samples[addr + shifts[i] * 2];
        if (count_next != 0 && count_this != 0 && count_next2 != 0)
          result[i] =
            (std::log(count_next2) * (-1.0) + std::log(count_next) * 4.0 -
             std::log(count_this) * 3.0) /
            (2.0 * binWidth);
      }
    } else if (addr == addr_last) {
      if (dimensions[i].isRealPeriodic() == true) {
        const unsigned long int &count_prev = samples[addr - shifts[i]];
        const unsigned long int &count_next = samples[addr_first];
        if (count_next != 0 && count_prev != 0)
          result[i] =
            (std::log(count_next) - std::log(count_prev)) / (2 * binWidth);
      } else {
        const unsigned long int &count_prev = samples[addr - shifts[i]];
        const unsigned long int &count_prev2 = samples[addr - shifts[i] * 2];
        if (count_prev != 0 && count_this != 0 && count_prev2 != 0)
          result[i] = (std::log(count_this) * 3.0 - std::log(count_prev) * 4.0 +
                       std::log(count_prev2)) /
                      (2.0 * binWidth);
      }
    } else {
      const unsigned long int &count_prev = samples[addr - shifts[i]];
      const unsigned long int &count_next = samples[addr + shifts[i]];
      if (count_next != 0 && count_prev != 0)
        result[i] =
          (std::log(count_next) - std::log(count_prev)) / (2 * binWidth);
    }
  }
  return result;
}

void DRRForceGrid::write1DPMF(string filename) const {
  filename += suffix + ".pmf";
  FILE *ppmf;
  ppmf = fopen(filename.c_str(), "w");
  const double w = dimensions[0].binWidth;
  double pmf = 0;
  fprintf(ppmf, "%.9f %.9f\n", endpoints[0], pmf);
  for (size_t i = 0; i < dimensions[0].nbins; ++i) {
    vector<double> pos(1, 0);
    pos[0] = table[0][i];
    const vector<double> f = getGradient(pos, true);
    pmf += f[0] * w / outputunit;
    fprintf(ppmf, "%.9f %.9f\n", endpoints[i + 1], pmf);
  }
  fclose(ppmf);
}

void DRRForceGrid::writeAll(const string &filename) const {
  string countname = filename + suffix + ".count";
  string gradname = filename + suffix + ".grad";
  vector<double> pos(ndims, 0);
  FILE *pGrad, *pCount;
  pGrad = fopen(gradname.c_str(), "w");
  pCount = fopen(countname.c_str(), "w");
  char *buffer1, *buffer2;
  buffer1 = (char *)malloc((sizeof(double)) * sampleSize * ndims);
  buffer2 = (char *)malloc((sizeof(double)) * sampleSize * ndims);
  setvbuf(pGrad, buffer1, _IOFBF, (sizeof(double)) * sampleSize * ndims);
  setvbuf(pCount, buffer2, _IOFBF, (sizeof(double)) * sampleSize * ndims);
  fwrite(headers.c_str(), sizeof(char), strlen(headers.c_str()), pGrad);
  fwrite(headers.c_str(), sizeof(char), strlen(headers.c_str()), pCount);
  for (size_t i = 0; i < sampleSize; ++i) {
    for (size_t j = 0; j < ndims; ++j) {
      pos[j] = table[j][i];
      fprintf(pGrad, " %.9f", table[j][i]);
      fprintf(pCount, " %.9f", table[j][i]);
    }
    fprintf(pCount, " %lu\n", getCount(pos, true));
    vector<double> f = getGradient(pos, true);
    for (size_t j = 0; j < ndims; ++j) {
      fprintf(pGrad, " %.9f", (f[j] / outputunit));
    }
    fprintf(pGrad, "\n");
  }
  fclose(pGrad);
  fclose(pCount);
  free(buffer1);
  free(buffer2);
  if (ndims == 1) {
    write1DPMF(filename);
  }
}

void DRRForceGrid::writeDivergence(const string &filename) const {
  string divname = filename + suffix + ".div";
  vector<double> pos(ndims, 0);
  FILE *pDiv;
  pDiv = fopen(divname.c_str(), "w");
  fwrite(headers.c_str(), sizeof(char), strlen(headers.c_str()), pDiv);
  for (size_t i = 0; i < sampleSize; ++i) {
    for (size_t j = 0; j < ndims; ++j) {
      pos[j] = table[j][i];
      fprintf(pDiv, " %.9f", table[j][i]);
    }
    const double divergence = getDivergence(pos);
    fprintf(pDiv, " %.9f", (divergence / outputunit));
    fprintf(pDiv, "\n");
  }
  fclose(pDiv);
}

bool ABF::store_getbias(const vector<double> &pos, const vector<double> &f,
                        vector<double> &fbias) {
  if (!isInBoundary(pos)) {
    std::fill(begin(fbias), end(fbias), 0.0);
    return false;
  }
  const size_t baseaddr = sampleAddress(pos);
  unsigned long int &count = samples[baseaddr];
  ++count;
  double factor = 2 * (static_cast<double>(count)) / mFullSamples - 1;
  auto it_fa = begin(forces) + baseaddr * ndims;
  auto it_fb = begin(fbias);
  auto it_f = begin(f);
  auto it_maxFactor = begin(mMaxFactors);
  do {
    // Clamp to [0,maxFactors]
    factor = factor < 0 ? 0 : factor > (*it_maxFactor) ? (*it_maxFactor) : factor;
    (*it_fa) += (*it_f); // Accumulate instantaneous force
    (*it_fb) = factor * (*it_fa) * (-1.0) /
               static_cast<double>(count); // Calculate bias force
    ++it_fa;
    ++it_fb;
    ++it_f;
    ++it_maxFactor;
  } while (it_f != end(f));

  return true;
}

ABF ABF::mergewindow(const ABF &aWA, const ABF &aWB) {
  const vector<DRRAxis> dR = merge(aWA.dimensions, aWB.dimensions);
  const string suffix = ".abf";
  ABF result(dR, suffix);
  const size_t nrows = result.sampleSize;
  const size_t ncols = result.ndims;
  vector<double> pos(ncols, 0);
  for (size_t i = 0; i < nrows; ++i) {
    for (size_t j = 0; j < ncols; ++j) {
      pos[j] = result.table[j][i];
    }
    const unsigned long int countA = aWA.getCount(pos);
    const unsigned long int countB = aWB.getCount(pos);
    const vector<double> aForceA = aWA.getAccumulatedForces(pos);
    const vector<double> aForceB = aWB.getAccumulatedForces(pos);
    result.store(pos, aForceA, countA);
    result.store(pos, aForceB, countB);
  }
  return result;
}

vector<double> CZAR::getGradient(const vector<double> &pos,
                                 bool SkipCheck) const {
  vector<double> result(ndims, 0);
  if (!SkipCheck) {
    if (!isInBoundary(pos)) {
      return result;
    }
  }
  if (kbt <= std::numeric_limits<double>::epsilon()) {
    std::cerr << "ERROR! The kbt shouldn't be zero when use CZAR estimator. "
              << '\n';
    std::abort();
  }
  const size_t baseaddr = sampleAddress(pos);
  const vector<double> log_deriv(getCountsLogDerivative(pos));
  const unsigned long int &count = samples[baseaddr];
  if (count == 0)
    return result;
  auto it_fa = begin(forces) + baseaddr * ndims;
  std::transform(it_fa, it_fa + ndims, begin(log_deriv), begin(result),
  [&count, this](double fa, double ld) {
    return fa * (-1.0) / count - kbt * ld;
  });
  return result;
}

CZAR CZAR::mergewindow(const CZAR &cWA, const CZAR &cWB) {
  const vector<DRRAxis> dR = merge(cWA.dimensions, cWB.dimensions);
  const double newkbt = cWA.kbt;
  const string suffix = ".czar";
  CZAR result(dR, suffix, newkbt);
  const size_t nrows = result.sampleSize;
  const size_t ncols = result.ndims;
  vector<double> pos(ncols, 0);
  for (size_t i = 0; i < nrows; ++i) {
    for (size_t j = 0; j < ncols; ++j) {
      pos[j] = result.table[j][i];
    }
    const unsigned long int countA = cWA.getCount(pos);
    const unsigned long int countB = cWB.getCount(pos);
    const vector<double> aForceA = cWA.getAccumulatedForces(pos);
    const vector<double> aForceB = cWB.getAccumulatedForces(pos);
    result.store(pos, aForceA, countA);
    result.store(pos, aForceB, countB);
  }
  return result;
}

void CZAR:: writeZCount(const string &filename) const {
  string countname = filename + ".zcount";
  vector<double> pos(ndims, 0);
  FILE *pCount;
  pCount = fopen(countname.c_str(), "w");
  char *buffer;
  buffer = (char *)malloc((sizeof(double)) * sampleSize * ndims);
  setvbuf(pCount, buffer, _IOFBF, (sizeof(double)) * sampleSize * ndims);
  fwrite(headers.c_str(), sizeof(char), strlen(headers.c_str()), pCount);
  for (size_t i = 0; i < sampleSize; ++i) {
    for (size_t j = 0; j < ndims; ++j) {
      pos[j] = table[j][i];
      fprintf(pCount, " %.9f", table[j][i]);
    }
    fprintf(pCount, " %lu\n", getCount(pos, true));
  }
  fclose(pCount);
  free(buffer);
}

}
}

#endif
