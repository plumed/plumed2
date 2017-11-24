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

bool DRRAxis::isInBoundary(double x) const {
  if (x < min || x > max)
    return false;
  else
    return true;
}

bool DRRAxis::isRealPeriodic() const {
  if (periodic == true) {
    if(std::abs(domainMax - max) < binWidth && std::abs(domainMin - min) < binWidth) {
      return true;
    }
    else {
      return false;
    }
  }
  else {
    return false;
  }
}

DRRAxis DRRAxis::merge(const DRRAxis &d1, const DRRAxis &d2) {
  const double newmin = std::min(d1.getMin(), d2.getMin());
  const double newmax = std::max(d1.getMax(), d2.getMax());
  const double newWidth = d1.getWidth();
  const size_t newbins = size_t(std::nearbyint((newmax - newmin) / newWidth));
  const bool newpbc = d1.isPeriodic();
  const double newdmin = std::min(d1.getDomainMin(), d2.getDomainMin());
  const double newdmax = std::max(d1.getDomainMax(), d2.getDomainMax());
  DRRAxis result(newmin, newmax, newbins, newpbc, newdmin, newdmax);
  return result;
}

std::vector<double> DRRAxis::getMiddlePoints() {
  std::vector<double> result(nbins, 0);
  const double width = getWidth();
  double temp = min - width / 2;
  std::generate(std::begin(result), std::end(result), [&]() {
    temp += width;
    return temp;
  });
  return result;
}

size_t DRRForceGrid::index1D(const DRRAxis &c, double x) {
#ifdef DEBUG_DRR
  if (x < c.min || x > c.max) {
    std::cerr << "This is a bug!" << '\n';
    std::cerr << "CV should be larger than minimal value or smaller than the "
              "maximum value of dimension."
              << '\n';
    std::cerr << "min = " << c.min << '\n';
    std::cerr << "max = " << c.max << '\n';
    std::cerr << "x = " << x << std::endl;
    std::abort();
  }
#endif
  size_t idx = size_t(std::floor((x - c.min) / c.binWidth));
  idx = (idx == c.nbins) ? (c.nbins - 1) : idx;
  return idx;
}

void DRRForceGrid::fillTable(const std::vector<std::vector<double>> &in) {
  table.resize(ndims, std::vector<double>(sampleSize, 0));
  for (size_t i = 0; i < ndims; ++i) {
    size_t repeatAll = 1, repeatOne = 1;
    for (size_t j = i + 1; j < ndims; ++j)
      repeatOne *= in[j].size();
    for (size_t j = 0; j < i; ++j)
      repeatAll *= in[j].size();
    size_t in_i_sz = in[i].size();
    for (size_t l = 0; l < in_i_sz; ++l)
      std::fill_n(std::begin(table[i]) + l * repeatOne, repeatOne, in[i][l]);
    for (size_t k = 0; k < repeatAll - 1; ++k)
      std::copy_n(std::begin(table[i]), repeatOne * in_i_sz,
                  std::begin(table[i]) + repeatOne * in_i_sz * (k + 1));
  }
}

DRRForceGrid::DRRForceGrid()
  : suffix(""), ndims(0), dimensions(0), sampleSize(0), forceSize(0),
    headers(""), table(0), forces(0), samples(0), endpoints(0), shifts(0), outputunit(1.0) {}

DRRForceGrid::DRRForceGrid(const std::vector<DRRAxis> &p_dimensions,
                           const std::string &p_suffix, bool initializeTable)
  : suffix(p_suffix), ndims(p_dimensions.size()), dimensions(p_dimensions) {
  sampleSize = 1;
  std::vector<std::vector<double>> mp(ndims);
  std::stringstream ss;
  ss << "# " << ndims << '\n';
  shifts.resize(ndims, 0);
  for (size_t i = 0; i < ndims; ++i) {
    sampleSize = dimensions[i].nbins * sampleSize;
    mp[i] = dimensions[i].getMiddlePoints();
    shifts[i] = std::accumulate(
                  std::begin(dimensions), std::begin(dimensions) + i, size_t(1),
    [](size_t k, const DRRAxis &d) { return k * d.getBins(); });
    ss.precision(std::numeric_limits<double>::max_digits10);
    ss << std::fixed << "# " << dimensions[i].min << ' '
       << dimensions[i].getWidth() << ' ' << dimensions[i].nbins;
    if (dimensions[i].isPeriodic())
      ss << " 1" << '\n';
    else
      ss << " 0" << '\n';
  }
  headers = ss.str();
  if (initializeTable)
    fillTable(mp);
  forceSize = sampleSize * ndims;
  forces.resize(forceSize, 0.0);
  samples.resize(sampleSize, 0);
  outputunit = 1.0;
  // For 1D pmf
  if (ndims == 1) {
    endpoints.resize(dimensions[0].nbins + 1, 0);
    double ep = dimensions[0].min;
    double stride = dimensions[0].getWidth();
    for (auto &i : endpoints) {
      i = ep;
      ep += stride;
    }
  }
}

bool DRRForceGrid::isInBoundary(const std::vector<double> &pos) const {
  bool result = true;
  for (size_t i = 0; i < ndims; ++i) {
    if (pos[i] < dimensions[i].min || pos[i] > dimensions[i].max)
      return false;
  }
  return result;
}

size_t DRRForceGrid::sampleAddress(const std::vector<double> &pos) const {
  size_t saddr = 0;
  for (size_t i = 0; i < ndims; ++i) {
    saddr += shifts[i] * index1D(dimensions[i], pos[i]);
  }
  return saddr;
}

bool DRRForceGrid::store(const std::vector<double> &pos,
                         const std::vector<double> &f,
                         unsigned long int nsamples) {
  if (isInBoundary(pos)) {
    if (nsamples == 0)
      return true;
    const size_t baseaddr = sampleAddress(pos) * ndims;
    samples[baseaddr / ndims] += nsamples;
    auto it_fa = std::begin(forces) + baseaddr;
    std::transform(std::begin(f), std::end(f), it_fa, it_fa,
                   std::plus<double>());
    return true;
  } else {
    return false;
  }
}

std::vector<DRRAxis> DRRForceGrid::merge(const std::vector<DRRAxis> &dA,
    const std::vector<DRRAxis> &dB) {
  std::vector<DRRAxis> dR(dA.size());
  std::transform(std::begin(dA), std::end(dA), std::begin(dB), std::begin(dR),
                 DRRAxis::merge);
  return dR;
}

std::vector<double>
DRRForceGrid::getAccumulatedForces(const std::vector<double> &pos) const {
  std::vector<double> result(ndims, 0);
  if (!isInBoundary(pos))
    return result;
  const size_t baseaddr = sampleAddress(pos) * ndims;
  std::copy(std::begin(forces) + baseaddr,
            std::begin(forces) + baseaddr + ndims, std::begin(result));
  return result;
}

unsigned long int DRRForceGrid::getCount(const std::vector<double> &pos,
    bool SkipCheck) const {
  if (!SkipCheck) {
    if (!isInBoundary(pos)) {
      return 0;
    }
  }
  return samples[sampleAddress(pos)];
}

std::vector<double> DRRForceGrid::getGradient(const std::vector<double> &pos,
    bool SkipCheck) const {
  std::vector<double> result(ndims, 0);
  if (!SkipCheck) {
    if (!isInBoundary(pos)) {
      return result;
    }
  }
  const size_t baseaddr = sampleAddress(pos) * ndims;
  if (samples[baseaddr / ndims] == 0)
    return result;
  auto it_fa = std::begin(forces) + baseaddr;
  std::transform(it_fa, it_fa + ndims, std::begin(result), [&](double fa) {
    return (-1.0) * fa / samples[baseaddr / ndims];
  });
  return result;
}

std::vector<double>
DRRForceGrid::getCountsLogDerivative(const std::vector<double> &pos) const {
  const size_t addr = sampleAddress(pos);
  const unsigned long int count_this = samples[addr];
  std::vector<double> result(ndims, 0);
  for (size_t i = 0; i < ndims; ++i) {
    const double binWidth = dimensions[i].getWidth();
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

// Write the gradients to a .grad file.
// void DRRForceGrid::writeGrad(std::string filename) const {
//   std::stringstream ssg;
//   std::fstream fsgrad;
//   filename += suffix + ".grad";
//   ssg << headers;
//   ssg << std::left << std::fixed << std::setprecision(OUTPUTPRECISION);
//   for (size_t i = 0; i < sampleSize; ++i) {
//     std::vector<double> pos(ndims, 0);
//     for (size_t j = 0; j < ndims; ++j) {
//       pos[j] = table[j][i];
//       ssg << ' ' << table[j][i];
//     }
//     const std::vector<double> f = getGradient(pos, true);
//     for (const auto &i : f)
//       ssg << ' ' << i;
//     ssg << '\n';
//   }
//   fsgrad.open(filename.c_str(), std::ios_base::out);
//   fsgrad.write(ssg.str().c_str(), ssg.str().length());
//   fsgrad.close();
// }

void DRRForceGrid::write1DPMF(std::string filename) const {
  filename += suffix + ".pmf";
  FILE *ppmf;
  ppmf = fopen(filename.c_str(), "w");
  const double w = dimensions[0].getWidth();
  double pmf = 0;
  fprintf(ppmf, "%.9f %.9f\n", endpoints[0], pmf);
  for (size_t i = 0; i < dimensions[0].nbins; ++i) {
    std::vector<double> pos(1, 0);
    pos[0] = table[0][i];
    const std::vector<double> f = getGradient(pos, true);
    pmf += f[0] * w;
    fprintf(ppmf, "%.9f %.9f\n", endpoints[i + 1], pmf);
  }
  fclose(ppmf);
}

// Write the gradients to a .count file.
// void DRRForceGrid::writeCount(std::string filename) const {
//   std::stringstream ssc;
//   std::fstream fscount;
//   filename += suffix + ".count";
//   ssc << headers;
//   ssc << std::left << std::fixed << std::setprecision(OUTPUTPRECISION);
//   std::vector<double> pos(ndims, 0);
//   for (size_t i = 0; i < sampleSize; ++i) {
//     for (size_t j = 0; j < ndims; ++j) {
//       pos[j] = table[j][i];
//       ssc << ' ' << table[j][i];
//     }
//     ssc << ' ' << getCount(pos, true) << '\n';
//   }
//   fscount.open(filename.c_str(), std::ios_base::out);
//   fscount.write(ssc.str().c_str(), ssc.str().length());
//   fscount.close();
// }

void DRRForceGrid::writeAll(const std::string &filename) const {
  std::string countname = filename + suffix + ".count";
  std::string gradname = filename + suffix + ".grad";
  std::vector<double> pos(ndims, 0);
  FILE *pGrad, *pCount;
  pGrad = fopen(gradname.c_str(), "w");
  pCount = fopen(countname.c_str(), "w");
  char *buffer1, *buffer2;
  buffer1 = (char*)malloc((sizeof(double))*sampleSize*ndims);
  buffer2 = (char*)malloc((sizeof(double))*sampleSize*ndims);
  setvbuf(pGrad, buffer1, _IOFBF, (sizeof(double))*sampleSize*ndims);
  setvbuf(pCount, buffer2, _IOFBF, (sizeof(double))*sampleSize*ndims);
  fwrite(headers.c_str(), sizeof(char), strlen(headers.c_str()), pGrad);
  fwrite(headers.c_str(), sizeof(char), strlen(headers.c_str()), pCount);
  for (size_t i = 0; i < sampleSize; ++i) {
    for (size_t j = 0; j < ndims; ++j) {
      pos[j] = table[j][i];
      fprintf(pGrad, " %.9f", table[j][i]);
      fprintf(pCount, " %.9f", table[j][i]);
    }
    fprintf(pCount, " %lu\n", getCount(pos, true));
    std::vector<double> f = getGradient(pos, true);
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

bool ABF::store_getbias(const std::vector<double> &pos,
                        const std::vector<double> &f,
                        std::vector<double> &fbias, double fullsamples) {
  if (!isInBoundary(pos)) {
    std::fill(std::begin(fbias), std::end(fbias), 0.0);
    return false;
  }
  const size_t baseaddr = sampleAddress(pos);
  unsigned long int &count = samples[baseaddr];
  ++count;
  double factor = 2 * (static_cast<double>(count)) / fullsamples - 1;
  factor = factor < 0 ? 0 : factor > 1 ? 1 : factor; // Clamp to [0,1]
  auto it_fa = std::begin(forces) + baseaddr * ndims;
  auto it_fb = std::begin(fbias);
  auto it_f = std::begin(f);
  do {
    (*it_fa) += (*it_f); // Accumulate instantaneous force
    (*it_fb) =
      factor * (*it_fa) * (-1.0) / static_cast<double>(count); // Calculate bias force
    ++it_fa;
    ++it_fb;
    ++it_f;
  } while (it_f != std::end(f));

  return true;
}

ABF ABF::mergewindow(const ABF &aWA, const ABF &aWB) {
  const std::vector<DRRAxis> dA = aWA.getDimensions();
  const std::vector<DRRAxis> dB = aWB.getDimensions();
  const std::vector<DRRAxis> dR = merge(dA, dB);
  const std::string suffix = ".abf";
  ABF result(dR, suffix);
  const std::vector<std::vector<double>> table = result.getTable();
  const size_t nrows = result.getSampleSize();
  const size_t ncols = result.getNumberOfDimension();
  std::vector<double> pos(ncols, 0);
  for (size_t i = 0; i < nrows; ++i) {
    for (size_t j = 0; j < ncols; ++j) {
      pos[j] = table[j][i];
    }
    const unsigned long int countA = aWA.getCount(pos, false);
    const unsigned long int countB = aWB.getCount(pos, false);
    const std::vector<double> aForceA = aWA.getAccumulatedForces(pos);
    const std::vector<double> aForceB = aWB.getAccumulatedForces(pos);
    result.store(pos, aForceA, countA);
    result.store(pos, aForceB, countB);
  }
  return result;
}

std::vector<double> CZAR::getGradient(const std::vector<double> &pos,
                                      bool SkipCheck) const {
  std::vector<double> result(ndims, 0);
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
  const size_t baseaddr = sampleAddress(pos) * ndims;
  const std::vector<double> log_deriv(getCountsLogDerivative(pos));
  if (samples[baseaddr / ndims] == 0)
    return result;
  auto it_fa = std::begin(forces) + baseaddr;
  std::transform(it_fa, it_fa + ndims, std::begin(log_deriv),
  std::begin(result), [&](double fa, double ld) {
    return fa * (-1.0) / samples[baseaddr / ndims] - kbt * ld;
  });
  return result;
}

CZAR CZAR::mergewindow(const CZAR &cWA, const CZAR &cWB) {
  const std::vector<DRRAxis> dA = cWA.getDimensions();
  const std::vector<DRRAxis> dB = cWB.getDimensions();
  const std::vector<DRRAxis> dR = merge(dA, dB);
  const double newkbt = cWA.getkbt();
  const std::string suffix = ".czar";
  CZAR result(dR, suffix, newkbt);
  const std::vector<std::vector<double>> table = result.getTable();
  const size_t nrows = result.getSampleSize();
  const size_t ncols = result.getNumberOfDimension();
  std::vector<double> pos(ncols, 0);
  for (size_t i = 0; i < nrows; ++i) {
    for (size_t j = 0; j < ncols; ++j) {
      pos[j] = table[j][i];
    }
    const unsigned long int countA = cWA.getCount(pos);
    const unsigned long int countB = cWB.getCount(pos);
    const std::vector<double> aForceA = cWA.getAccumulatedForces(pos);
    const std::vector<double> aForceB = cWB.getAccumulatedForces(pos);
    result.store(pos, aForceA, countA);
    result.store(pos, aForceB, countB);
  }
  return result;
}

}
}

#endif
