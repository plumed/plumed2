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
#include "MetaDGrid.h"
#include "MetaDHills.h"
#include "core/Action.h"
#include "core/Value.h"
#include "tools/File.h"
#include "tools/OpenMP.h"
#include "tools/Tools.h"
#include <algorithm>
#include <cmath>
#include <limits>

namespace PLMD {
namespace bias {

bool MetaDBiasGrid::setupBins(Action& action, const std::vector<Value*>& args, std::vector<std::string>& gmin, std::vector<std::string>& gmax,
                              std::vector<unsigned>& gbin, std::vector<double>& gspacing, bool adaptive, const std::vector<double>& sigma0,
                              const std::vector<double>& sigma0min) {
  bool nobounds=false;
  if(gmin.size()!=gmax.size()) {
    action.error("GRID_MAX and GRID_MIN should be either present or absent");
  }
  if(gbin.size()!=0     && gmin.size()==0) {
    action.error("If GRID_BIN is present also GRID_MIN and GRID_MAX should be present");
  }
  for(const auto s : gspacing) {
    if(!(s>0.0)) {
      action.error("GRID_SPACING should be positive");
    }
  }
  if(gmin.size()!=0) {
    if(gbin.size()==0 && gspacing.size()==0) {
      if(!adaptive) {
        action.log<<"  Binsize not specified, 1/5 of sigma will be be used\n";
        plumed_assert(sigma0.size()==args.size());
        gspacing.resize(args.size());
        for(unsigned i=0; i<gspacing.size(); i++) {
          gspacing[i]=0.2*sigma0[i];
        }
      } else {
        // with adaptive hills and grid a sigma min must be specified
        for(unsigned i=0; i<sigma0min.size(); i++)
          if(sigma0min[i]<=0) {
            action.error("When using ADAPTIVE Gaussians on a grid SIGMA_MIN must be specified");
          }
        action.log<<"  Binsize not specified, 1/5 of sigma_min will be be used\n";
        gspacing.resize(args.size());
        for(unsigned i=0; i<gspacing.size(); i++) {
          gspacing[i]=0.2*sigma0min[i];
        }
      }
    } else if(gspacing.size()!=0 && gbin.size()==0) {
      action.log<<"  The number of bins will be estimated from GRID_SPACING\n";
    } else if(gspacing.size()!=0 && gbin.size()!=0) {
      action.log<<"  You specified both GRID_BIN and GRID_SPACING\n";
      action.log<<"  The more conservative (highest) number of bins will be used for each variable\n";
    }
    if(gbin.size()==0) {
      gbin.assign(args.size(),1);
    }
    if(gspacing.size()!=0)
      for(unsigned i=0; i<args.size(); i++) {
        double a,b;
        Tools::convert(gmin[i],a);
        Tools::convert(gmax[i],b);
        unsigned n=std::ceil(((b-a)/gspacing[i]));
        if(gbin[i]<n) {
          gbin[i]=n;
        }
      }
  }
  if(gspacing.size()!=0 && gmin.size()==0) {
    // no range: the grid is placed around the first hill, on the lattice with a point at 0
    // (a first range of one bin is used until then); periodic variables are covered over the whole period
    nobounds=true;
    gmin.resize(args.size());
    gmax.resize(args.size());
    gbin.resize(args.size());
    for(unsigned i=0; i<args.size(); i++) {
      if(!(gspacing[i]>0.0)) {
        action.error("GRID_SPACING should be positive");
      }
      if(args[i]->isPeriodic()) {
        args[i]->getDomain(gmin[i],gmax[i]);
        double a,b;
        Tools::convert(gmin[i],a);
        Tools::convert(gmax[i],b);
        gbin[i]=std::ceil((b-a)/gspacing[i]);
      } else {
        gmin[i]="0";
        gmax[i]=exactString(gspacing[i]);
        gbin[i]=1;
      }
    }
  }
  return nobounds;
}

void MetaDBiasGrid::checkSpacing(Action& action, const std::vector<std::string>& gmin, const std::vector<std::string>& gmax,
                                 const std::vector<unsigned>& gbin, bool adaptive, const std::vector<double>& sigma0,
                                 const std::vector<double>& sigma0min) {
  for(unsigned i=0; i<gbin.size(); i++) {
    double a,b;
    Tools::convert(gmin[i],a);
    Tools::convert(gmax[i],b);
    const double mesh=(b-a)/static_cast<double>(gbin[i]);
    if(!adaptive) {
      if(mesh>0.5*sigma0[i]) {
        action.log<<"  WARNING: Using a "<<action.getName()<<" with a Grid Spacing larger than half of the Gaussians width (SIGMA) can produce artifacts\n";
      }
    } else {
      if(!(sigma0min[i]>0.0)) {
        action.error("When using ADAPTIVE Gaussians on a grid SIGMA_MIN must be specified");
      }
      if(mesh>0.5*sigma0min[i]) {
        action.log<<"  WARNING: to use a "<<action.getName()<<" with a GRID and ADAPTIVE you need to set a Grid Spacing lower than half of the Gaussians (SIGMA_MIN) \n";
      }
    }
  }
}

void MetaDBiasGrid::checkSizeLimit(Action& action, double limit) {
  if(!(limit>=1.0)) {
    action.error("GRID_SIZE_LIMIT should be at least 1");
  }
  // grid points are indexed with unsigned integers (see GridBase)
  if(limit>static_cast<double>(std::numeric_limits<unsigned>::max())) {
    action.error("GRID_SIZE_LIMIT cannot be larger than "+exactString(std::numeric_limits<unsigned>::max()));
  }
}

MetaDBiasGrid::MetaDBiasGrid(Action& action):
  action_(action) {
}

void MetaDBiasGrid::create(const std::string& funcl, const std::vector<Value*>& args, const std::vector<std::string>& gmin,
                           const std::vector<std::string>& gmax, const std::vector<unsigned>& gbin, bool sparse, bool spline, bool nobounds) {
  sparse_=sparse;
  if(!sparse) {
    grid_=Tools::make_unique<Grid>(funcl,args,gmin,gmax,gbin,spline,true);
  } else {
    grid_=Tools::make_unique<SparseGrid>(funcl,args,gmin,gmax,gbin,spline,true);
  }
  pristine_=nobounds;
  std::vector<std::string> actualmin=grid_->getMin();
  std::vector<std::string> actualmax=grid_->getMax();
  for(unsigned i=0; i<args.size(); i++) {
    std::string is;
    Tools::convert(i,is);
    if(gmin[i]!=actualmin[i]) {
      action_.error("GRID_MIN["+is+"] must be adjusted to "+actualmin[i]+" to fit periodicity");
    }
    if(gmax[i]!=actualmax[i]) {
      action_.error("GRID_MAX["+is+"] must be adjusted to "+actualmax[i]+" to fit periodicity");
    }
  }
}

void MetaDBiasGrid::read(IFile& file, const std::string& filename, const std::string& funcl, const std::vector<Value*>& args,
                         const std::vector<std::string>& gmin, const std::vector<std::string>& gmax, const std::vector<unsigned>& gbin,
                         bool sparse, bool spline, bool nobounds) {
  sparse_=sparse;
  grid_=GridBase::create(funcl, args, file, sparse, spline, true);
  if(grid_->getDimension()!=args.size()) {
    action_.error("mismatch between dimensionality of input grid and number of arguments");
  }
  for(unsigned i=0; i<args.size(); ++i) {
    if( args[i]->isPeriodic()!=grid_->getIsPeriodic()[i] ) {
      action_.error("periodicity mismatch between arguments and input bias");
    }
    // the grid in the file can have a different range (e.g. it was extended), but it must be on the same
    // lattice as the grid in input: same spacing, points at the same positions
    std::string is;
    Tools::convert(i,is);
    double min,max,fmin;
    Tools::convert(gmin[i],min);
    Tools::convert(gmax[i],max);
    Tools::convert(grid_->getMin()[i],fmin);
    const double dx=(max-min)/gbin[i];
    const double fdx=grid_->getDx(i);
    if(args[i]->isPeriodic()) {
      if(grid_->getNbin()[i]!=gbin[i]) {
        action_.error("the number of bins of the grid in "+filename+" for the periodic variable "+is+" does not match the input");
      }
    } else {
      const double shift=(fmin-min)/dx;
      // (the bounds of an extended grid are written with 12 significant digits)
      if(std::fabs(fdx-dx)>1e-6*dx || std::fabs(shift-std::round(shift))>1e-4) {
        action_.error("the grid in "+filename+" for variable "+is+" is not on the lattice of the grid in input: "
                      "its spacing should be "+exactString(dx)+" and its points should be at "+gmin[i]+" plus multiples of the spacing");
      }
    }
  }
  // a grid without initial range read from a file written before the first hill: it is placed around the first hill
  pristine_=false;
  if(nobounds) {
    bool empty=true;
    grid_->forEachStoredPoint([&](GridBase::index_t,double value) {
      empty=empty && value==0.0;
    });
    pristine_=empty;
  }
}

void MetaDBiasGrid::setSizeLimit(double limit) {
  sizelimit_=limit;
  // an initial grid larger than GRID_SIZE_LIMIT (e.g. given with GRID_MIN, GRID_MAX, GRID_BIN) can still be extended
  double size=1.0;
  for(const auto n : grid_->getNbin()) {
    size*=n;
  }
  if(size>sizelimit_) {
    sizelimit_=std::min(2.0*size,static_cast<double>(std::numeric_limits<unsigned>::max()));
    action_.log.printf("  The initial grid has more points than GRID_SIZE_LIMIT: the limit is set to %s points\n",exactString(sizelimit_).c_str());
  }
}

void MetaDBiasGrid::setInterval(double lower, double upper, double time) {
  interval_=true;
  lower_=lower;
  upper_=upper;
  // the bias beyond the limits is constant: the grid contains the interval and at least one point
  // beyond each limit, which holds this constant value (see extend)
  const std::vector<double> center(1,0.5*(lower_+upper_));
  const std::vector<double> halfwidth(1,0.5*(upper_-lower_)+grid_->getDx(0));
  extend(center,halfwidth,time);
}

void MetaDBiasGrid::extend(const std::vector<double>& x, const std::vector<double>& halfwidth, double time) {
  const unsigned ncv=grid_->getDimension();
  const std::vector<bool> pbc=grid_->getIsPeriodic();
  const std::vector<unsigned> nbin=grid_->getNbin();
  // points to add at the lower and upper end along each CV (for a pristine grid: the change of the range, also negative)
  std::vector<long long> lower(ncv,0), upper(ncv,0);
  // as above, with some more points so that the grid is not extended at every hill
  std::vector<long long> lowergrow(ncv,0), uppergrow(ncv,0);
  bool extend=false;
  for(unsigned i=0; i<ncv; i++) {
    if(pbc[i]) {
      continue;
    }
    // points within halfwidth on each side (as in MetaDKernel::getSupport, which adds the hill to the points from
    // the one below x minus this number to the one below x plus this number), plus one more on each side,
    // so that the bias and its derivatives are zero at the boundary of the grid
    double n=std::ceil(halfwidth[i]/grid_->getDx(i));
    if(!std::isfinite(n) || n<0.0) {
      // support not finite (degenerate hill): the hill is added to the whole grid, x only must be in it
      n=0.0;
    }
    const double k=std::floor(grid_->getLatticeCoordinate(i,x[i]));
    const double lo=k-n-1.0;
    const double hi=k+n+1.0;
    if(!std::isfinite(lo) || !std::isfinite(hi) || std::fabs(lo)>1e15 || std::fabs(hi)>1e15) {
      std::string is;
      Tools::convert(i,is);
      action_.error("cannot extend the grid along variable "+is+" to contain the point "+exactString(x[i]));
    }
    const double last=static_cast<double>(nbin[i])-1.0;
    if(pristine_) {
      lower[i]=lowergrow[i]=static_cast<long long>(-lo);
      upper[i]=uppergrow[i]=static_cast<long long>(hi-last);
      extend=true;
    } else {
      // geometric growth: at least a quarter of the current number of points
      const long long grow=static_cast<long long>(std::ceil(0.25*nbin[i]));
      if(lo<0.0) {
        lower[i]=static_cast<long long>(-lo);
        lowergrow[i]=std::max(lower[i],grow);
        extend=true;
      }
      if(hi>last) {
        upper[i]=static_cast<long long>(hi-last);
        uppergrow[i]=std::max(upper[i],grow);
        extend=true;
      }
    }
  }
  const bool pristine=pristine_;
  pristine_=false;
  if(!extend) {
    return;
  }
  auto size=[&](const std::vector<long long>& lo,const std::vector<long long>& up) {
    double n=1.0;
    for(unsigned i=0; i<ncv; i++) {
      n*=static_cast<double>(nbin[i])+lo[i]+up[i];
    }
    return n;
  };
  if(size(lowergrow,uppergrow)>sizelimit_) {
    lowergrow=lower;
    uppergrow=upper;
    if(size(lower,upper)>sizelimit_) {
      std::string point;
      for(unsigned i=0; i<ncv; i++) {
        point+=(i>0?",":"")+exactString(x[i]);
      }
      action_.error("a Gaussian at "+point+" requires a grid with "+exactString(size(lower,upper))+" points, more than GRID_SIZE_LIMIT="
                    +exactString(sizelimit_)+": increase GRID_SIZE_LIMIT, or check the variables (e.g. a distance that grows without limit) "
                    "and the width of the Gaussians (e.g. SIGMA_MAX with ADAPTIVE)");
    }
  }
  grid_->resize(lowergrow,uppergrow);
  if(interval_ && !pristine) {
    // With INTERVAL, all the points beyond a limit have the same value (the hills are computed at the limit there) and
    // zero derivative: the new points beyond a limit take the value of the old point at that end of the grid.
    const unsigned n=grid_->getNbin()[0];
    std::vector<double> der(1);
    auto fill=[&](unsigned from,unsigned first,unsigned last,bool beyond) {
      if(!beyond) {
        return;
      }
      const double value=grid_->getValueAndDerivatives(from,der);
      der[0]=0.0;
      for(unsigned k=first; k<last; k++) {
        grid_->setValueAndDerivatives(k,value,der);
      }
    };
    const unsigned first=static_cast<unsigned>(lowergrow[0]);
    const unsigned last=n-1-static_cast<unsigned>(uppergrow[0]);
    // (a point at the limit, e.g. in a grid file written before grids were extended, has the same value)
    fill(first,0,first,lowergrow[0]>0 && grid_->getPoint(first)[0]<=lower_);
    fill(last,last+1,n,uppergrow[0]>0 && grid_->getPoint(last)[0]>=upper_);
  }
  const std::vector<std::string> gmin=grid_->getMin();
  const std::vector<std::string> gmax=grid_->getMax();
  const std::vector<unsigned> newnbin=grid_->getNbin();
  Log& log=action_.log;
  if(midline_) {
    log.printf("\n");
    midline_=false;
  }
  log.printf("  %s grid: min",pristine ? "Placing the" : "Extending the");
  for(unsigned i=0; i<ncv; i++) {
    log.printf(" %s",gmin[i].c_str());
  }
  log.printf(" max");
  for(unsigned i=0; i<ncv; i++) {
    log.printf(" %s",gmax[i].c_str());
  }
  log.printf(" points");
  for(unsigned i=0; i<ncv; i++) {
    log.printf(" %u",newnbin[i]);
  }
  log.printf(" at time %f\n",time);
}

void MetaDBiasGrid::addHill(const MetaDGaussian& hill, const MetaDKernel& kernel, double time) {
  // the grid must contain all the points where the hill is not zero
  extend(hill.center,kernel.getCutoff(hill),time);
  if(!interval_) {
    addHillToPoints(hill,kernel);
  } else {
    size_t ncv=grid_->getDimension();
    std::vector<unsigned> nneighb=kernel.getSupport(hill,*grid_);
    std::vector<Grid::index_t> neighbors=grid_->getNeighbors(hill.center,nneighb);
    std::vector<double> der(ncv);
    std::vector<double> xx(ncv);
    // every MPI process adds the whole hill to its copy of the grid (see addHillToPoints)
    std::vector<double> dp(ncv);
    for(size_t i=0; i<neighbors.size(); ++i) {
      Grid::index_t ineigh=neighbors[i];
      for(unsigned j=0; j<ncv; ++j) {
        der[j]=0.0;
      }
      grid_->getPoint(ineigh,xx);
      double bias=kernel.evaluateWithDerivatives(xx,hill,der.data(),dp.data());
      grid_->addValueAndDerivatives(ineigh,bias,der);
    }
  }
}

void MetaDBiasGrid::addHillToPoints(const MetaDGaussian& hill, const MetaDKernel& kernel) {
  // The hill is evaluated on the grid points within its support. The difference of a grid point from the
  // center of the hill along CV d only depends on the index of the point along d, so the differences are
  // computed once per direction; the value at each point is then computed from them
  // (MetaDKernel::evaluateFromDifferences, same arithmetic as MetaDKernel::evaluateWithDerivatives).
  const unsigned ncv=grid_->getDimension();
  const std::vector<unsigned> nneigh=kernel.getSupport(hill,*grid_);
  const std::vector<unsigned> nbin=grid_->getNbin();
  const std::vector<bool> pbc=grid_->getIsPeriodic();
  std::vector<unsigned> center(ncv);
  grid_->getIndices(hill.center,center);
  // grid indices and differences along each direction (as in GridBase::getNeighbors)
  std::vector<std::vector<unsigned>> indices(ncv);
  std::vector<std::vector<double>> differences(ncv);
  std::vector<unsigned> tmpindices(center);
  std::vector<double> point(ncv);
  for(unsigned d=0; d<ncv; d++) {
    unsigned nb=2*nneigh[d]+1;
    if(pbc[d] && nb>nbin[d]) {
      nb=nbin[d];
    }
    for(unsigned k=0; k<nb; k++) {
      int i0=static_cast<int>(k)-static_cast<int>(nneigh[d])+static_cast<int>(center[d]);
      const int n=nbin[d];
      if(pbc[d]) {
        i0=((i0%n)+n)%n;
      } else if(i0<0 || i0>=n) {
        continue;
      }
      tmpindices[d]=i0;
      grid_->getPoint(tmpindices.data(),ncv,point.data(),ncv);
      double dd=kernel.difference(d,hill.center[d],point[d]);
      if(!hill.multivariate) {
        dd*=hill.invsigma[d];
      }
      indices[d].push_back(i0);
      differences[d].push_back(dd);
    }
    tmpindices[d]=center[d];
  }
  // linear index of a point (column-major order, as in GridBase::getIndex)
  std::vector<Grid::index_t> stride(ncv);
  Grid::index_t npoints=1;
  for(unsigned d=0; d<ncv; d++) {
    stride[d]=(d==0 ? 1 : stride[d-1]*nbin[d-1]);
    npoints*=indices[d].size();
  }
  if(npoints==0) {
    return;
  }
  // Loop over all the points of the box, the first direction being the fastest.
  // With MPI every process adds the whole hill to its copy of the grid: the hills are the same on all
  // processes, so the grids stay identical without communication (splitting the points among the processes
  // requires summing values and derivatives of the whole box at every hill, which costs more than computing them).
  // With OpenMP the box is split in blocks of consecutive points; each point is computed and added by one
  // thread only, with the same arithmetic as the serial loop, so the result does not depend on the number of
  // threads. Sparse grids are filled serially, as new points are inserted.
  unsigned nt=OpenMP::getNumThreads();
  if(sparse_ || npoints<static_cast<Grid::index_t>(1024*nt)) {
    nt=1;
  }
  ThreadWorkArrays<unsigned> kwork(nt,ncv);
  ThreadWorkArrays<double> dpwork(nt,ncv);
  ThreadWorkArrays<double> derwork(nt,ncv);
  GridBase* grid=grid_.get();
  const bool sparse=sparse_;
  #pragma omp parallel for num_threads(nt) schedule(static,1)
  for(unsigned block=0; block<nt; block++) {
    const Grid::index_t first=(npoints*block)/nt;
    const Grid::index_t last=(npoints*(block+1))/nt;
    // indices of the first point of the block along each direction
    unsigned* k=kwork.get(block);
    Grid::index_t m=first;
    for(unsigned d=0; d<ncv; d++) {
      k[d]=m%indices[d].size();
      m/=indices[d].size();
    }
    double* dp=dpwork.get(block);
    double* der=derwork.get(block);
    for(Grid::index_t n=first; n<last; n++) {
      for(unsigned d=0; d<ncv; d++) {
        dp[d]=differences[d][k[d]];
        der[d]=0.0;
      }
      bool inside;
      const double bias=kernel.evaluateFromDifferences(dp,hill,der,inside);
      if(inside || sparse) {
        // adding zero outside the cutoff only matters for sparse grids, where it creates the point
        Grid::index_t index=0;
        for(unsigned d=0; d<ncv; d++) {
          index+=indices[d][k[d]]*stride[d];
        }
        grid->addValueAndDerivatives(index,bias,der,ncv);
      }
      // next point
      for(unsigned d=0; d<ncv; d++) {
        if(++k[d]<indices[d].size()) {
          break;
        }
        k[d]=0;
      }
    }
  }
}

double MetaDBiasGrid::getBias(const std::vector<double>& cv) const {
  double bias=0.0;
  // the grid contains all the points where the hills are not zero
  if(interval_ && isIntervalBeyondGrid(cv[0])) {
    bias = getIntervalBiasOutsideGrid(cv[0]);
  } else if(grid_->isInside(cv)) {
    bias = grid_->getValue(cv);
  }
  return bias;
}

double MetaDBiasGrid::getBiasAndDerivatives(const std::vector<double>& cv, std::vector<double>& der) const {
  const unsigned ncv=grid_->getDimension();
  double bias=0.0;
  // the grid contains all the points where the hills are not zero
  if(grid_->isInside(cv) && !(interval_ && isIntervalBeyondGrid(cv[0]))) {
    std::vector<double> vder(ncv);
    bias=grid_->getValueAndDerivatives(cv,vder);
    for(unsigned i=0; i<ncv; i++) {
      der[i]=vder[i];
    }
  } else {
    if(interval_) {
      bias=getIntervalBiasOutsideGrid(cv[0]);
    }
    for(unsigned i=0; i<ncv; i++) {
      der[i]=0.0;
    }
  }
  return bias;
}

bool MetaDBiasGrid::isIntervalBeyondGrid(double cv) const {
  // the grid contains the interval and at least one point beyond each limit, and the bias is constant beyond the limits:
  // beyond the last point of the grid (where the spline would use this point only) the bias is the value at this point
  const double k=grid_->getLatticeCoordinate(0,cv);
  return !(k>=0.0 && k<static_cast<double>(grid_->getNbin()[0])-1.0);
}

double MetaDBiasGrid::getIntervalBiasOutsideGrid(double cv) const {
  const unsigned n=grid_->getNbin()[0];
  const bool below=!(grid_->getLatticeCoordinate(0,cv)>=0.0);
  return grid_->getValue(below ? 0 : n-1);
}

}
}
