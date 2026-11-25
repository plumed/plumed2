/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2023 The plumed team
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
#ifndef __PLUMED_tools_Grid_h
#define __PLUMED_tools_Grid_h

#include <vector>
#include <array>
#include <string>
#include <map>
#include <cmath>
#include <memory>
#include <cstddef>

#include "Exception.h"

namespace PLMD {


// simple function to enable various weighting

class WeightBase {
public:
  virtual double projectInnerLoop(double &input, double &v)=0;
  virtual double projectOuterLoop(double &v)=0;
  virtual ~WeightBase() {}
};

class BiasWeight:public WeightBase {
public:
  double beta,invbeta;
  double shift=0.0;
  explicit BiasWeight(double v) {
    beta=v;
    invbeta=1./beta;
  }
  //double projectInnerLoop(double &input, double &v) override {return  input+exp(beta*v);}
  //double projectOuterLoop(double &v) override {return -invbeta*std::log(v);}
  double projectInnerLoop(double &input, double &v) override {
    auto betav=beta*v;
    auto x=betav-shift;
    if(x>0) {
      shift=betav;
      return input*std::exp(-x)+1;
    } else {
      return input+std::exp(x);
    }
  }
  double projectOuterLoop(double &v) override {
    auto res=-invbeta*(std::log(v)+shift);
    shift=0.0;
    return res;
  }
};

class ProbWeight:public WeightBase {
public:
  double beta,invbeta;
  explicit ProbWeight(double v) {
    beta=v;
    invbeta=1./beta;
  }
  double projectInnerLoop(double &input, double &v) override {
    return  input+v;
  }
  double projectOuterLoop(double &v) override {
    return -invbeta*std::log(v);
  }
};






class Value;
class IFile;
class OFile;
class KernelFunctions;
class Communicator;

/// \ingroup TOOLBOX
class GridBase {
public:
// we use a size_t here
// should be 8 bytes on all 64-bit machines
// and more portable than "unsigned long long"
  typedef std::size_t index_t;
// to restore old implementation (unsigned) use the following instead:
// typedef unsigned index_t;
/// Maximum dimension (exaggerated value).
/// Can be used to replace local std::vectors with std::arrays (allocated on stack).
  static constexpr std::size_t maxdim=16;


  /**
  Auxiliary class for accelerating grid calculations.

  Many functions in the Grid class have loops running over the dimensions.
  These loops are always with a little number of iterations, so would likely
  benefit unrolling. Even better, if the compiler knew the dimension a priori
  it could remove the loops in most common cases (e.g. dimension=1 or 2).

  To obtain this, I added this AcceleratorBase class. Specifically:
  - Any method in the Grid class that is a bottleneck in performance and
    has a loop can be moved to this auxiliary class.
  - The method should be declared here as a pure virtual function
    of class AcceleratorBase. The implementation instead goes in
    class Accelerator
  - Class Accelerator is derived from AcceleratorBase but depends on
    dimension as a template variable.
  - The function AcceleratorBase::create takes case of constructing
    accelerators for each possible value of dimension.

  The trick works because we don't expect dimension to ever be very
  high. By restricting it to <=16, we can have a limited number of
  compiled instances. The compiler will choose at runtime which one should
  be used.

  This class, as well as the derived classes, are expected to be
  empty. In other words, they only contain the pointer to a virtual
  table that the compiler will use to call the proper implementation
  variant, with the correct dimension.

  \warning
  Interface might change at any time.
  Do not use this outside of class GridBase and children.
  */
  class AcceleratorBase {
  public:
    /// Creates an accelerator with proper dimension
    static std::unique_ptr<AcceleratorBase> create(unsigned dim);
    virtual ~AcceleratorBase() = default;
    /// Can be used to check which value of dimension was used
    virtual unsigned getDimension() const=0;
    virtual std::vector<GridBase::index_t> getNeighbors(const GridBase& grid, const std::vector<unsigned> & nbin_,const std::vector<bool> & pbc_,const unsigned* indices,std::size_t indices_size, const std::vector<unsigned> &nneigh) const=0;
    virtual GridBase::index_t getIndex(const GridBase& grid, const std::vector<unsigned> & nbin_, const unsigned* indices,std::size_t indices_size) const=0;
    virtual void getPoint(const std::vector<double> & min_,const std::vector<double> & dx_, const unsigned* indices,std::size_t indices_size,double* point,std::size_t point_size) const=0;
    virtual void getIndices(const std::vector<unsigned> & nbin_, GridBase::index_t index, unsigned* indices, std::size_t indices_size) const=0;
    virtual void getIndices(const std::vector<double> & min_,const std::vector<double> & dx_, const std::vector<double> & x, unsigned* rindex_data,std::size_t rindex_size) const=0;
  };

  /**
  Auxiliary class for managing AcceleratorBase.

  class GridBase should contain a std::unique_ptr<AcceleratorBase>, which
  can be used to access the specialized versions. However, this would
  make a GridBase object not copyable. Instead of redefining
  copy constructor and copy assignment for GridBase, which has a lot of
  members, we use this wrapper class that just manages the lifetime
  of the underlying Accelerator object.

  The underlying object is made accessible through the -> operator,
  so that one can simply call functions as accelerator->function.

  \warning
  Interface might change at any time.
  Do not use this outside of class GridBase and children.
  */
  class AcceleratorHandler {
    /// This is the underlying pointer.
    std::unique_ptr<AcceleratorBase> ptr;
  public:
    /// Enable access to methods of AcceleratorBase
    const AcceleratorBase* operator->() const {
      return ptr.get();
    }
    /// Enable access to methods of AcceleratorBase
    AcceleratorBase* operator->() {
      return ptr.get();
    }
    /// Conversion to bool allows to check if the ptr has been set
    explicit operator bool() const noexcept {
      return bool(ptr);
    }
    /// Move constructor
    AcceleratorHandler(const AcceleratorHandler& other):
      ptr((other.ptr?AcceleratorBase::create(other->getDimension()):nullptr))
    {}
    /// Move assignment
    AcceleratorHandler & operator=(const AcceleratorHandler & other) {
      if(this!=&other) {
        ptr.reset();
        if(other.ptr) {
          ptr=AcceleratorBase::create(other->getDimension());
        }
      }
      return *this;
    }
    /// Constructor without arguments result in a non-usable accelerator (dimension is unspecified)
    AcceleratorHandler() = default;
    /// Constructor with an argument creates an accelerator with a fixed dimensionality
    AcceleratorHandler(unsigned dimension):
      ptr(AcceleratorBase::create(dimension))
    {}
  };

protected:
  AcceleratorHandler accelerator;
  std::string funcname;
  std::vector<std::string> argnames;
  std::vector<std::string> str_min_, str_max_;
  std::vector<double> min_,max_,dx_;
  std::vector<unsigned> nbin_;
  std::vector<bool> pbc_;
  index_t maxsize_;
  unsigned dimension_;
  bool dospline_, usederiv_;
  std::string fmt_; // format for output
/// get "neighbors" for spline
  unsigned getSplineNeighbors(const unsigned* indices, std::size_t indices_size, index_t* neighbors, std::size_t neighbors_size)const;
// std::vector<index_t> getSplineNeighbors(const std::vector<unsigned> & indices)const;


public:
/// this constructor here is Value-aware
  GridBase(const std::string& funcl, const std::vector<Value*> & args, const std::vector<std::string> & gmin,
           const std::vector<std::string> & gmax, const std::vector<unsigned> & nbin, bool dospline,
           bool usederiv);
/// this constructor here is not Value-aware
  GridBase(const std::string& funcl, const std::vector<std::string> &names, const std::vector<std::string> & gmin,
           const std::vector<std::string> & gmax, const std::vector<unsigned> & nbin, bool dospline,
           bool usederiv, const std::vector<bool> &isperiodic, const std::vector<std::string> &pmin,
           const std::vector<std::string> &pmax );
/// this is the real initializator
  void Init(const std::string & funcl, const std::vector<std::string> &names, const std::vector<std::string> & gmin,
            const std::vector<std::string> & gmax, const std::vector<unsigned> & nbin, bool dospline, bool usederiv,
            const std::vector<bool> &isperiodic, const std::vector<std::string> &pmin, const std::vector<std::string> &pmax);
/// get lower boundary
  std::vector<std::string> getMin() const;
/// get upper boundary
  std::vector<std::string> getMax() const;
/// get bin size
  std::vector<double> getDx() const;
  double getDx(index_t j) const ;
/// get bin volume
  double getBinVolume() const;
/// get number of bins
  std::vector<unsigned> getNbin() const;
/// get if periodic
  std::vector<bool> getIsPeriodic() const;
/// get grid dimension
  unsigned getDimension() const;
/// get argument names  of this grid
  std::vector<std::string> getArgNames() const;
/// get if the grid has derivatives
  bool hasDerivatives() const {
    return usederiv_;
  }

/// methods to handle grid indices
  void getIndices(index_t index, std::vector<unsigned>& rindex) const;
  void getIndices(index_t index, unsigned* rindex_data, std::size_t rindex_size) const;
  void getIndices(const std::vector<double> & x, std::vector<unsigned>& rindex) const;
  void getIndices(const std::vector<double> & x, unsigned* rindex_data,std::size_t rindex_size) const;
  std::vector<unsigned> getIndices(index_t index) const;
  std::vector<unsigned> getIndices(const std::vector<double> & x) const;
  index_t getIndex(const unsigned* indices,std::size_t indices_size) const;
  index_t getIndex(const std::vector<unsigned> & indices) const;
  index_t getIndex(const std::vector<double> & x) const;
  std::vector<double> getPoint(index_t index) const;
  std::vector<double> getPoint(const std::vector<unsigned> & indices) const;
  std::vector<double> getPoint(const std::vector<double> & x) const;
/// faster versions relying on preallocated vectors
  void getPoint(index_t index,std::vector<double> & point) const;
  void getPoint(const std::vector<unsigned> & indices,std::vector<double> & point) const;
  void getPoint(const std::vector<double> & x,std::vector<double> & point) const;
  void getPoint(const unsigned* indices_data,std::size_t indices_size,std::vector<double> & point) const;
  void getPoint(const unsigned* indices_data,std::size_t indices_size,double* point,std::size_t point_size) const;

/// get neighbors
  std::vector<index_t> getNeighbors(index_t index,const std::vector<unsigned> & neigh) const;
  std::vector<index_t> getNeighbors(const std::vector<unsigned> & indices,const std::vector<unsigned> & neigh) const;
  std::vector<index_t> getNeighbors(const std::vector<double> & x,const std::vector<unsigned> & neigh) const;
/// get nearest neighbors (those separated by exactly one lattice unit)
  std::vector<index_t> getNearestNeighbors(const index_t index) const;
  std::vector<index_t> getNearestNeighbors(const std::vector<unsigned> &indices) const;

/// write header for grid file
  void writeHeader(OFile& file);

/// read grid from file
  static std::unique_ptr<GridBase> create(const std::string&,const std::vector<Value*>&,IFile&,bool,bool,bool);
/// read grid from file and check boundaries are what is expected from input
  static std::unique_ptr<GridBase> create(const std::string&,const std::vector<Value*>&, IFile&,
                                          const std::vector<std::string>&,const std::vector<std::string>&,
                                          const std::vector<unsigned>&,bool,bool,bool);
/// get grid size
  virtual index_t getSize() const=0;
/// get grid value
  virtual double getValue(index_t index) const=0;
  double getValue(const std::vector<unsigned> & indices) const;
  double getValue(const std::vector<double> & x) const;
/// get grid value and derivatives
  virtual double getValueAndDerivatives(index_t index, double* der, std::size_t der_size) const=0;
  double getValueAndDerivatives(index_t index, std::vector<double>& der) const;
  double getValueAndDerivatives(const std::vector<unsigned> & indices, std::vector<double>& der) const;
  double getValueAndDerivatives(const std::vector<double> & x, std::vector<double>& der) const;

/// set grid value
  virtual void setValue(index_t index, double value)=0;
  void setValue(const std::vector<unsigned> & indices, double value);
/// set grid value and derivatives
  virtual void setValueAndDerivatives(index_t index, double value, std::vector<double>& der)=0;
  void setValueAndDerivatives(const std::vector<unsigned> & indices, double value, std::vector<double>& der);
/// add to grid value
  virtual void addValue(index_t index, double value)=0;
  void addValue(const std::vector<unsigned> & indices, double value);
/// add to grid value and derivatives
  virtual void addValueAndDerivatives(index_t index, double value, std::vector<double>& der)=0;
  void addValueAndDerivatives(const std::vector<unsigned> & indices, double value, std::vector<double>& der);
/// add a kernel function to the grid
  void addKernel( const KernelFunctions& kernel );

/// get minimum value
  virtual double getMinValue() const = 0;
/// get maximum value
  virtual double getMaxValue() const = 0;

/// dump grid on file
  virtual void writeToFile(OFile&)=0;
/// dump grid to gaussian cube file
  void writeCubeFile(OFile&, const double& lunit);

  virtual ~GridBase() = default;

/// set output format
  void setOutputFmt(const std::string & ss) {
    fmt_=ss;
  }
/// reset output format to the default %14.9f format
  void resetToDefaultOutputFmt() {
    fmt_="%14.9f";
  }
///
/// Find the maximum over paths of the minimum value of the gridded function along the paths
/// for all paths of neighboring grid lattice points from a source point to a sink point.
  double findMaximalPathMinimum(const std::vector<double> &source, const std::vector<double> &sink);
};

class Grid : public GridBase {
  std::vector<double> grid_;
  std::vector<double> der_;
  double contour_location=0.0;
public:
  Grid(const std::string& funcl, const std::vector<Value*> & args, const std::vector<std::string> & gmin,
       const std::vector<std::string> & gmax,
       const std::vector<unsigned> & nbin, bool dospline, bool usederiv):
    GridBase(funcl,args,gmin,gmax,nbin,dospline,usederiv) {
    grid_.assign(maxsize_,0.0);
    if(usederiv_) {
      der_.assign(maxsize_*dimension_,0.0);
    }
  }
/// this constructor here is not Value-aware
  Grid(const std::string& funcl, const std::vector<std::string> &names, const std::vector<std::string> & gmin,
       const std::vector<std::string> & gmax, const std::vector<unsigned> & nbin, bool dospline,
       bool usederiv, const std::vector<bool> &isperiodic, const std::vector<std::string> &pmin,
       const std::vector<std::string> &pmax ):
    GridBase(funcl,names,gmin,gmax,nbin,dospline,usederiv,isperiodic,pmin,pmax) {
    grid_.assign(maxsize_,0.0);
    if(usederiv_) {
      der_.assign(maxsize_*dimension_,0.0);
    }
  }
  index_t getSize() const override;
/// this is to access to Grid:: version of these methods (allowing overloading of virtual methods)
  using GridBase::getValue;
  using GridBase::getValueAndDerivatives;
  using GridBase::setValue;
  using GridBase::setValueAndDerivatives;
  using GridBase::addValue;
  using GridBase::addValueAndDerivatives;
/// get grid value
  double getValue(index_t index) const override;
/// get grid value and derivatives
  double getValueAndDerivatives(index_t index, double* der, std::size_t der_size) const override;
/// set grid value
  void setValue(index_t index, double value) override;
/// set grid value and derivatives
  void setValueAndDerivatives(index_t index, double value, std::vector<double>& der) override;
/// add to grid value
  void addValue(index_t index, double value) override;
/// add to grid value and derivatives
  void addValueAndDerivatives(index_t index, double value, std::vector<double>& der) override;

/// get minimum value
  double getMinValue() const override;
/// get maximum value
  double getMaxValue() const override;
/// Scale all grid values and derivatives by a constant factor
  void scaleAllValuesAndDerivatives( const double& scalef );
/// Takes the scalef times the logarithm of all grid values and derivatives
  void logAllValuesAndDerivatives( const double& scalef );
/// dump grid on file
  void writeToFile(OFile&) override;

/// Set the minimum value of the grid to zero and translates accordingly
  void setMinToZero();
/// apply function: takes  pointer to  function that accepts a double and apply
  void applyFunctionAllValuesAndDerivatives( double (*func)(double val), double (*funcder)(double valder) );
/// Get the difference from the contour
  double getDifferenceFromContour(const std::vector<double> & x, std::vector<double>& der) const ;
/// Find a set of points on a contour in the function
  void findSetOfPointsOnContour(const double& target, const std::vector<bool>& nosearch, unsigned& npoints, std::vector<std::vector<double> >& points );
/// Since this method returns a concrete Grid, it should be here and not in GridBase - GB
/// project a high dimensional grid onto a low dimensional one: this should be changed at some time
/// to enable many types of weighting
  Grid project( const std::vector<std::string> & proj, WeightBase *ptr2obj  );
  void projectOnLowDimension(double &val, std::vector<int> &varHigh, WeightBase* ptr2obj );
  void mpiSumValuesAndDerivatives( Communicator& comm );
/// Integrate the function calculated on the grid
  double integrate( std::vector<unsigned>& npoints );
  void clear();
};


class SparseGrid : public GridBase {

  std::map<index_t,double> map_;
  std::map< index_t,std::vector<double> > der_;

public:
  SparseGrid(const std::string& funcl, const std::vector<Value*> & args, const std::vector<std::string> & gmin,
             const std::vector<std::string> & gmax,
             const std::vector<unsigned> & nbin, bool dospline, bool usederiv):
    GridBase(funcl,args,gmin,gmax,nbin,dospline,usederiv) {}

  index_t getSize() const override;
  index_t getMaxSize() const;

/// this is to access to Grid:: version of these methods (allowing overloading of virtual methods)
  using GridBase::getValue;
  using GridBase::getValueAndDerivatives;
  using GridBase::setValue;
  using GridBase::setValueAndDerivatives;
  using GridBase::addValue;
  using GridBase::addValueAndDerivatives;

/// get grid value
  double getValue(index_t index) const override;
/// get grid value and derivatives
  double getValueAndDerivatives(index_t index, double* der, std::size_t der_size) const override;

/// set grid value
  void setValue(index_t index, double value) override;
/// set grid value and derivatives
  void setValueAndDerivatives(index_t index, double value, std::vector<double>& der) override;
/// add to grid value
  void addValue(index_t index, double value) override;
/// add to grid value and derivatives
  void addValueAndDerivatives(index_t index, double value, std::vector<double>& der) override;

/// get minimum value
  double getMinValue() const override;
/// get maximum value
  double getMaxValue() const override;
/// dump grid on file
  void writeToFile(OFile&) override;

  virtual ~SparseGrid() = default;
};


inline
GridBase::index_t GridBase::getIndex(const unsigned* indices,std::size_t indices_size) const {
  plumed_dbg_assert(accelerator);
  return accelerator->getIndex(*this,nbin_,indices,indices_size);
}

inline
GridBase::index_t GridBase::getIndex(const std::vector<unsigned> & indices) const {
  plumed_dbg_assert(indices.size()==dimension_);
  return getIndex(indices.data(),indices.size());
}

inline
GridBase::index_t GridBase::getIndex(const std::vector<double> & x) const {
  plumed_dbg_assert(x.size()==dimension_);
  std::array<unsigned,maxdim> indices;
  getIndices(x,indices.data(),dimension_);
  return getIndex(indices.data(),dimension_);
}

inline
std::vector<unsigned> GridBase::getIndices(index_t index) const {
  std::vector<unsigned> indices(dimension_);
  getIndices(index,indices.data(),indices.size());
  return indices;
}

inline
void GridBase::getIndices(index_t index, unsigned* indices, std::size_t indices_size) const {
  plumed_assert(indices_size==dimension_);
  plumed_assert(accelerator);
  accelerator->getIndices(nbin_,index,indices,dimension_);
}

inline
void GridBase::getIndices(index_t index, std::vector<unsigned>& indices) const {
  if (indices.size()!=dimension_) {
    indices.resize(dimension_);
  }
  getIndices(index,indices.data(),indices.size());
}

inline
std::vector<unsigned> GridBase::getIndices(const std::vector<double> & x) const {
  plumed_dbg_assert(x.size()==dimension_);
  std::vector<unsigned> indices(dimension_);
  getIndices(x,indices);
  return indices;
}

inline
void GridBase::getIndices(const std::vector<double> & x, std::vector<unsigned>& indices) const {
  indices.resize(dimension_);
  getIndices(x,indices.data(),indices.size());
}

inline
void GridBase::getIndices(const std::vector<double> & x, unsigned* rindex_data,std::size_t rindex_size) const {
  plumed_dbg_assert(x.size()==dimension_);
  plumed_dbg_assert(rindex_size==dimension_);
  accelerator->getIndices(min_,dx_,x,rindex_data,rindex_size);
}

inline
void GridBase::getPoint(const unsigned* indices,std::size_t indices_size,double* point,std::size_t point_size) const {
  plumed_dbg_assert(indices_size==dimension_);
  plumed_dbg_assert(point_size==dimension_);
  plumed_dbg_assert(accelerator);
  accelerator->getPoint(min_,dx_,indices,indices_size,point,point_size);
}

inline
std::vector<double> GridBase::getPoint(const std::vector<unsigned> & indices) const {
  plumed_dbg_assert(indices.size()==dimension_);
  std::vector<double> x(dimension_);
  getPoint(indices,x);
  return x;
}

inline
std::vector<double> GridBase::getPoint(index_t index) const {
  plumed_dbg_assert(index<maxsize_);
  std::vector<double> x(dimension_);
  getPoint(index,x);
  return x;
}

inline
std::vector<double> GridBase::getPoint(const std::vector<double> & x) const {
  plumed_dbg_assert(x.size()==dimension_);
  std::vector<double> r(dimension_);
  getPoint(x,r);
  return r;
}

inline
void GridBase::getPoint(index_t index,std::vector<double> & point) const {
  plumed_dbg_assert(index<maxsize_);
  std::array<unsigned,maxdim> indices;
  getIndices(index,indices.data(),dimension_);
  getPoint(indices.data(),dimension_,point);
}

inline
void GridBase::getPoint(const std::vector<unsigned> & indices,std::vector<double> & point) const {
  plumed_dbg_assert(indices.size()==dimension_);
  plumed_dbg_assert(point.size()==dimension_);
  getPoint(indices.data(),indices.size(),point.data(),point.size());
}

inline
void GridBase::getPoint(const unsigned* indices_data,std::size_t indices_size,std::vector<double> & point) const {
  plumed_dbg_assert(indices_size==dimension_);
  plumed_dbg_assert(point.size()==dimension_);
  getPoint(indices_data,indices_size,point.data(),point.size());
}

inline
void GridBase::getPoint(const std::vector<double> & x,std::vector<double> & point) const {
  plumed_dbg_assert(x.size()==dimension_);
  std::array<unsigned,maxdim> indices;
  getIndices(x,indices.data(),dimension_);
  getPoint(indices.data(),dimension_,point.data(),point.size());
}

inline
std::vector<GridBase::index_t> GridBase::getNeighbors(const std::vector<unsigned> &indices,const std::vector<unsigned> &nneigh)const {
  plumed_dbg_assert(accelerator);
  return accelerator->getNeighbors(*this,nbin_,pbc_,indices.data(),indices.size(),nneigh);
}

inline
std::vector<GridBase::index_t> GridBase::getNeighbors(const std::vector<double> & x,const std::vector<unsigned> & nneigh)const {
  plumed_dbg_assert(x.size()==dimension_ && nneigh.size()==dimension_);
  std::array<unsigned,maxdim> indices;
  plumed_dbg_assert(accelerator);
  accelerator->getIndices(min_,dx_,x,indices.data(),dimension_);
  return accelerator->getNeighbors(*this,nbin_,pbc_,indices.data(),dimension_,nneigh);
}

inline
std::vector<GridBase::index_t> GridBase::getNeighbors(index_t index,const std::vector<unsigned> & nneigh)const {
  plumed_dbg_assert(index<maxsize_ && nneigh.size()==dimension_);
  std::array<unsigned,maxdim> indices;
  plumed_dbg_assert(accelerator);
  accelerator->getIndices(nbin_,index,indices.data(),dimension_);
  return accelerator->getNeighbors(*this,nbin_,pbc_,indices.data(),dimension_,nneigh);
}





}

#endif
