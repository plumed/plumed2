/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2023 The plumed team
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
#ifndef __PLUMED_tools_BiasRepresentation_h
#define __PLUMED_tools_BiasRepresentation_h

#include "Exception.h"
#include <limits>
#include <memory>
#include <vector>

namespace PLMD {

class Value;
class Grid;
class IFile;
class KernelFunctions;
class Communicator;

//+PLUMEDOC INTERNAL biasrepresentation
/*

*/
//+ENDPLUMEDOC

/// this class implements a general purpose class that aims to
/// provide a Grid/list
/// transparently add gaussians to a bias

class BiasRepresentation {
public:
  /// create a bias representation from a list of pointer to values
  BiasRepresentation(const std::vector<Value*> & tmpvalues, Communicator &cc  );
  /// create a bias using explicit sigma in input (needed for histogram building)
  BiasRepresentation(const std::vector<Value*> & tmpvalues, Communicator &cc,  const std::vector<double> & sigma);
  /// create a bias containing a grid representation
  BiasRepresentation(const std::vector<Value*> & tmpvalues, Communicator &cc, const std::vector<std::string> &  gmin, const std::vector<std::string> & gmax,
                     const std::vector<unsigned> & nbin, bool doInt, double lowI_, double uppI_);
  /// create a histogram with grid representation and sigmas in input
  BiasRepresentation(const std::vector<Value*> & tmpvalues, Communicator &cc, const std::vector<std::string> & gmin, const std::vector<std::string> & gmax, const std::vector<unsigned> & nbin, const std::vector<double> & sigma);
  /// retrieve the number of dimension of the representation
  unsigned 	getNumberOfDimensions();
  /// add the grid to the representation
  void 		addGrid(const std::vector<std::string> & gmin, const std::vector<std::string> & gmax, const std::vector<unsigned> & nbin );
  /// push a kernel on the representation (includes widths and height)
  void 		pushKernel( IFile * ff, double time=0.0 );
  /// read a kernel from a file without adding it; bf is the bias factor of the hill (1 if not in the file)
  std::unique_ptr<KernelFunctions> readKernel( IFile * ff, double & bf );
  /// add a kernel previously read with readKernel
  void 		addKernel( std::unique_ptr<KernelFunctions> kk, double bf, double time=0.0 );
  /// time average of the bias between t0 and t1 (if hasend is false: until the last hill), obtained by scaling the hills
  /// deposited between t0 and t1 by (t1-t)/(t1-t0) and ignoring those after t1
  void setTimeAverage(double t0, bool hasend, double t1);
  /// the end of the time average (t1, or the time of the last hill if earlier or if t1 was not given)
  double getTimeAverageEnd() const;
  /// the largest time of the kernels added
  double getLastTime() const;
  /// set the flag that rescales the free energy to the bias
  void 		setRescaledToBias(bool rescaled);
  /// check if the representation is rescaled to the bias
  const bool & 	isRescaledToBias();
  /// check if the sigma values are already provided (in case of a histogram representation with input sigmas)
  bool  	hasSigmaInInput();
  /// get the names of the variables
  std::vector<std::string> getNames();
  /// get the pointer to the values
  const std::vector<Value*> & getPtrToValues();
  /// get the number of kernels contained in the representation
  int 		getNumberOfKernels();
  /// get the name of the i-th value
  const std::string & getName(unsigned i);
  /// get a pointer to a specific value
  Value* 	getPtrToValue(unsigned i);
  /// get the pointer to the grid
  Grid* 	getGridPtr();
  /// get a new histogram point from a file
  std::unique_ptr<KernelFunctions> readFromPoint(IFile *ifile);
  /// get an automatic min/max from the set so to know how to configure the grid
  void getMinMaxBin(std::vector<double> &vmin, std::vector<double> &vmax, std::vector<unsigned> &vbin);
  /// clear the representation (grid included)
  void clear();
private:
  int ndim;
  bool hasgrid;
  bool rescaledToBias;
  bool doInt_;
  double lowI_;
  double uppI_;
  std::vector<Value*> values;
  std::vector<std::string> names;
  /// number of kernels added
  int nkernels=0;
  /// support of the kernels added without a grid (see getMinMaxBin)
  std::vector<double> supmin,supmax,minbinsize;
  void resetSupport();
  std::vector<double> histosigma;
  Communicator& mycomm;
  std::unique_ptr<Grid> BiasGrid_;
  /// with MPI, the sum of the grids of all the processes
  std::unique_ptr<Grid> reducedGrid_;
  /// geometry of the grid
  std::vector<unsigned> gridnbin;
  std::vector<bool> gridpbc;
  std::vector<double> griddx,gridmin;
  std::vector<std::size_t> gridstride;
  /// work arrays of addSeparableGaussian: for each dimension the grid indices, (x-c)/w, its square and exp(-((x-c)/w)^2/2)
  std::vector<std::vector<unsigned>> sepIndex;
  std::vector<std::vector<double>> sepT,sepT2,sepExp;
  /// time average
  bool timeaverage=false;
  double avt0=0.0;
  double avt1=0.0;
  bool avhasend=false;
  double lasttime=-std::numeric_limits<double>::max();
  /// hills after t0, summed with weights 1 and t
  std::unique_ptr<Grid> avS1_,avS2_;
  void addToGrid( Grid & g, const KernelFunctions & kk, double f );
  void addSeparableGaussian( Grid & g, const KernelFunctions & kk, double h, double a, double b, double r2cut, const std::vector<double> & w );
  void addGenericKernel( Grid & g, const KernelFunctions & kk, double f );
};

}

#endif
