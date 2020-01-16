/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2019 The plumed team
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
  void 		pushKernel( IFile * ff);
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
  std::vector<std::unique_ptr<KernelFunctions>> hills;
  std::vector<double> biasf;
  std::vector<double> histosigma;
  Communicator& mycomm;
  std::unique_ptr<Grid> BiasGrid_;
};

}

#endif
