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
#ifndef __PLUMED_tools_KernelFunctions_h
#define __PLUMED_tools_KernelFunctions_h

#include "Matrix.h"
#include "core/Value.h"
#include <vector>
#include <memory>

namespace PLMD {

class KernelFunctions {
private:
/// Is the metric matrix diagonal
  enum {diagonal,multi,vonmises} dtype;
/// What type of kernel are we using
  enum {gaussian,truncatedgaussian,uniform,triangular} ktype;
/// The center of the kernel function
  std::vector<double> center;
/// The width of the kernel
  std::vector<double> width;
/// The height of the kernel
  double height;
/// Used to set all the data in the kernel during construction - avoids double coding as this has two constructors
  void setData( const std::vector<double>& at, const std::vector<double>& sig, const std::string& type, const std::string& mtype, const double& w );
/// Convert the width into matrix form
  Matrix<double> getMatrix() const;
public:
  explicit KernelFunctions( const std::string& input );
  KernelFunctions( const std::vector<double>& at, const std::vector<double>& sig, const std::string& type, const std::string& mtype, const double& w );
  explicit KernelFunctions( const KernelFunctions* in );
/// Normalise the function and scale the height accordingly
  void normalize( const std::vector<Value*>& myvals );
/// Get the dimensionality of the kernel
  unsigned ndim() const;
/// Get the cutoff for a kernel
  double getCutoff( const double& width ) const ;
/// Get the position of the center
  std::vector<double> getCenter() const;
/// Get the support
  std::vector<unsigned> getSupport( const std::vector<double>& dx ) const;
/// get it in continuous form
  std::vector<double> getContinuousSupport( ) const;
/// Evaluate the kernel function with constant intervals
  double evaluate( const std::vector<Value*>& pos, std::vector<double>& derivatives, bool usederiv=true, bool doInt=false, double lowI_=-1, double uppI_=-1 ) const;
/// Read a kernel function from a file
  static std::unique_ptr<KernelFunctions> read( IFile* ifile, const bool& cholesky, const std::vector<std::string>& valnames );
};

inline
Matrix<double> KernelFunctions::getMatrix() const {
  unsigned k=0, ncv=ndim(); Matrix<double> mymatrix(ncv,ncv);
  for(unsigned i=0; i<ncv; i++) {
    for(unsigned j=i; j<ncv; j++) {
      mymatrix(i,j)=mymatrix(j,i)=width[k]; // recompose the full inverse matrix
      k++;
    }
  }
  return mymatrix;
}

inline
unsigned KernelFunctions::ndim() const {
  return center.size();
}

inline
std::vector<double> KernelFunctions::getCenter() const {
  return center;
}

}
#endif
