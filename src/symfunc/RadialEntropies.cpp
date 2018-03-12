/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2017 The plumed team
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
#include "SymmetryFunctionBase.h"
#include "multicolvar/MultiColvarBase.h"
#include "core/ActionRegister.h"
#include "tools/KernelFunctions.h"
#include <string>
#include <cmath>

using namespace std;

namespace PLMD {
namespace symfunc {

//+PLUMEDOC MCOLVAR RADIAL_ENTROPY
/*

\par Examples

*/
//+ENDPLUMEDOC


class RadialEntropy : public SymmetryFunctionBase {
  int support;
  unsigned nbins;
  double int_max, density, spacing;
  std::vector<KernelFunctions> mykernel;
public:
  static void registerKeywords( Keywords& keys );
  explicit RadialEntropy(const ActionOptions&);
  void compute( const double& weight, const Vector& vec, MultiValue& myvals ) const { plumed_error(); }
  void computeSymmetryFunction( const unsigned& current, MultiValue& myvals ) const ;
};

PLUMED_REGISTER_ACTION(RadialEntropy,"RADIAL_ENTROPY")

void RadialEntropy::registerKeywords( Keywords& keys ) {
  SymmetryFunctionBase::registerKeywords( keys ); keys.remove("ONESHOT");
  keys.add("compulsory","UPPER_LIMIT","the upper limit for the integration of the radial distribution function");
  keys.add("compulsory","NBINS","the number of bins to use for the numerical integration of the rdf");
  keys.add("compulsory","KERNEL","gaussian","the type of kernel to use in the kernel density estimation");
  keys.add("compulsory","BANDWIDTH","the bandwidth to use for the gaussian kernels");
  keys.add("compulsory","DENSITY","the density of the system");
}

RadialEntropy::RadialEntropy(const ActionOptions&ao):
  Action(ao),
  SymmetryFunctionBase(ao)
{
  std::string ktype; parse("KERNEL",ktype); double bw; parse("BANDWIDTH",bw);
  std::vector<double> zero(1,0), bandw(1); bandw[0]=bw;
  mykernel.push_back( KernelFunctions( zero, bandw, ktype, "DIAGONAL", 1.0 ) );
  std::vector<Value*> fvals; mykernel[0].normalize( fvals );
  log.printf("  construcing rdf using %s kernels with bandwidth %f \n", ktype.c_str(), bw );
  parse("UPPER_LIMIT",int_max); parse("NBINS",nbins);
  log.printf("  computing entropy by integrating from 0 to %f using %d integration bins \n",int_max, nbins );
  std::vector<double> bin_size_vec(1); bin_size_vec[0] = spacing = int_max / static_cast<double>( nbins );
  std::vector<unsigned> sup( mykernel[0].getSupport( bin_size_vec ) ); support=sup[0];
  parse("DENSITY",density); log.printf("  density of system is set at %f \n",density );
  addValueWithDerivatives();
}

void RadialEntropy::computeSymmetryFunction( const unsigned& current, MultiValue& myvals ) const {
  unsigned matind = getPntrToArgument(0)->getPositionInMatrixStash();
  unsigned matind_x = getPntrToArgument(1)->getPositionInMatrixStash();
  unsigned matind_y = getPntrToArgument(2)->getPositionInMatrixStash();
  unsigned matind_z = getPntrToArgument(3)->getPositionInMatrixStash();
  std::vector<Value*> pval; pval.push_back( new Value() ); pval[0]->setNotPeriodic();
  std::vector<double> der(1); std::vector<double> grid( nbins, 0 );
  for(unsigned i=0; i<myvals.getNumberOfStashedMatrixElements(matind); ++i) {
    unsigned iind = myvals.getStashedMatrixIndex(matind,i); Vector dist;
    double weight=myvals.getStashedMatrixElement( matind, iind );
    dist[0]=myvals.getStashedMatrixElement( matind_x, iind );
    dist[1]=myvals.getStashedMatrixElement( matind_y, iind );
    dist[2]=myvals.getStashedMatrixElement( matind_z, iind );
    unsigned ind = std::floor( dist.modulo() / spacing );
    for(int j=-support; j<=support; ++j) {
      pval[0]->set( j*spacing ); grid[ind+j] += weight*mykernel[0].evaluate( pval, der );
    }
  }
  delete pval[0];
  // Now do integral
  for(unsigned i=1; i<nbins; ++i) {
    double r = i*spacing; grid[i] = grid[i] / (4*pi*density*r*r );
    addToValue( 0, -2*pi*density*spacing*r*r*(grid[i]*( std::log( grid[i] ) - 1 ) + 1 ), myvals );
  }
}

}
}
