/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2014 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

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
#include "core/ActionRegister.h"
#include "tools/HistogramBead.h"
#include "Gradient.h"

namespace PLMD {
namespace crystallisation {

PLUMED_REGISTER_ACTION(Gradient,"GRADIENT")

void Gradient::registerKeywords( Keywords& keys ){
  VolumeGradientBase::registerKeywords( keys );
  keys.add("atoms","ORIGIN","we will use the position of this atom as the origin in our calculation");
  keys.add("compulsory","DIR","xyz","the directions in which we are calculating the graident.  Should be x, y, z, xy, xz, yz or xyz");
  keys.add("compulsory","NBINS","number of bins to use in each direction for the calculation of the gradient");
  keys.add("compulsory","SIGMA","1.0","the width of the function to be used for kernel density estimation");
  keys.add("compulsory","KERNEL","gaussian","the type of kernel function to be used");
}

Gradient::Gradient(const ActionOptions&ao):
Action(ao),
VolumeGradientBase(ao),
nbins(3)
{
  std::vector<AtomNumber> atom;
  parseAtomList("ORIGIN",atom);
  if( atom.size()!=1 ) error("should only be one atom specified");
  log.printf("  origin is at position of atom : %d\n",atom[0].serial() );

  std::string direction; parse("DIR",direction);
  std::vector<unsigned> tbins; parseVector("NBINS",tbins);
  for(unsigned i=0;i<tbins.size();++i){
     if( tbins[i]<2 ) error("Number of grid points should be greater than 1");
  }

  if( direction=="x" ){
      if( tbins.size()!=1 ) error("mismatch between number of bins and direction");
      nbins[0]=tbins[0]; nbins[1]=0; nbins[2]=0;
  } else if( direction=="y" ){
      if( tbins.size()!=1 ) error("mismatch between number of bins and direction");
      nbins[0]=0; nbins[1]=tbins[0]; nbins[2]=0;
  } else if( direction=="z" ){
      if( tbins.size()!=1 ) error("mismatch between number of bins and direction");
      nbins[0]=0; nbins[1]=0; nbins[2]=tbins[0];
  } else if( direction=="xy" ){
      if( tbins.size()!=2 ) error("mismatch between number of bins and direction");
      nbins[0]=tbins[0]; nbins[1]=tbins[1]; nbins[2]=0;
  } else if( direction=="xz" ){
      if( tbins.size()!=2 ) error("mismatch between number of bins and direction");
      nbins[0]=tbins[0]; nbins[1]=0; nbins[2]=tbins[1];
  } else if( direction=="yz" ){
      if( tbins.size()!=2 ) error("mismatch between number of bins and direction");
      nbins[0]=0; nbins[1]=tbins[0]; nbins[2]=tbins[1];
  } else if( direction=="xyz" ){
      if( tbins.size()!=3 ) error("mismatch between number of bins and direction");
      nbins[0]=tbins[0]; nbins[1]=tbins[1]; nbins[2]=tbins[2];
  } else {
     error( direction + " is not valid gradient direction");
  }

  // Find number of quantities
  if( getPntrToMultiColvar()->isDensity() ) vend=1;
  else if( getPntrToMultiColvar()->getNumberOfQuantities()==5 ) vend=1;
  else vend= 1 + getPntrToMultiColvar()->getNumberOfQuantities()-5; // +1 is for weight 
  nquantities = vend + nbins[0] + nbins[1] + nbins[2];

  // Output some nice information
  std::string functype=getPntrToMultiColvar()->getName();
  std::transform( functype.begin(), functype.end(), functype.begin(), tolower );
  log.printf("  calculating gradient of %s in %s direction \n",functype.c_str(), direction.c_str() ); 

  parse("SIGMA",sigma); 
  bead.isNotPeriodic(); std::string kerneltype; 
  parse("KERNEL",kerneltype); bead.setKernelType( kerneltype );
  checkRead(); requestAtoms(atom);

  // And setup the vessel
  std::string input; addVessel( "GRADIENT", input, -1 );
  // And resize atoms
  finishTaskListUpdate();
}

unsigned Gradient::getCentralAtomElementIndex(){
  plumed_error();
}

void Gradient::setupRegions(){
//  if( !getPbc().isOrthorombic() ) error("cell must be orthorhombic when using gradient");
}

void Gradient::calculateAllVolumes(){
  Vector cpos = pbcDistance( getPosition(0), getPntrToMultiColvar()->retrieveCentralAtomPos() );
  Vector oderiv, fpos = getPbc().realToScaled( cpos );  
  
  Vector deriv; unsigned nbase=vend;
  for(unsigned idir=0;idir<3;++idir){
      deriv[0]=deriv[1]=deriv[2]=0.0;
      double delx = 1.0 / static_cast<double>( nbins[idir] );
      for(unsigned jbead=0;jbead<nbins[idir];++jbead){
          // Calculate what box we are in
          bead.set( -0.5+jbead*delx, -0.5+(jbead+1)*delx, sigma );
          double weight=bead.calculate( fpos[0], deriv[idir] );
          oderiv = getPbc().realToScaled( deriv );          
          // Set and derivatives
          setNumberInVolume( nbase+jbead, weight, oderiv ); 
          addReferenceAtomDerivatives( nbase+jbead, 0, -oderiv );
//          addBoxDerivatives( nbase+jbead, -Tensor(cpos,oderiv) );
      }
      nbase += nbins[idir];
  }
}

}
}
