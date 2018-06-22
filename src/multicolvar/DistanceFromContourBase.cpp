/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2018 The plumed team
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
#include "tools/KernelFunctions.h"
#include "DistanceFromContourBase.h"

namespace PLMD {
namespace multicolvar {

void DistanceFromContourBase::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys ); ActionWithValue::registerKeywords( keys );
  ActionAtomistic::registerKeywords( keys ); ActionWithArguments::registerKeywords( keys );
  keys.remove("NUMERICAL_DERIVATIVES");
  keys.add("atoms","POSITIONS","the positions of the atoms that we are calculating the contour from");
  keys.add("atoms","ATOM","The atom whose perpendicular distance we are calculating from the contour");
  keys.add("compulsory","BANDWIDTH","the bandwidths for kernel density esimtation");
  keys.add("compulsory","KERNEL","gaussian","the kernel function you are using.  More details on  the kernels available "
           "in plumed plumed can be found in \\ref kernelfunctions.");
  keys.add("compulsory","CONTOUR","the value we would like for the contour");
}

DistanceFromContourBase::DistanceFromContourBase( const ActionOptions& ao ):
  Action(ao),
  ActionWithValue(ao),
  ActionAtomistic(ao),
  ActionWithArguments(ao),
  nactive(0)
{
  if( getNumberOfArguments()>1 ) error("should only use one argument for this action");
  if( getNumberOfArguments()==1 ) {
    if( getPntrToArgument(0)->getRank()!=1 ) error("ARG for distance from contour should be rank one");
  }
  // Read in the multicolvar/atoms
  std::vector<AtomNumber> atoms; parseAtomList("POSITIONS",atoms);
  std::vector<AtomNumber> origin; parseAtomList("ATOM",origin);
  if( origin.size()!=1 ) error("should only specify one atom for origin keyword");
  std::vector<AtomNumber> center; 
  if( keywords.exists("ORIGIN") ) {
    parseAtomList("ORIGIN",center);
    if( center.size()!=1 ) error("should only specify one atom for center keyword");
  }

  if( center.size()==1 ) log.printf("  calculating distance between atom %d and contour along vector connecting it to atom %d \n", origin[0].serial(),center[0].serial() ); 
  else log.printf("  calculating distance between atom %d and contour \n", origin[0].serial() );

  log.printf("  contour is in field constructed from positions of atoms : ");
  for(unsigned i=0; i<atoms.size(); ++i) log.printf("%d ",atoms[i].serial() );
  if( getNumberOfArguments()==1 ) {
    if( getPntrToArgument(0)->getShape()[0]!=atoms.size() ) error("mismatch between number of atoms and size of vector specified using ARG keyword");
    log.printf("\n  and weights from %s \n", getPntrToArgument(0)->getName().c_str() );
  } else {
    log.printf("\n  all weights are set equal to one \n");
  }
  // Request everything we need
  active_list.resize( atoms.size(), 0 ); 
  std::vector<Value*> args( getArguments() ); atoms.push_back( origin[0] );
  if( center.size()==1 ) atoms.push_back( center[0] );
  requestAtoms( atoms ); requestArguments( args, false );

  // Read in details of phase field construction
  parseVector("BANDWIDTH",bw); parse("KERNEL",kerneltype); parse("CONTOUR",contour);
  log.printf("  constructing phase field using %s kernels with bandwidth (%f, %f, %f) \n",kerneltype.c_str(), bw[0], bw[1], bw[2] );

  // And a cutoff
  std::vector<double> pp( bw.size(),0 );
  KernelFunctions kernel( pp, bw, kerneltype, "DIAGONAL", 1.0 );
  double rcut = kernel.getCutoff( bw[0] );
  for(unsigned j=1; j<bw.size(); ++j) {
    if( kernel.getCutoff(bw[j])>rcut ) rcut=kernel.getCutoff(bw[j]);
  }
  rcut2=rcut*rcut;

  // Create the vector of values that holds the position
  for(unsigned i=0; i<3; ++i) pval.push_back( new Value() ); 
  forcesToApply.resize( 3*getNumberOfAtoms() + 9 );
}

DistanceFromContourBase::~DistanceFromContourBase() {
  for(unsigned i=0; i<3; ++i) delete pval[i];
}

void DistanceFromContourBase::lockRequests() {
  ActionWithArguments::lockRequests();
  ActionAtomistic::lockRequests();
}

void DistanceFromContourBase::unlockRequests() {
  ActionWithArguments::unlockRequests();
  ActionAtomistic::unlockRequests();
}

double DistanceFromContourBase::getDifferenceFromContour( const std::vector<double>& x, std::vector<double>& der ) {
  std::string min, max;
  for(unsigned j=0; j<3; ++j) {
    Tools::convert( -0.5*getBox()(j,j), min );
    Tools::convert( +0.5*getBox()(j,j), max );
    pval[j]->setDomain( min, max ); pval[j]->set( x[j] );
  }
  double sumk = 0, sumd = 0; std::vector<double> pp(3), ddd(3,0);
  for(unsigned i=0; i<nactive; ++i) {
    Vector distance = pbcDistance( getPosition(getNumberOfAtoms()-1), getPosition(active_list[i]) );
    for(unsigned j=0; j<3; ++j) pp[j] = distance[j];
  
    // Now create the kernel and evaluate
    KernelFunctions kernel( pp, bw, kerneltype, "DIAGONAL", 1.0 );
    std::vector<Value*> fvals; kernel.normalize( fvals );
    double newval = kernel.evaluate( pval, ddd, true );
    if( getNumberOfArguments()==1 ) {
      sumk += getPntrToArgument(0)->get(active_list[i])*newval;
      sumd += newval;
    } else sumk += newval;
  }
  if( getNumberOfArguments()==0 ) return sumk - contour;
  return (sumk/sumd) - contour;
}

void DistanceFromContourBase::apply() {
  if( doNotCalculateDerivatives() ) return ;
  std::vector<Vector>&   f(modifyForces());
  Tensor&           v(modifyVirial());
  const unsigned    nat=getNumberOfAtoms();

  std::fill(forcesToApply.begin(),forcesToApply.end(),0);
  if(getPntrToComponent(3)->applyForce(forcesToApply)) {
    for(unsigned j=0; j<nat; ++j) {
      f[j][0]+=forcesToApply[3*j+0];
      f[j][1]+=forcesToApply[3*j+1];
      f[j][2]+=forcesToApply[3*j+2];
    }
    v(0,0)+=forcesToApply[3*nat+0];
    v(0,1)+=forcesToApply[3*nat+1];
    v(0,2)+=forcesToApply[3*nat+2];
    v(1,0)+=forcesToApply[3*nat+3];
    v(1,1)+=forcesToApply[3*nat+4];
    v(1,2)+=forcesToApply[3*nat+5];
    v(2,0)+=forcesToApply[3*nat+6];
    v(2,1)+=forcesToApply[3*nat+7];
    v(2,2)+=forcesToApply[3*nat+8];
  }
}

}
}
