/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2023 The plumed team
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
#include "DistanceFromContourBase.h"

namespace PLMD {
namespace contour {

void DistanceFromContourBase::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys );
  ActionWithValue::registerKeywords( keys );
  ActionAtomistic::registerKeywords( keys );
  ActionWithArguments::registerKeywords( keys );
  keys.remove("NUMERICAL_DERIVATIVES");
  keys.addInputKeyword("optional","ARG","vector","the label of the weights to use when constructing the density.  If this keyword is not here the weights are assumed to be one.");
  keys.add("atoms","POSITIONS","the positions of the atoms that we are calculating the contour from");
  keys.add("atoms","ATOM","The atom whose perpendicular distance we are calculating from the contour");
  keys.add("compulsory","BANDWIDTH","the bandwidths for kernel density esimtation");
  keys.add("compulsory","KERNEL","GAUSSIAN","the kernel function you are using.  More details on  the kernels available "
           "in plumed plumed can be found in \\ref kernelfunctions.");
  keys.add("compulsory","CUTOFF","6.25","the cutoff at which to stop evaluating the kernel functions is set equal to sqrt(2*x)*bandwidth in each direction where x is this number");
  keys.add("compulsory","CONTOUR","the value we would like for the contour");
}

DistanceFromContourBase::DistanceFromContourBase( const ActionOptions& ao ):
  Action(ao),
  ActionWithValue(ao),
  ActionAtomistic(ao),
  ActionWithArguments(ao),
  mymin(this),
  nactive(0) {
  if( getNumberOfArguments()>1 ) {
    error("should only use one argument for this action");
  }
  if( getNumberOfArguments()==1 &&  getPntrToArgument(0)->getRank()!=1 ) {
    error("ARG for distance from contour should be rank one");
  }

  // Read in the multicolvar/atoms
  std::vector<AtomNumber> atoms;
  parseAtomList("POSITIONS",atoms);
  std::vector<AtomNumber> origin;
  parseAtomList("ATOM",origin);
  if( origin.size()!=1 ) {
    error("should only specify one atom for origin keyword");
  }
  std::vector<AtomNumber> center;
  if( keywords.exists("ORIGIN") ) {
    parseAtomList("ORIGIN",center);
    if( center.size()!=1 ) {
      error("should only specify one atom for center keyword");
    }
  }

  if( center.size()==1 ) {
    log.printf("  calculating distance between atom %d and contour along vector connecting it to atom %d \n", origin[0].serial(),center[0].serial() );
  } else {
    log.printf("  calculating distance between atom %d and contour \n", origin[0].serial() );
  }

  log.printf("  contour is in field constructed from positions of atoms : ");
  for(unsigned i=0; i<atoms.size(); ++i) {
    log.printf("%d ",atoms[i].serial() );
  }
  if( getNumberOfArguments()==1 ) {
    if( getPntrToArgument(0)->getShape()[0]!=atoms.size() ) {
      error("mismatch between number of atoms and size of vector specified using ARG keyword");
    }
    log.printf("\n  and weights from %s \n", getPntrToArgument(0)->getName().c_str() );
  } else {
    log.printf("\n  all weights are set equal to one \n");
  }
  // Request everything we need
  active_list.resize( atoms.size(), 0 );
  std::vector<Value*> args( getArguments() );
  atoms.push_back( origin[0] );
  if( center.size()==1 ) {
    atoms.push_back( center[0] );
  }
  requestArguments( args );
  requestAtoms( atoms );
  // Fix to request arguments
  if( args.size()==1 ) {
    addDependency( args[0]->getPntrToAction() );
  }

  // Read in details of phase field construction
  parseVector("BANDWIDTH",bw);
  parse("KERNEL",kerneltype);
  parse("CONTOUR",contour);
  std::string errors;
  switchingFunction.set( kerneltype + " R_0=1.0 NOSTRETCH", errors );
  if( errors.length()!=0 ) {
    error("problem reading switching function description " + errors);
  }
  double det=1;
  for(unsigned i=0; i<bw.size(); ++i) {
    det*=bw[i]*bw[i];
  }
  gvol=1.0;
  if( kerneltype=="GAUSSIAN" ) {
    gvol=pow( 2*pi, 0.5*bw.size() ) * pow( det, 0.5 );
  }
  log.printf("  constructing phase field using %s kernels with bandwidth (%f, %f, %f) \n",kerneltype.c_str(), bw[0], bw[1], bw[2] );

  // And a cutoff
  std::vector<double> pp( bw.size(),0 );
  double dp2cutoff;
  parse("CUTOFF",dp2cutoff);
  double rcut =  sqrt(2*dp2cutoff)*bw[0];
  for(unsigned j=1; j<bw.size(); ++j) {
    if( sqrt(2*dp2cutoff)*bw[j]>rcut ) {
      rcut=sqrt(2*dp2cutoff)*bw[j];
    }
  }
  rcut2=rcut*rcut;

  // Create the vector of values that holds the position
  forcesToApply.resize( 3*getNumberOfAtoms() + 9 );
}

void DistanceFromContourBase::lockRequests() {
  ActionWithArguments::lockRequests();
  ActionAtomistic::lockRequests();
}

void DistanceFromContourBase::unlockRequests() {
  ActionWithArguments::unlockRequests();
  ActionAtomistic::unlockRequests();
}

double DistanceFromContourBase::evaluateKernel( const Vector& cpos, const Vector& apos, std::vector<double>& der ) const {
  Vector distance = pbcDistance( getPosition(getNumberOfAtoms()-1), cpos );
  Vector dist2 = pbcDistance( distance, apos );
  double dval=0;
  for(unsigned j=0; j<3; ++j) {
    der[j] = dist2[j]/bw[j];
    dval += der[j]*der[j];
  }
  double dfunc, newval = switchingFunction.calculateSqr( dval, dfunc ) / gvol;
  double tmp = dfunc / gvol;
  for(unsigned j=0; j<3; ++j) {
    der[j] *= -tmp;
  }
  return newval;
}

double DistanceFromContourBase::getDifferenceFromContour( const std::vector<double>& x, std::vector<double>& der ) {
  // Transer the position to the local Vector
  for(unsigned j=0; j<3; ++j) {
    pval[j] = x[j];
  }
  // Now find the contour
  double sumk = 0, sumd = 0;
  std::vector<double> pp(3), ddd(3,0);
  for(unsigned i=0; i<nactive; ++i) {
    double newval = evaluateKernel( getPosition(active_list[i]), pval, ddd );

    if( getNumberOfArguments()==1 ) {
      sumk += getPntrToArgument(0)->get(active_list[i])*newval;
      sumd += newval;
    } else {
      sumk += newval;
    }
  }
  if( getNumberOfArguments()==0 ) {
    return sumk - contour;
  }
  if( fabs(sumk)<epsilon ) {
    return -contour;
  }
  return (sumk/sumd) - contour;
}

void DistanceFromContourBase::apply() {
  if( !checkForForces() ) {
    return ;
  }
  unsigned ind=0;
  setForcesOnAtoms( getForcesToApply(), ind );
}

}
}
