/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

   This file is part of plumed, version 2.0.

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
#include "SketchMapBase.h"

namespace PLMD {
namespace dimred {

void SketchMapBase::registerKeywords( Keywords& keys ){
  DimensionalityReductionBase::registerKeywords( keys );
  keys.remove("NLOW_DIM");
  keys.add("compulsory","HIGH_DIM_FUNCTION","the parameters of the switching function in the high dimensional space");
  keys.add("compulsory","LOW_DIM_FUNCTION","the parameters of the switching function in the low dimensional space");
  keys.add("compulsory","MIXPARAM","0.0","the ammount of the pure distances to mix into the stress function");
}

SketchMapBase::SketchMapBase( const ActionOptions& ao ):
Action(ao),
DimensionalityReductionBase(ao)
{
  // Read in the switching functions
  std::string linput,hinput, errors;
  parse("HIGH_DIM_FUNCTION",hinput);
  highdf.set(hinput,errors);
  if(errors.length()>0) error(errors);
  parse("LOW_DIM_FUNCTION",linput);
  lowdf.set(hinput,errors);
  if(errors.length()>0) error(errors);
  log.printf("  filter function for dissimilarities in high dimensional space has cutoff %s \n",highdf.description().c_str() );
  log.printf("  filter function for distances in low dimensionality space has cutoff %s \n",highdf.description().c_str() );

  // Read the mixing parameter
  parse("MIXPARAM",mixparam);
  if( mixparam<0 || mixparam>1 ) error("mixing parameter must be between 0 and 1");
  log.printf("  mixing %f of pure distances with %f of filtered distances \n",mixparam,1.-mixparam);
}

void SketchMapBase::calculateProjections( const Matrix<double>& targets, Matrix<double>& projections ){
  Matrix<double> transformed( targets.nrows(), targets.ncols() );
  Matrix<double> distances( targets.nrows(), targets.ncols() ); 

  // Transform the high dimensional distances
  double df;
  for(unsigned i=1;i<distances.ncols();++i){
      for(unsigned j=0;j<i;++j){
          distances(i,j)=distances(j,i)=sqrt( targets(i,j) );
          transformed(i,j)=transformed(j,i)=1.0 - highdf.calculate( distances(i,j), df );
      }
  }
  // And minimse
  minimise( transformed, distances, projections );
}

}
}


