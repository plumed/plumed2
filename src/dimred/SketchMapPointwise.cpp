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
#include "core/ActionRegister.h"
#include "SketchMapBase.h"
#include "tools/ConjugateGradient.h"
#include "tools/GridSearch.h"

//+PLUMEDOC DIMRED SKETCHMAP_POINTWISE
/*
Optimise the sketch-map stress function using a pointwise global optimisation algorithm.

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace dimred {

class SketchMapPointwise : public SketchMapBase {
private:
  unsigned ncycles;
  bool nointerpolation;
  double cgtol, gbuf;
  std::vector<unsigned> npoints;
public:
  static void registerKeywords( Keywords& keys ); 
  SketchMapPointwise( const ActionOptions& ao );
  void minimise( const Matrix<double>& , const Matrix<double>& , Matrix<double>& );
};

PLUMED_REGISTER_ACTION(SketchMapPointwise,"SKETCHMAP_POINTWISE")

void SketchMapPointwise::registerKeywords( Keywords& keys ){
  SketchMapBase::registerKeywords( keys );
  keys.add("compulsory","NCYCLES","5","the number of cycles of global optimisation to attempt");
  keys.add("compulsory","CGTOL","1E-6","the tolerance for the conjugate gradient minimisation");
  keys.add("compulsory","BUFFER","1.1","grid extent for search is (max projection - minimum projection) multiplied by this value");
  keys.add("compulsory","NGRID","10","number of points to use in each grid direction");
  keys.addFlag("NOINTERPOLATION",false,"do not use any interpolation to make the optimisation faster");
}

SketchMapPointwise::SketchMapPointwise( const ActionOptions& ao ):
Action(ao),
SketchMapBase(ao),
npoints(nlow)
{
  parseVector("NGRID",npoints); parse("BUFFER",gbuf);
  if( npoints.size()!=nlow ) error("vector giving number of grid point in each direction has wrong size");
  parse("NCYCLES",ncycles); parse("CGTOL",cgtol); 
  parseFlag("NOINTERPOLATION",nointerpolation);

  log.printf("  doing %d cycles of global optimisation sweeps\n",ncycles);
  log.printf("  using grid of points that is %d",npoints[0]);
  for(unsigned j=1;j<npoints.size();++j) log.printf(" by %d",npoints[j]);
  log.printf(" and that is %f larger than the difference between the position of the minimum and maximum projection \n",gbuf);
  log.printf("  tolerance for conjugate gradient algorithm equals %f \n",cgtol);
  if( nointerpolation ) log.printf("  not using any interpolation in grid search\n");
}

void SketchMapPointwise::minimise( const Matrix<double>& transformed, const Matrix<double>& distances, Matrix<double>& projections ){
  std::vector<double> gmin( nlow ), gmax( nlow ), mypoint( nlow ); 

  // Find the extent of the grid
  for(unsigned j=0;j<nlow;++j) gmin[j]=gmax[j]=projections(0,j); 
  for(unsigned i=1;i<getNumberOfDataPoints();++i){
      for(unsigned j=0;j<nlow;++j){
          if( projections(i,j) < gmin[j] ) gmin[j] = projections(i,j); 
          if( projections(i,j) > gmax[j] ) gmax[j] = projections(i,j); 
      }
  }
  for(unsigned j=0;j<nlow;++j){ 
     double gbuffer = 0.5*gbuf*( gmax[j]-gmin[j] ) - 0.5*( gmax[j]- gmin[j] ); 
     gmin[j]-=gbuffer; gmax[j]+=gbuffer; 
  }

  // And do the search
  ConjugateGradient<SketchMapPointwise> mycgminimise( this );
  GridSearch<SketchMapPointwise> mygridsearch( gmin, gmax, npoints, this );
  // Run multiple loops over all projections
  for(unsigned i=0;i<ncycles;++i){
      for(unsigned j=0;j<getNumberOfDataPoints();++j){
          // Setup target distances and target functions for calculate stress 
          for(unsigned k=0;k<getNumberOfDataPoints();++k) setTargetDistance( k, distances(j,k)  ); 

          if( nointerpolation ){
              // Minimise using grid search
              mygridsearch.minimise( mypoint, &SketchMapPointwise::calculateStress );
              // Minimise output using conjugate gradient
              mycgminimise.minimise( cgtol, mypoint, &SketchMapPointwise::calculateStress ); 
          } else {
              plumed_merror("use a class here that Sandip will implement that allows you to do this with some interpolation");
          }
          // Save new projection 
          for(unsigned k=0;k<nlow;++k) projections(j,k)=mypoint[k];
      }
  }
}

}
}
