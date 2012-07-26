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
#include "FieldVessel.h"
#include "ActionWithDistribution.h"
#include "CubicInterpolation.h"

namespace PLMD {

std::string FieldVessel::documentation(){
   std::ostringstream ostr;
   ostr<<"Fields provide an alternative to collective variables \\cite field-cvs in which the instantaneous state of the system ";
   ostr<<"is represented by a function calculated on on a 1 or 2 D grid rather than by a vector of collective variables. ";
   ostr<<"Within this scheme biasing methods such as metadynamics (\\ref FIELD_METAD) can be re-formulated in terms of overlap integrals. ";
   ostr<<"Alternatively one can simply output the instaneous value of the FIELD using \\ref DUMPFIELD. ";
   ostr<<"To make calculations involving fields efficient we calculate the value of the field at a small number of grid points and then use either ";
   ostr<<"cubic or bicubic interpolation to interpolate the function onto a finer grid for the eventual integration. The input for a field looks ";
   ostr<<"something like FIELD=(MIN=\\f$r_0\\f$ MAX=\\f$x_1\\f$ NSPLINE=\\f$n\\f$ SIGMA=\\f$\\sigma\\f$) where \\f$r_0\\f$ is the lower limit for the ";
   ostr<<"integrals, \\f$x_1\\f$ is the upper limit for the integrals, \\f$n\\f$ is the number of points at which the function should be explicitally ";
   ostr<<"calculated and \\f$\\sigma\\f$ is the value of the \\f$\\sigma\\f$ parameter in the equation above.";
   return ostr.str();
}

void FieldVessel::getKeywords( Keywords& field_keys ){
  field_keys.add("compulsory","NPSPLINE","the number of points in each direction at which to calculate the value of the field");
  field_keys.add("compulsory","MIN","the minimum value at which the value of the field should be calculated in each direction");
  field_keys.add("compulsory","MAX","the maximum value at which the value of the field should be calculated in each direction");
}

FieldVessel::FieldVessel( const VesselOptions& da ):
VesselStoreAllValues(da)
{
  std::vector<std::string> data=Tools::getWords(da.parameters);
  std::vector<unsigned> nspline; std::vector<double> min,max;

  bool found_s=Tools::parseVector( data,"NSPLINE", nspline );
  if(!found_s) error("did not find NPLINE keyword");
  ndx=nspline.size(); 
  bool found_min=Tools::parseVector( data, "MIN", min );
  if(!found_min) error("did not find MIN keyword");
  if(min.size()!=ndx){
     std::string ll,ww;
     Tools::convert(ndx,ll);
     Tools::convert(min.size(),ww);
     error("found " + ww + " values for MIN when expecting only " + ll);
  }
  bool found_max=Tools::parseVector( data, "MAX", max );
  if(!found_max) error("did not find MAX keyword");
  if(max.size()!=ndx){
     std::string ll,ww;
     Tools::convert(ndx,ll);
     Tools::convert(max.size(),ww);
     error("found " + ww + " values for MAX when expecting only " + ll);
  }

  // Make sure we deal with periodicity
  if( getAction()->isPeriodic() ){
      double min, max;
      getAction()->retrieveDomain( min, max );
      tmpvalue.setDomain( min, max );
  } else {
      tmpvalue.setNotPeriodic();
  }

  // Setup the interpolators
  if( ndx==1 ) f_interpolator=new InterpolateCubic( nspline, min, max );
  else if( ndx==2 ) f_interpolator=new InterpolateBicubic( nspline, min, max );
  else plumed_assert(0);
}

FieldVessel::~FieldVessel(){
   delete f_interpolator;
   for(unsigned i=0;i<df_interpolators.size();++i) delete df_interpolators[i];
}

void FieldVessel::local_resizing(){
   // Delete the old derivative interpolators
   for(unsigned i=0;i<df_interpolators.size();++i) delete df_interpolators[i];
   df_interpolators.resize(0);

   // Find what is best to interpolate
   unsigned nvals=getAction()->getNumberOfActiveMembers();
   unsigned nder=getAction()->getNumberOfDerivatives();
   if( nvals>nder ){ ndX=nder; mergeBeforeInterpol=true; } 
   else { ndX=nvals; mergeBeforeInterpol=false; }

   // Resize the grid
   nper=(ndX+1)*(ndx+1); 
   grid_buffer.resize( f_interpolator->getNumberOfSplinePoints()*nper );
   
   // Retrieve data from the interpolator
   std::vector<unsigned> nspline(ndx); std::vector<double> min(ndx), max(ndx);
   f_interpolator->getNumbersOfPoints( nspline );
   f_interpolator->getSplinePoint( 0, min );
   f_interpolator->getSplinePoint( f_interpolator->getNumberOfSplinePoints()-1, max );
 
   // Create new derivative interpolators
   if( ndx==1 ){
       for(unsigned i=0;i<ndX;++i) df_interpolators.push_back( new InterpolateCubic( nspline, min, max ) );
   } else if( ndx==2 ){
       for(unsigned i=0;i<ndX;++i) df_interpolators.push_back( new InterpolateBicubic( nspline, min, max ) );
   } else {
       plumed_assert(0);
   }
   forces.resize(ndX);
}

std::vector<unsigned> FieldVessel::get_nspline() const {
  std::vector<unsigned> nsp(ndx);
  f_interpolator->getNumbersOfPoints( nsp );
  return nsp;
}

std::vector<double> FieldVessel::get_min() const {
  std::vector<double> min(ndx);
  f_interpolator->getSplinePoint( 0, min );
  return min; 
}

std::vector<double> FieldVessel::get_max() const {
  std::vector<double> max(ndx);
  f_interpolator->getSplinePoint( f_interpolator->getNumberOfSplinePoints()-1, max);
  return max;
}

void FieldVessel::finish( const double& tolerance ){
   unsigned stride=comm.Get_size();
   unsigned rank=comm.Get_rank();
   if(serial){ stride=1; rank=0; }
   ActionWithDistribution* aa=getAction();

   // Clear stuff from prior calculations
   grid_buffer.assign( grid_buffer.size(),0.0 );
   forces.assign( forces.size(), 0.0 );
   wasForced=false;

   unsigned kk,ik=0; bool keep; 
   std::vector<double> thisp(ndx);
   Value tmpstress; std::vector<Value> tmpder(ndX);
   tmpvalue.resizeDerivatives(ndx); tmpstress.resizeDerivatives(ndx);    
   for(unsigned nhigh=0;nhigh<ndX;++nhigh) tmpder[nhigh].resizeDerivatives(ndx);

   for(unsigned i=0;i<f_interpolator->getNumberOfSplinePoints();++i){
       f_interpolator->getSplinePoint( i, thisp );
       // Calculate the contributions of all the active colvars
       if( aa->isTimeForNeighborListUpdate() ){ keep=false; } else { keep=true; }
       for(unsigned j=0;j<aa->getNumberOfActiveMembers();++j){
           if( (ik++)%stride!=rank ) continue;  // Ensures we parallelize the double loop over nodes

           kk=aa->getActiveMember(j); getValue( kk, tmpvalue );

           // Calculate the field at point i that arises because of the jth component of the field
           calculateEnergy( kk, thisp, tmpvalue, tmpstress, tmpder );
           if( fabs( tmpstress.get() )>tolerance ){
               keep=true; unsigned nn=i*nper; 
               grid_buffer[nn]+=tmpstress.get(); nn++;
               for(unsigned nlow=0;nlow<ndx;++nlow){ grid_buffer[nn]+=tmpstress.getDerivative(nlow); nn++; } 
               for(unsigned nhigh=0;nhigh<ndX;++nhigh){
                   grid_buffer[nn]+=tmpder[nhigh].get(); nn++;
                   for(unsigned nlow=0;nlow<ndx;++nlow){ grid_buffer[nn]+=tmpder[nhigh].getDerivative(nlow); nn++; }
               }
               plumed_assert( nn==(i+1)*nper );
           }
       }
       // If the contribution of this quantity is very small at neighbour list time ignore it
       // untill next neighbour list time
       if( aa->isTimeForNeighborListUpdate() && !keep ){ aa->deactivateValue(kk); }
   }
   // Update the dynamic list 
   if( aa->isTimeForNeighborListUpdate() ){ aa->updateActiveMembers(); }
   // Accumulate the field
   if(!serial) comm.Sum( &grid_buffer[0], grid_buffer.size() ); 

   /// Transfer everything to the interpolators
   std::vector<Value> vv( f_interpolator->getNumberOfSplinePoints() );
   for(unsigned i=0;i<f_interpolator->getNumberOfSplinePoints();++i){
      vv[i].set( grid_buffer[i*nper] );
      vv[i].resizeDerivatives( ndx );
      for(unsigned j=0;j<ndx;++j) vv[i].addDerivative( j, grid_buffer[ i*nper + j + 1 ] );
   }
   f_interpolator->set_table( vv );

   unsigned nn;
   for(unsigned n=0;n<df_interpolators.size();++n){
      for(unsigned i=0;i<f_interpolator->getNumberOfSplinePoints();++i){
         vv[i].clearDerivatives();
         nn=i*nper + (n+1)*(1+ndx); vv[i].set( grid_buffer[nn] ); nn++;
         for(unsigned j=0;j<ndx;++j){ vv[i].addDerivative( j, grid_buffer[nn] ); nn++; }
      }
      df_interpolators[n]->set_table( vv );
   }
}

void FieldVessel::mergeFieldDerivatives( const std::vector<double>& inder, Value* outval ){
   plumed_assert( inder.size()==ndX );

   if(mergeBeforeInterpol){
       for(unsigned i=0;i<ndX;++i) outval->addDerivative( i, inder[i] ); 
   } else {
       Value tmpder, tmpder2; 
       ActionWithDistribution* aa=getAction();
       for(unsigned i=0;i<ndX;++i){
           getValue( i, tmpder ); aa->mergeDerivatives( i, tmpder, inder[i], tmpder2 );
           for(unsigned j=0;j<tmpder2.getNumberOfDerivatives();++j) outval->addDerivative( j, tmpder2.getDerivative(j) );
       }
   }
}

double FieldVessel::getFieldAt( const std::vector<double>& pp ){
   plumed_assert( pp.size()==ndx );
   return f_interpolator->get_fdf(pp);
}

void FieldVessel::getDerivativesAt( const std::vector<double>& pp, std::vector<double>& der ){
  plumed_assert( pp.size()==ndx && der.size()==ndX );
  for(unsigned i=0;i<ndX;++i) der[i]=df_interpolators[i]->get_fdf(pp);
}

void FieldVessel::addForces( const std::vector<double>& inforces ){
  plumed_assert( inforces.size()==forces.size() );
  wasForced=true;
  for(unsigned i=0;i<ndX;++i) forces[i]=inforces[i];
}

bool FieldVessel::applyForce( std::vector<double>& outforces ){
  if(!wasForced) return false;

  if( mergeBeforeInterpol ){
      plumed_assert( outforces.size()==ndX );
      for(unsigned i=0;i<forces.size();++i) outforces[i]=forces[i];
  } else {
      Value tmpder, tmpder2; 
      ActionWithDistribution* aa=getAction();
      for(unsigned i=0;i<ndX;++i){
          getValue( i, tmpder ); aa->mergeDerivatives( i, tmpder, outforces[i], tmpder2 );
          plumed_assert( outforces.size()==tmpder2.getNumberOfDerivatives() );
          for(unsigned j=0;j<tmpder2.getNumberOfDerivatives();++j) outforces[j]=tmpder2.getDerivative(j);
      }
  }
  return true;
}

}
