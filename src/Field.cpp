#include "Field.h"
#include "PlumedCommunicator.h"
#include "Keywords.h"


using namespace PLMD;

std::string Field::documentation(){
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

Field::Field( const std::string ftype, const unsigned d ){
  if( ftype=="identity") fstyle=identity;
  else if( ftype=="gaussian") fstyle=gaussian;
  else errormsg = "type " + ftype + " is not implemented in field";
  ndx=d;
}

Field::~Field(){
  delete f_interpolator;
  for(unsigned i=0;i<df_interpolators.size();++i) delete df_interpolators[i];
}

void Field::read( const std::string& input, const unsigned nfunc, std::string& report ){
  std::vector<std::string> data=Tools::getWords(input);
  std::vector<unsigned> nspline; std::vector<double> min,max; 

  bool found_s=Tools::parseVector( data,"NSPLINE", nspline );
  if(!found_s) error("did not find NPLINE keyword");
  if(nspline.size()!=ndx){
     std::string ll,ww;
     Tools::convert(ndx,ll);
     Tools::convert(nspline.size(),ww);
     error("found " + ww + " values for NSPLINE when expecting only " + ll);
  }
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
  bool found_sigma=Tools::parse( data, "SIGMA", sigma );
  if(!found_sigma) error("did not find SIGMA keyword");

  if( data.size()!=0 ){
      std::string err="found the following rogue keywords : ";
      for(unsigned i=0;i<data.size();++i) err=err + data[i] + ", ";
      error(err);
  }
  std::ostringstream ostr;
  ostr<<"generating "<<ndx<<" dimensional field min : ";
  for(unsigned i=0;i<min.size();++i) ostr<<min[i]<<" ";
  ostr<<"max : ";
  for(unsigned i=0;i<max.size();++i) ostr<<max[i]<<" ";
  ostr<<"nspline : ";
  for(unsigned i=0;i<nspline.size();++i) ostr<<nspline[i];
  report=ostr.str();

// Set everything for base quantities 
  baseq_nder.resize(nfunc); baseq_starts.resize(nfunc);
// Set everything for grid
  unsigned np=1; for(unsigned i=0;i<nspline.size();++i) np*=nspline[i];
  npoints=np; 

  // Setup the interpolators
  if( ndx==1 ) f_interpolator=new InterpolateCubic( nspline, min, max );
  else if( ndx==2 ) f_interpolator=new InterpolateBicubic( nspline, min, max );
  else plumed_assert(0);
}

void Field::resizeDerivatives( const unsigned D ){
  ndX=D; nper=(ndX+1)*(ndx+1); grid_buffer.resize( npoints*nper );

  for(unsigned i=0;i<df_interpolators.size();++i) delete df_interpolators[i];
  df_interpolators.resize(0);
  std::vector<unsigned> nspline(ndx); std::vector<double> min(ndx), max(ndx);
  f_interpolator->getNumbersOfPoints( nspline ); 
  f_interpolator->getSplinePoint( 0, min );
  f_interpolator->getSplinePoint( f_interpolator->getNumberOfSplinePoints()-1, max );
  if( ndx==1 ){
      for(unsigned i=0;i<D;++i) df_interpolators.push_back( new InterpolateCubic( nspline, min, max ) ); 
  } else if( ndx==2 ){ 
      for(unsigned i=0;i<D;++i) df_interpolators.push_back( new InterpolateBicubic( nspline, min, max ) );
  } else {
      plumed_assert(0);
  }
  forces.resize(D);
}

void Field::retrieveBoundaries( std::vector<double>& min, std::vector<double>& max ){
  min.resize(ndx); max.resize(ndx);
  f_interpolator->getSplinePoint( 0, min ); 
  f_interpolator->getSplinePoint( f_interpolator->getNumberOfSplinePoints()-1, max );
}

void Field::clear(){
  wasforced=false;
  baseq_buffer.assign( baseq_buffer.size(), 0.0 );
  grid_buffer.assign( grid_buffer.size(), 0.0 );
}

void Field::get_nspline( std::vector<unsigned>& nspline ) const {
  plumed_assert( ndx==1 );
  nspline[0]=npoints;
}

void Field::resizeBaseQuantityBuffers( const std::vector<unsigned>& cv_sizes ){
  baseq_nder.resize( cv_sizes.size() ); baseq_starts.resize( cv_sizes.size() );
  unsigned nn=0;
  for(unsigned i=0;i<cv_sizes.size();++i){
      baseq_starts[i]=nn; baseq_nder[i]=cv_sizes[i]; nn+=cv_sizes[i]+1;
  }
  baseq_buffer.resize(nn); // And resize the buffers
}

void Field::setBaseQuantity( const unsigned nn, Value* val ){
  plumed_assert( nn<baseq_nder.size() );
  plumed_assert( val->getNumberOfDerivatives()==baseq_nder[nn] );

  unsigned kk=baseq_starts[nn];
  baseq_buffer[kk]=val->get(); kk++;
  for(unsigned i=0;i<val->getNumberOfDerivatives();++i){ baseq_buffer[kk]=val->getDerivative(i); kk++; }
  if( (nn+1)==baseq_starts.size() ){ plumed_assert( kk==baseq_buffer.size() ); }
  else{ plumed_assert( kk==baseq_starts[nn+1] ); }
}

void Field::gatherBaseQuantities( PlumedCommunicator& comm ){
  comm.Sum( &baseq_buffer[0],baseq_buffer.size() );
}

void Field::getSplinePoint( const unsigned nn, std::vector<double>& pp ) const {
  f_interpolator->getSplinePoint( nn, pp );
}

void Field::extractBaseQuantity( const unsigned nn, Value* val ){
  plumed_assert( nn<baseq_nder.size() );
  if( baseq_nder[nn]!=val->getNumberOfDerivatives() ) val->resizeDerivatives(baseq_nder[nn]);

  val->clearDerivatives();
  unsigned kk=baseq_starts[nn];
  val->set( baseq_buffer[kk] ); kk++;
  for(unsigned i=0;i<val->getNumberOfDerivatives();++i){ val->addDerivative( i, baseq_buffer[kk] ); kk++; }
  if( (nn+1)==baseq_starts.size() ){ plumed_assert( kk==baseq_buffer.size() ); }
  else{ plumed_assert( kk==baseq_starts[nn+1] ); }
}

double Field::calculateField( const std::vector<double>& pp ) const {
  if( fstyle==identity ){
     return f_interpolator->get_fdf(pp);
  } else if( fstyle==gaussian ) {
     double tmp=f_interpolator->get_fdf(pp);
     return exp( -tmp/(2*sigma*sigma) );
  }
  plumed_massert(0, "no field calculate style defined");
  return 0;
}

void Field::calculateFieldDerivatives( const std::vector<double>& pp, std::vector<double>& der ) const {
  plumed_assert( der.size()==ndX );
  if( fstyle==identity ){
     for(unsigned i=0;i<der.size();++i) der[i]=df_interpolators[i]->get_fdf(pp);
  } else if( fstyle==gaussian ) {
     double tmp=f_interpolator->get_fdf(pp); double pref=-exp( -tmp/(2*sigma*sigma) ) / ( 2*sigma*sigma );  
     for(unsigned i=0;i<der.size();++i) der[i]=pref*df_interpolators[i]->get_fdf(pp);
  }
}

void Field::gatherField( PlumedCommunicator& comm ){
  comm.Sum( &grid_buffer[0],grid_buffer.size() );
}

void Field::set_tables(){
  std::vector<Value> vv(npoints);
  for(unsigned i=0;i<npoints;++i){
     vv[i].set( grid_buffer[i*nper] );
     vv[i].resizeDerivatives( ndx );
     for(unsigned j=0;j<ndx;++j) vv[i].addDerivative( j, grid_buffer[ i*nper + j + 1 ] );
  }
  f_interpolator->set_table( vv );

  unsigned kk;
  for(unsigned n=0;n<df_interpolators.size();++n){
     for(unsigned i=0;i<npoints;++i){
        vv[i].clearDerivatives();
        kk=i*nper + (n+1)*(1+ndx); vv[i].set( grid_buffer[kk] ); kk++;
        for(unsigned j=0;j<ndx;++j){ vv[i].addDerivative( j, grid_buffer[kk] ); kk++; }
     }
     df_interpolators[n]->set_table( vv );
  } 
}

bool Field::applyForces( std::vector<double>& outforces ) const {
  plumed_assert( outforces.size()==forces.size() );
  for(unsigned i=0;i<forces.size();++i) outforces[i]=forces[i];
  return wasforced;
}

void Field::addForces( std::vector<double>& inforces ){
  plumed_assert( inforces.size()==forces.size() );
  wasforced=true;
  for(unsigned i=0;i<forces.size();++i) forces[i]=inforces[i];
}


void Field::error( const std::string & msg){
  errormsg=msg; 
}

std::string Field::errorMessage() const {
  return errormsg;
}

bool Field::check() const {
  if( errormsg.size()!=0 ) return false;
  return true;
}

void Field::printKeywords( Log& log ){
  Keywords field_keys;
  field_keys.add("compulsory","NPSPLINE","the number of points in each direction at which to calculate the value of the field");
  field_keys.add("compulsory","MIN","the minimum value at which the value of the field should be calculated in each direction");
  field_keys.add("compulsory","MAX","the maximum value at which the value of the field should be calculated in each direction");
  field_keys.add("compulsory","SIGMA","the value of the sigma parameter in the field");
  field_keys.print( log );
}
