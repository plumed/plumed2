#include "DistributionFunctions.h"

namespace PLMD {

void cvdens::writeDocs( std::string& docs ){
  std::ostringstream ostr;
  ostr<<"\\par CV_DENSITY_X/CV_DENSITY_Y/CV_DENSITY_Z \n\n";
  docs=ostr.str();
}

cvdens::cvdens( const std::vector<std::string>& parameters ) :
DistributionFunction(parameters)
{
  Tools::convert(parameters[0],ax);
  Tools::convert(parameters[1],bx);
  Tools::convert(parameters[2],xsig);
  Tools::convert(parameters[3],ay);
  Tools::convert(parameters[4],by);
  Tools::convert(parameters[5],ysig);
  Tools::convert(parameters[6],az);
  Tools::convert(parameters[7],bz);
  Tools::convert(parameters[8],zsig);
  beads.resize(3);
  if( !(ax==0.0 && bx==1.0) ){ beads[0].set(ax,bx,xsig); dir.push_back(0); }
  if( !(ay==0.0 && by==1.0) ){ beads[1].set(ay,by,ysig); dir.push_back(1); } 
  if( !(az==0.0 && bz==1.0) ){ beads[2].set(az,bz,zsig); dir.push_back(2); }
  plumed_assert( dir.size()>0 && dir.size()<=3 );
  addAccumulator( true );    // This holds the numerator - value of cv times "location in box" 
  addAccumulator( true );    // This holds the denominator - number of atoms in box
  isDensity=( parameters[9]=="isdensity" );
}

std::string cvdens::message(){
  std::ostringstream ostr;
  ostr<<"average value of cv for ";
  for(unsigned i=0;i<dir.size();++i){
     if( dir[i]==0 ) ostr<<"x between "<<ax<<" and "<< bx<<" ";
     if( dir[i]==1 ) ostr<<"y between "<<ay<<" and "<< by<<" ";
     if( dir[i]==2 ) ostr<<"z between "<<az<<" and "<< bz<<" ";
  }
  return ostr.str();
}

void cvdens::calculate( Value* value_in, std::vector<Value>& aux ){
  plumed_assert( aux.size()==3 );

  double f, df;
  copyValue( 0, value_in );
  for(unsigned i=0;i<dir.size();++i){
     f=beads[ dir[i] ].calculate( aux[ dir[i] ].get(), df );
     aux[ dir[i] ].chainRule(df); aux[ dir[i] ].set(f);
     multiplyValue( 0, &aux[ dir[i] ] );
     if( i==0 ){ copyValue( 1, &aux[ dir[i] ] ); }
     else{ multiplyValue( 1, &aux[ dir[i] ] ); }
  }
}

void cvdens::finish( Value* value_out ){
  if( !isDensity ){ quotient( getPntrToAccumulator(0), getPntrToAccumulator(1), value_out ); }
  else{ extractDerivatives( 1, value_out ); value_out->set( getPntrToAccumulator(1)->get() ); }
}

}
