#include "DistributionFunctions.h"
#include "Value.h"
#include "Keywords.h"

namespace PLMD {

std::string cvdens::documentation(){
  std::ostringstream ostr;
  ostr<<"To make this quantity continuous it is calculated using \\f$ a = \\frac{\\sum_{i=1}^N s_i w(x_i)w(y_i)w(z_i)}{\\sum_{i=1}^N w(x_i)w(y_i)w(z_i)}\\f$ ";
  ostr<<"where the sum is over the collective variables \\f$ s_i \\f$ which can be thought to be at \\f$ (x_i,y_i,z_i)\\f$. The three \\f$w(x)\\f$ functions ";
  ostr<<"specify the extent of the subregion in which we are calculating the average CV in each direction. To make these assignment functions continuous they ";
  ostr<<"are calculated using "<<HistogramBead::documentation(true)<<".  If no parameters for the \\f$ w(d) \\f$ function in a particular direction is specified ";
  ostr<<"then the subcell is assumed to incorporate the entirity of the box in that direction."; 
  return ostr.str();
}

cvdens::cvdens( const std::string& parameters ) :
DistributionFunction(parameters)
{
  std::string errors;
  HistogramBead xbin; xbin.set(parameters, "X", errors); 
  if ( xbin.hasBeenSet() ){ beads.push_back(xbin); dir.push_back(0); }
  HistogramBead ybin; ybin.set(parameters, "Y", errors);
  if ( ybin.hasBeenSet() ){ beads.push_back(ybin); dir.push_back(1); }
  HistogramBead zbin; zbin.set(parameters, "Z", errors);
  if ( zbin.hasBeenSet() ){ beads.push_back(zbin); dir.push_back(2); }  
  if( beads.size()==0 ) error("no subcell has been specified");

  isDensity=(parameters.find("density")!=std::string::npos);
  addAccumulator( true );    // This holds the numerator - value of cv times "location in box" 
  addAccumulator( true );    // This holds the denominator - number of atoms in box
  addAccumulator( true );    // This holds a tempory value during the calculation 
}

std::string cvdens::message(){
  std::ostringstream ostr;
  ostr<<"average value of cv for ";
  for(unsigned i=0;i<dir.size();++i){
     if( dir[i]==0 ) ostr<<"x "<<beads[i].description();
     if( dir[i]==1 ) ostr<<"y "<<beads[i].description();
     if( dir[i]==2 ) ostr<<"z "<<beads[i].description();
     if( dir.size()>1 && i!=(dir.size()-1) ) ostr<<", "; 
  }
  return ostr.str();
}

void cvdens::printKeywords( Log& log ){
  Keywords dkeys;
  dkeys.add("optional","XLOWER","the lower boundary in x of the subcell in fractional coordinates.  If this keyword is absent then the lower boundary is placed at 0.0"); 
  dkeys.add("optional","XUPPER","the upper boundary in x of the subcell in fractional coordinates.  If this keyword is absent then the upper boundary is placed at 1.0");
  dkeys.add("optional","XSMEAR","(default=0.5) the width of the Gaussian that is used to smear the density in the x direction");
  dkeys.add("optional","YLOWER","the lower boundary in y of the subcell in fractional coordinates.  If this keyword is absent then the lower boundary is placed at 0.0");
  dkeys.add("optional","YUPPER","the upper boundary in y of the subcell in fractional coordinates.  If this keyword is absent then the upper boundary is placed at 1.0");
  dkeys.add("optional","YSMEAR","(default=0.5) the width of the Gaussian that is used to smear the density in the y direction");
  dkeys.add("optional","ZLOWER","the lower boundary in z of the subcell in fractional coordinates.  If this keyword is absent then the lower boundary is placed at 0.0");
  dkeys.add("optional","ZUPPER","the upper boundary in z of the subcell in fractional coordinates.  If this keyword is absent then the upper boundary is placed at 1.0");
  dkeys.add("optional","ZSMEAR","(default=0.5) the width of the Gaussian that is used to smear the density in the z direction"); 
  dkeys.print(log);
}

std::string cvdens::getLabel(){
  std::string lb,ub,lab;
  if(isDensity) lab = "densityFor";
  else lab = "averageFor";
  for(unsigned i=0;i<dir.size();++i){
     Tools::convert( beads[i].getlowb(),  lb ); 
     Tools::convert( beads[i].getbigb(), ub );
     if(dir[i]==0) lab=lab + "Xin" + lb +"&" + ub; 
     if(dir[i]==1) lab=lab + "Yin" + lb +"&" + ub;
     if(dir[i]==2) lab=lab + "Zin" + lb +"&" + ub;
  }
  return lab;
}

void cvdens::calculate( Value* value_in, std::vector<Value>& aux ){
  plumed_assert( aux.size()==3 );

  double f, df;
  copyValue( 0, value_in );
  for(unsigned i=0;i<dir.size();++i){
     copyValue( 2, &aux[ dir[i] ] ); 
     f=beads[i].calculate( aux[ dir[i] ].get(), df );
     getPntrToValue(2)->chainRule(df); getPntrToValue(2)->set(f);
     multiplyValue( 0, getPntrToValue(2) );
     if(i==0){ copyValue( 1, getPntrToValue(2) ); } 
     else{ multiplyValue( 1, getPntrToValue(2) ); }
     getPntrToValue(2)->clearDerivatives();
  }
}

void cvdens::finish( Value* value_out ){
  if( !isDensity ){ quotient( getPntrToAccumulator(0), getPntrToAccumulator(1), value_out ); }
  else{ extractDerivatives( 1, value_out ); value_out->set( getPntrToAccumulator(1)->get() ); }
}

}
