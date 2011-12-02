#include "ColvarModifier.h"
#include "ColvarModifierFunctions.h"

namespace PLMD {

//+PLUMEDOC MODIFIER HISTOGRAM 
/**

Calculates the number of colvars that are within a certain range.  To make this quantity continuous it is calculated using:

\f[
S = \sum_i \int_a^b G( s_i, \sigma )
\f]

where \f$G( s_i, \sigma )\f$ is a normalized Gaussian function of width \f$\sigma\f$ centered at the value of the colvar \f$s_i\f$.  \f$a\f$ and \f$b\f$ are
the lower and upper bounds of the range of interest respectively.  The values of \f$a\f$ and \f$b\f$ can be specified using the RANGE keyword and the final
value of the function is stored as <label>.between<\f$a\f$>&<\f$b\f$>.

You can also use multiple ranges with this histogram keyword.  This can be done in one of two ways either you can use multiple instances of the RANGE keyword - i.e.
RANGE1, RANGE2, .... or, if you can specify that there should be a number of equally sized bins between some upper and lower bounds by using RANGE in combination with
the NBINS keyword.

*/
//+ENDPLUMEDOC


ColvarModifierHistogram::ColvarModifierHistogram(ColvarWithModifiers* act):
ColvarModifier(act)
{
  
  double sm=0.5; parse("SMEAR", true, sm);

  std::pair<double,double> range;
  std::vector<std::pair<double,double> > ranges;
  std::string rang="none"; parse("RANGE", true, rang);
  if( rang=="none"){
     std::string newrang, num; 
     for(int i=1;; ++i ){
        newrang="none"; Tools::convert(i,num);
        parse("RANGE" + num, true, newrang);
        if( newrang=="none") break;
        interpretRange( newrang, range );
        ranges.push_back( range );
     }
  } else {
     int nbins=1; parse("NBINS",true, nbins); 
     if( nbins<0 ) error("the number of bins is negative");
     interpretRange( rang, range );
     double delr = ( range.second - range.first ) / static_cast<double>(nbins);
     for(unsigned i=0;i<nbins;++i){
         ranges.push_back( std::pair<double,double>( range.first + i*delr, range.first + (i+1)*delr ) );
     }
  }
  if( ranges.size()==0 ) error("You must define a range/some ranges to use the histogram keyword");

  HistogramBead histo; std::string lb, ub, smear; 
  Tools::convert( sm, smear);
  for(unsigned i=0;i<ranges.size();++i){ 
      histo.set( ranges[i].first, ranges[i].second, sm );
      hf.push_back( histo ); totals.push_back(0);
      Tools::convert( ranges[i].first, lb ); 
      Tools::convert( ranges[i].second, ub );
      report("  calculating number of values between " + lb + " and " + ub + ".  Smear parameter equals " + smear ); 
      addValue("between" + lb + "&" + ub);    
  }
}

void ColvarModifierHistogram::finishCalculation(){
  for (unsigned j=0; j<hf.size(); ++j){ setValue( j, totals[j], 1.0 ); totals[j]=0.0; }
}

double ColvarModifierHistogram::differential( const unsigned& ival, const double& value ){
  assert( ival<hf.size() ); double df;
  totals[ival]+=hf[ival].calculate( value, df ); 
  return df;
} 

}
