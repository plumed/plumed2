#include "DistributionFunctions.h"

namespace PLMD {

std::string gradient::documentation(){
  std::ostringstream ostr;
  ostr<<"Gradient collective coordinates can be used to drive the formation of interfaces between phases. ";
  ostr<<"The gradient of a CV along an axis is calculated using \\f$ g = \\sum_{i=1}^{N} ( n_{i_1} - n_{i} )^2 \\f$ ";
  ostr<<"where \\f$N\\f$ is the number of bins you devide the axis into and \\f$n_i\\f$ is the average value of the CV in the \\f$i\\f$th bin. ";
  ostr<<"The average value of the CV in a bin is calculated using \\f$ n_i = \\frac{\\sum_{j=1}^M s_j w_i(x_j)}{\\sum_{j=1}^M w_i(x_j)}\\f$. ";
  ostr<<"The sum in this expression runs over all the cvs you are calculating as part of the MultiColvar and the value of \\f$ w_i(x_j) \\f$ is ";
  ostr<<"calculated using \\f$ w_i(x_j) = \\int_a^b \\frac{1}{\\sqrt{2\\pi}\\sigma} \\exp\\left( -\\frac{(x-x_j)^2}{2\\sigma^2} \\right) \\textrm{d}x \\f$ ";
  ostr<<"You can either specify that you would like to calculate the total gradient \\f$ g_x + g_y + g_z \\f$ by using ";
  ostr<<"(XBINS=\\f$N_x\\f$ YBINS=\\f$N_y\\f$ ZBINS=\\f$N_z\\f$ XSMEAR=\\f$\\sigma_x N_x\\f$ YSMEAR=\\f$\\sigma_y N_y\\f$ ZSMEAR=\\f$\\sigma_z N_z\\f$) ";
  ostr<<"Alternatively, specifying only XBINS will give you the gradient along the x-axis, while specifying XBINS and YBINS will give you \\f$ g_x + g_y \\f$ ";
  ostr<<"and so on. By default all SMEAR parameters are set equal to 0.5."; 
  return ostr.str();
}

gradient::gradient( const std::string& parameters ) :
DistributionFunction(parameters),
isDensity(false)
{
    unsigned xbins,ybins,zbins; HistogramBead tmpbead;
    std::vector<std::string> data=Tools::getWords(parameters);
    bool in_x=Tools::parse(data,"XBINS",xbins);
    bool in_y=Tools::parse(data,"YBINS",ybins);
    bool in_z=Tools::parse(data,"ZBINS",zbins);
    if( !in_x && !in_y && !in_z ) error("use XBINS/YBINS or ZBINS otherwise I do nothing");
    if(in_x){
       double smear=0.5; Tools::parse(data,"XSMEAR",smear);
       bounds.push_back( beads.size() );
       for(unsigned i=0;i<xbins;++i){
           tmpbead.set( (1.0/xbins)*i, (1.0/xbins)*(i+1), smear );
           beads.push_back( tmpbead );
           addAccumulator( true ); addAccumulator( true );
       }
       dir.push_back(0);
    }
    if(in_y){
       double smear=0.5; Tools::parse(data,"YSMEAR",smear);
       bounds.push_back( beads.size() ); 
       for(unsigned i=0;i<ybins;++i){
           tmpbead.set( (1.0/ybins)*i, (1.0/ybins)*(i+1), smear ); 
           beads.push_back( tmpbead );
           addAccumulator( true ); addAccumulator( true );
       } 
       dir.push_back(1);
    }
    if(in_z){
       double smear=0.5; Tools::parse(data,"ZSMEAR",smear);
       bounds.push_back( beads.size() );
       for(unsigned i=0;i<zbins;++i){
           tmpbead.set( (1.0/zbins)*i, (1.0/zbins)*(i+1), smear ); 
           beads.push_back( tmpbead );
           addAccumulator( true ); addAccumulator( true );
       } 
       dir.push_back(2);
    }
    bounds.push_back( beads.size() );
    final_bin.resize( beads.size() );
    bool idens=Tools::parseFlag(data,"density",isDensity);
    if( ! data.empty() ){
        std::string msg="found the following rogue keywords in switching function input : ";
        for(unsigned i=0;i<data.size();++i) msg = msg + data[i] + " "; 
        error(msg);
    } 
}

std::string gradient::message(){
  std::ostringstream ostr;
  if( dir.size()==3 ){
      ostr<<"gradient of the average cv value";
  } else if (dir.size()==2) {
      ostr<<"gradient in the ";
      if( dir[0]==0 && dir[1]==1 ) ostr<<"x and y directions";
      if( dir[0]==0 && dir[1]==2 ) ostr<<"x and z directions";
      if( dir[0]==1 && dir[1]==2 ) ostr<<"y and z directions"; 
  } else if (dir.size()==1) {
      ostr<<"gradient in the ";
      if( dir[0]==0 ) ostr<<"x direction";
      if( dir[0]==1 ) ostr<<"y direction";
      if( dir[0]==2 ) ostr<<"z direction"; 
  } else plumed_assert(0);
  return ostr.str();
}

void gradient::printKeywords( Log& log ){
  Keywords dkeys;
  dkeys.add("optional","XBINS","the number of bins to use in the x direction");
  dkeys.add("optional","XSMEAR","(default=0.5) the amount to smear the positions by in the x direction");
  dkeys.add("optional","YBINS","the number of bins to use in the y direction");
  dkeys.add("optional","YSMEAR","(default=0.5) the amount to smear the positions by in the y direction");
  dkeys.add("optional","ZBINS","the number of bins to use in the z direction");
  dkeys.add("optional","ZSMEAR","(default=0.5) the amount to smear the positions by in the z direction");
  dkeys.print(log);
}

std::string gradient::getLabel(){
  std::string lab;
  if( dir.size()==3 ){
      lab="gradient";
  } else if (dir.size()==2) {
      if( dir[0]==0 && dir[1]==1 ) lab="xygradient";
      if( dir[0]==0 && dir[1]==2 ) lab="xzgradient";
      if( dir[0]==1 && dir[1]==2 ) lab="yzgradient";
  } else if (dir.size()==1) {
      if( dir[0]==0 ) lab="xgradient";
      if( dir[0]==1 ) lab="ygradient";
      if( dir[0]==2 ) lab="zgradient";
  } else plumed_assert(0);
  return lab;
}

void gradient::calculate( Value* value_in, std::vector<Value>& aux ){
  plumed_assert( aux.size()==3 );

  unsigned nn=0;  double f, df;
  for(unsigned i=0;i<dir.size();++i){
     for(unsigned j=bounds[i];j<bounds[i+1];++j){
        copyValue( nn, value_in ); 
        copyValue( nn+1, &aux[ dir[i] ] );
        f=beads[j].calculate( aux[ dir[i] ].get(), df );
        getPntrToValue( nn+1 )->chainRule( df ); getPntrToValue( nn+1 )->set( f );
        multiplyValue( nn, getPntrToValue( nn+1 ) ); 
        nn+=2;
     }
  }
  plumed_assert( getNumberOfAccumulators()==nn );
}

void gradient::finish( Value* value_out ){
  if( !isDensity ){ 
     unsigned nn=0;
     for(unsigned i=0;i<dir.size();++i){
        for(unsigned j=bounds[i];j<bounds[i+1];++j){
            unsigned nder=getPntrToAccumulator(nn)->getNumberOfDerivatives();
            if( final_bin[j].getNumberOfDerivatives()!=nder ) final_bin[j].resizeDerivatives( nder );
            quotient( getPntrToAccumulator(nn), getPntrToAccumulator(nn+1), &final_bin[j] );
            nn+=2;
        }
     }
     plumed_assert( getNumberOfAccumulators()==nn );
  } else{ 
     unsigned nn=0;
     for(unsigned i=0;i<dir.size();++i){
        for(unsigned j=bounds[i];j<bounds[i+1];++j){
            unsigned nder=getPntrToAccumulator(nn)->getNumberOfDerivatives();
            if( final_bin[j].getNumberOfDerivatives()!=nder ) final_bin[j].resizeDerivatives( nder );
            extractDerivatives( nn+1, &final_bin[j] ); final_bin[j].set( getPntrToAccumulator(nn+1)->get() );
            nn+=2;
        }
     }
     plumed_assert( getNumberOfAccumulators()==nn ); 
  }
  double tmp,tval=0; unsigned jj;
  for(unsigned i=0;i<dir.size();++i){
      for(unsigned j=bounds[i]+1;j<bounds[i+1];++j){
         tmp=final_bin[j-1].get() - final_bin[j].get();
         for(unsigned k=0;k<final_bin[j].getNumberOfDerivatives();++k){
             value_out->addDerivative( k, 2*tmp*(final_bin[j-1].getDerivative(k) - final_bin[j].getDerivative(k) ) );
         }
         tval+=tmp*tmp; 
      }
      jj=bounds[i+1]-1; 
      tmp=final_bin[jj].get() - final_bin[ bounds[i] ].get();
      for(unsigned k=0;k<final_bin[jj].getNumberOfDerivatives();++k){
          value_out->addDerivative( k, 2*tmp*(final_bin[jj].getDerivative(k) - final_bin[bounds[i]].getDerivative(k) ) );
      }
      tval+=tmp*tmp;
  }
  value_out->set(tval);
}

}
