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

#include "VesselValueAccess.h"
#include "MultiColvar.h"
#include "HistogramBead.h"

namespace PLMD {

class gradient : public VesselAccumulator {
private:
  bool isDensity;
  MultiColvar* mycolv;
  std::vector<unsigned> dir, bounds;
  Value myval, tmpval, tmpval2, tmpweight, tmpweight2;
  std::vector<Value> catom_pos,final_bin;
  std::vector<HistogramBead> beads;
  std::string getLabel();
public:
  static void reserveKeyword( Keywords& keys );
  gradient( const VesselOptions& da );
  bool calculate( const unsigned& icv, const double& tolerance );
  void finish( const double& tolerance );
  void printKeywords();
};

PLUMED_REGISTER_VESSEL(gradient,"GRADIENT")

void gradient::reserveKeyword( Keywords& keys ){
  keys.reserve("numbered","GRADIENT","calcualte the gradient of a CV across the box and store it in value called gradient. "
  "Gradient collective coordinates can be used to drive the formation of interfaces between phases. "
  "The gradient of a CV along an axis is calculated using \\f$ g = \\sum_{i=1}^{N} ( n_{i_1} - n_{i} )^2 \\f$ "
  "where \\f$N\\f$ is the number of bins you devide the axis into and \\f$n_i\\f$ is the average value of the CV in the \\f$i\\f$th bin. "
  "The average value of the CV in a bin is calculated using \\f$ n_i = \\frac{\\sum_{j=1}^M s_j w_i(x_j)}{\\sum_{j=1}^M w_i(x_j)}\\f$. "
  "The sum in this expression runs over all the cvs you are calculating as part of the MultiColvar and the value of \\f$ w_i(x_j) \\f$ is "
  "calculated using \\f$ w_i(x_j) = \\int_a^b \\frac{1}{\\sqrt{2\\pi}\\sigma} \\exp\\left( -\\frac{(x-x_j)^2}{2\\sigma^2} \\right) \\textrm{d}x \\f$ "
  "You can either specify that you would like to calculate the total gradient \\f$ g_x + g_y + g_z \\f$ by using "
  "(XBINS=\\f$N_x\\f$ YBINS=\\f$N_y\\f$ ZBINS=\\f$N_z\\f$ XSMEAR=\\f$\\sigma_x N_x\\f$ YSMEAR=\\f$\\sigma_y N_y\\f$ ZSMEAR=\\f$\\sigma_z N_z\\f$) "
  "Alternatively, specifying only XBINS will give you the gradient along the x-axis, while specifying XBINS and YBINS will give you \\f$ g_x + g_y \\f$ "
  "and so on. By default all SMEAR parameters are set equal to 0.5."); 
}

gradient::gradient( const VesselOptions& da ) :
VesselAccumulator(da),
isDensity(false),
catom_pos(3)
{
    mycolv=dynamic_cast<MultiColvar*>( getAction() );
    plumed_massert( mycolv, "cvdens can only be used with MultiColvars");
  
    mycolv->useCentralAtom();
    isDensity=mycolv->isDensity();

    unsigned xbins,ybins,zbins; HistogramBead tmpbead;
    std::vector<std::string> data=Tools::getWords(da.parameters);
    bool in_x=Tools::parse(data,"XBINS",xbins);
    bool in_y=Tools::parse(data,"YBINS",ybins);
    bool in_z=Tools::parse(data,"ZBINS",zbins);
    if( !in_x && !in_y && !in_z ) error("use XBINS/YBINS or ZBINS otherwise I do nothing");
    if(in_x){
       double smear=0.5; Tools::parse(data,"XSMEAR",smear);
       bounds.push_back( beads.size() );
       for(unsigned i=0;i<xbins;++i){
           tmpbead.set( (1.0/xbins)*i, (1.0/xbins)*(i+1), smear );
           tmpbead.isPeriodic( 0., 1.0 );
           beads.push_back( tmpbead );
           addBufferedValue(); addBufferedValue();
       }
       dir.push_back(0);
    }
    if(in_y){
       double smear=0.5; Tools::parse(data,"YSMEAR",smear);
       bounds.push_back( beads.size() ); 
       for(unsigned i=0;i<ybins;++i){
           tmpbead.set( (1.0/ybins)*i, (1.0/ybins)*(i+1), smear ); 
           tmpbead.isPeriodic( 0., 1.0 );
           beads.push_back( tmpbead );
           addBufferedValue(); addBufferedValue();
       } 
       dir.push_back(1);
    }
    if(in_z){
       double smear=0.5; Tools::parse(data,"ZSMEAR",smear);
       bounds.push_back( beads.size() );
       for(unsigned i=0;i<zbins;++i){
           tmpbead.set( (1.0/zbins)*i, (1.0/zbins)*(i+1), smear ); 
           tmpbead.isPeriodic( 0., 1.0 );
           beads.push_back( tmpbead );
           addBufferedValue(); addBufferedValue();
       } 
       dir.push_back(2);
    }
    bounds.push_back( beads.size() );
    final_bin.resize( beads.size() );

    if( !data.empty() ){
        std::string msg="found the following rogue keywords in switching function input : ";
        for(unsigned i=0;i<data.size();++i) msg = msg + data[i] + " "; 
        error(msg);
    } 
 
   addOutput(getLabel()); 
   log.printf("  value %s.%s contains the ",(getAction()->getLabel()).c_str(),getLabel().c_str());
   if( dir.size()==3 ){
       log.printf("gradient of the average cv value\n");;
   } else if (dir.size()==2) {
       log.printf("gradient in the ");
       if( dir[0]==0 && dir[1]==1 ) log.printf("x and y directions\n");
       if( dir[0]==0 && dir[1]==2 ) log.printf("x and z directions\n");
       if( dir[0]==1 && dir[1]==2 ) log.printf("y and z directions\n");
   } else if (dir.size()==1) {
       log.printf("gradient in the ");
       if( dir[0]==0 ) log.printf("x direction\n");
       if( dir[0]==1 ) log.printf("y direction\n");
       if( dir[0]==2 ) log.printf("z direction\n"); 
   } else plumed_assert(0);
}

void gradient::printKeywords(){
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

bool gradient::calculate( const unsigned& icv, const double& tolerance ){
  bool keep=false; double f, df; 
  mycolv->retreiveLastCalculatedValue( myval );
  mycolv->retrieveCentralAtomPos( catom_pos );

  unsigned nn=0; double pos; 
  for(unsigned i=0;i<dir.size();++i){
     pos=catom_pos[ dir[i] ].get();
     for(unsigned j=bounds[i];j<bounds[i+1];++j){
         f=beads[j].calculate( pos , df );
         if( f>tolerance ){
             keep=true;
             copy( catom_pos[ dir[i] ], tmpweight );
             tmpweight.chainRule(df); tmpweight.set(f);

             tmpweight2.set( tmpweight.get() );
             mycolv->mergeDerivatives( icv, tmpweight, 1.0, tmpweight2 );
             addValue( nn+1, tmpweight2 );
            
             product( myval, tmpweight, tmpval );
             if( fabs( tmpval.get() )>tolerance ){ 
               tmpval2.set( tmpval.get() );
               mycolv->mergeDerivatives( icv, tmpval, 1.0, tmpval2 );
               addValue( nn, tmpval2 ); 
             }
         }
         nn+=2;
     }
  }
  return keep;
}

void gradient::finish( const double& tolerance ){
  if( !isDensity ){ 
     unsigned nn=0;
     for(unsigned i=0;i<dir.size();++i){
        for(unsigned j=bounds[i];j<bounds[i+1];++j){
            getValue( nn, tmpval2 ); getValue( nn+1, tmpweight2 );
            quotient( tmpval2, tmpweight2, &final_bin[j] );
            nn+=2;
        }
     }
  } else{ 
     unsigned nn=0;
     for(unsigned i=0;i<dir.size();++i){
        for(unsigned j=bounds[i];j<bounds[i+1];++j){
            getValue( nn+1, tmpweight2 );
            copy( tmpweight2, final_bin[j] );
            nn+=2;
        }
     }
  }

  Value* value_out=getPntrToOutput(0);
  double tmp, tval=0; unsigned jj;
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
