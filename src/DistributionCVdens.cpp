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

#include "FunctionVessel.h"
#include "MultiColvar.h"
#include "HistogramBead.h"

namespace PLMD {

//+PLUMEDOC INTERNAL subcell
/*
Imagine we have a collection of collective variables that can all be assigned to a particular point in three 
dimensional space. For example, we could have the values of the coordination numbers for all the atoms in the
system. Because each CV value can be assigned to a particular point in space we can calculate the average
value of the cv in a particular part of the box using:

\f[
\overline{s}_{\tau} = \frac{ \sum_i s_i w(x_i)w(y_i)w(z_i) }{ \sum_i w(x_i)w(y_i)w(z_i) }  
\f]  

where the sum is over the collective variables, \f$s_i\f$, each of which can be thought to be at \f$ (x_i,y_i,z_i)\f$.
Because we want to calculate the average cv in the subregion \f$\tau\f$ we introduce the three \f$w()\f$ to tell us
whether or not \f$i\f$ is inside \f$\tau\f$. These functions are given by:

\f[
w(x_i) = \int_{a_x}^{b_x} \textrm{d}x K\left( \frac{x - x_i}{w_x} \right)
\f]

where \f$K\f$ is one of the kernel functions described on \ref histogrambead.

All the input to calculate these quantities is provided through a single keyword that will have the following form

KEYWORD={TYPE XLOWER=\f$a_x\f$ XUPPER=\f$b_x\f$ XSMEAR=\f$\frac{w_x}{b_x-a_x}\f$ YLOWER=\f$a_y\f$ YUPPER=\f$b_y\f$ YSMEAR=\f$\frac{w_y}{b_y-a_y}\f$ ZLOWER=\f$a_z\f$ ZUPPER=\f$b_z\f$ ZSMEAR=\f$\frac{w_z}{b_z-a_z}\f$}

The \f$a_x\f$, \f$b_x\f$, \f$a_y\f$ ... parameters are specified as a fraction of the box length.  In directions for which upper and lower bounds are not
specified the \f$a\f$ and \f$b\f$ values are assumed equal to 0 and 1. If SMEAR values are not specified \f$\frac{w_x}{b_x-a_x}=0.5\f$.    
*/
//+ENDPLUMEDOC


class cvdens : public NormedSumVessel {
private:
  bool isDensity;
  Value tmpval, tmpweight;
  MultiColvar* mycolv;
  std::vector<Value> catom_pos;
  std::vector<unsigned> dir;
  std::vector<HistogramBead> beads;
  std::string getLabel();
  void calculateDensity( Value& weight );
public:
  static void reserveKeyword( Keywords& keys );
  cvdens( const VesselOptions& da );
  void getWeight( const unsigned& i, Value& weight );
  void compute( const unsigned& i, const unsigned& j, Value& theval );
  void printKeywords();
};

PLUMED_REGISTER_VESSEL(cvdens,"SUBCELL")

void cvdens::reserveKeyword( Keywords& keys ){
  keys.reserve("numbered","SUBCELL","calculate the average value of the CV within a portion of the box and store it in a value called subcell. " 
                                    "For more details on how this quantity is calculated see \\ref subcell.");
}

cvdens::cvdens( const VesselOptions& da ) :
NormedSumVessel(da),
catom_pos(3)
{
  mycolv=dynamic_cast<MultiColvar*>( getAction() );
  plumed_massert( mycolv, "cvdens can only be used with MultiColvars");

  mycolv->useCentralAtom();
  isDensity=mycolv->isDensity();
  if(!isDensity) useNorm(); 

  std::string errors;
  HistogramBead xbin; xbin.set(da.parameters, "X", errors); xbin.isPeriodic( 0, 1.0 ); 
  if ( xbin.hasBeenSet() ){ beads.push_back(xbin); dir.push_back(0); }
  HistogramBead ybin; ybin.set(da.parameters, "Y", errors); ybin.isPeriodic( 0., 1.0 );
  if ( ybin.hasBeenSet() ){ beads.push_back(ybin); dir.push_back(1); }
  HistogramBead zbin; zbin.set(da.parameters, "Z", errors); zbin.isPeriodic( 0, 1.0 );
  if ( zbin.hasBeenSet() ){ beads.push_back(zbin); dir.push_back(2); }  
  if( beads.size()==0 ) error("no subcell has been specified");

  addOutput( getLabel() );
  if(isDensity) log.printf("  value %s.%s contains density of atoms in box with ",(getAction()->getLabel()).c_str(),getLabel().c_str());
  else log.printf("  value %s.%s contains average value of cv for ",(getAction()->getLabel()).c_str(),getLabel().c_str());
  for(unsigned i=0;i<dir.size();++i){
     if( dir[i]==0 ) log.printf("x %s",(beads[i].description()).c_str());
     if( dir[i]==1 ) log.printf("y %s",(beads[i].description()).c_str());
     if( dir[i]==2 ) log.printf("z %s",(beads[i].description()).c_str());
     if( dir.size()>1 && i!=(dir.size()-1) ) log.printf(", "); 
  }
  log.printf("\n");
}

void cvdens::printKeywords(){
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

void cvdens::calculateDensity( Value& weight ){
  mycolv->retrieveCentralAtomPos( catom_pos ); double f, df;
  for(unsigned i=0;i<dir.size();++i){
     f=beads[i].calculate( catom_pos[ dir[i] ].get(), df );
     catom_pos[ dir[i] ].chainRule(df); catom_pos[ dir[i] ].set(f);
     if(i==0){ 
       copy( catom_pos[ dir[i] ], tmpweight ); 
     } else{ 
       product( tmpweight, catom_pos[ dir[i] ], weight );
       copy( weight, tmpweight );
     }
  }
  copy( tmpweight, weight );
}

void cvdens::getWeight( const unsigned& i, Value& weight ){
  calculateDensity( weight );
}

void cvdens::compute( const unsigned& i, const unsigned& j, Value& theval ){
  plumed_assert( j==0 );
  if(isDensity){
     calculateDensity( theval );
  } else {
     mycolv->retreiveLastCalculatedValue( tmpval );
     product( tmpval, tmpweight, theval ); 
  }
}

}
