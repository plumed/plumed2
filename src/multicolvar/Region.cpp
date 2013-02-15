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

#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "vesselbase/VesselRegister.h"
#include "vesselbase/FunctionVessel.h"
#include "ActionVolume.h"
#include "MultiColvar.h"
#include "tools/HistogramBead.h"

namespace PLMD {
namespace multicolvar {

//+PLUMEDOC INTERNAL region 
/*

Imagine we have a collection of collective variables that can all be assigned to a particular point in three 
dimensional space. For example, we could have the values of the coordination numbers for all the atoms in the
system. Because each CV value can be assigned to a particular point in space we can calculate the average
value of the cv in a particular part of the box using:

\f[
\overline{s}_{\tau} = \frac{ \sum_i s_i w(x_i,y_i,z_i) }{ \sum_i w(x_i,y_i,z_i) }  
\f]  

where the sum is over the collective variables, \f$s_i\f$, each of which can be thought to be at \f$ (x_i,y_i,z_i)\f$.
The function \f$ w(x_i,y_i,z_i) \f$ measures whether or not the system is in the subregion of interest. It
is equal to:

\f[
w(x_i,y_i,z_i) = \int_\tau \textrm{d}x\textrm{d}y\textrm{d}z K\left( \frac{x - x_i}{\sigma} \right)K\left( \frac{y - y_i}{\sigma} \right)K\left( \frac{z - z_i}{\sigma} \right) 
\f]

where \f$K\f$ is one of the kernel functions described on \ref histogrambead and \f$\sigma\f$ is a bandwidth parameter.
By default a Gaussian kernel is used by you can change the kernel using the KERNEL keyword.  
The value of \f$\sigma\f$ is then specified using the keyword SIGMA.
To define the region to integrate over, \f$\tau\f$, you use a separate action and reference it here
using its label and the VOLUME keyword.  The following actions can be used:  

<table align=center frame=void width=95%% cellpadding=5%%>
<tr> <td width=5%> \subpage  SUBCELL </td> <td>
Calculate the average value of the cv in a user-specified fraction of the box 
</td> </tr>
</table>

If you use the name of the action alone the quantity is calculated for the atoms inside the region.
When the name of the action is preceded by an ! then the quantity is calculated for all the atoms that are
not inside the region of specified by the using volume.
When REGION is used with the \ref DENSITY action the number of atoms in the specified region is calculated  

\par Examples

The following commands tell plumed to calculate the average coordination number for the atoms
that have x (in fractional coordinates) less than 0.5. The final value will be labeled r1_av.
\verbatim
SUBCELL XLOWER=0.0 XUPPER=0.5 LABEL=r1
COORDINATIONNUMBER SPECIES=1-100 R_0=1.0 REGION={SIGMA=0.1 VOLUME=r1}
\endverbatim

The following commands tell plumed to calculate the number of atoms in an ion chanel in a protein.
The extent of the chanel is calculated from the positions of atoms 1 4 5 10. The final value will be labeled r2_num.
\verbatim
CAVITY ATOMS=1,4,5,10 LABEL=r2
DENSITY SPECIES=20-500 REGION={SIGMA=0.1 VOLUME=r2}
\endverbatim

*/
//+ENDPLUMEDOC


class Region : public vesselbase::FunctionVessel {
private:
  std::string volname;
  double sigma;
  bool isDensity;
  Vector wdf;
  ActionVolume* myvol;
  MultiColvar* mycolv;
  bool not_in;
  HistogramBead bead;
public:
  static void registerKeywords( Keywords& keys );
  static void reserveKeyword( Keywords& keys );
  Region( const vesselbase::VesselOptions& da );
  unsigned getNumberOfTerms(){ return 3; }
  std::string function_description();
  bool calculate();
  void finish();
};

PLUMED_REGISTER_VESSEL(Region,"REGION")

void Region::registerKeywords( Keywords& keys ){
  FunctionVessel::registerKeywords( keys );
  keys.add("compulsory","VOLUME","the label for the action that describes the volume of interest");
  keys.add("compulsory","SIGMA","the width of the function to be used for kernel density estimation");
  keys.add("compulsory","KERNEL","gaussian","the type of kernel function to be used");
}

void Region::reserveKeyword( Keywords& keys ){
  keys.reserve("numbered","REGION","calculate the average value of the CV within a portion of the box. " 
                                    "For more details on how this quantity is calculated see \\ref region.",true);
}

Region::Region( const vesselbase::VesselOptions& da ) :
FunctionVessel(da),
not_in(false)
{
  mycolv=dynamic_cast<MultiColvar*>( getAction() );
  plumed_massert( mycolv, "REGION can only be used with MultiColvars");

  isDensity=mycolv->isDensity();

  parse("VOLUME",volname);
  if( volname.substr(0,1)=="!" ){ volname=volname.substr(1); not_in=true; }
  myvol = (mycolv->plumed).getActionSet().selectWithLabel<ActionVolume*>(volname);
  if(!myvol){ error( "in REGION " + volname + " is not a valid volume element"); return; }
  mycolv->centralAtomDerivativesAreInFractional=myvol->derivativesOfFractionalCoordinates();
  mycolv->addDependency( myvol );

  parse("SIGMA",sigma); myvol->setSigma( sigma );

  std::string kerneltype; parse("KERNEL",kerneltype);
  bead.isPeriodic( 0, 1.0 ); bead.setKernelType( kerneltype );
}

std::string Region::function_description(){
  std::string str_sigma; Tools::convert( sigma, str_sigma );
  std::string notin=""; if(not_in) notin=" not ";
  if( isDensity ) return "the number of atoms" + notin + "in " + volname + ". Sigma is equal to " + str_sigma;
  return "the average value of the cv" + notin + "in " + volname + ". Sigma is equal to " + str_sigma; 
}

bool Region::calculate(){
  Vector catom_pos=mycolv->retrieveCentralAtomPos();

  double weight; Vector wdf;
  weight=myvol->calculateNumberInside( catom_pos, bead, wdf );
  if( not_in ){ weight = 1.0 - weight; wdf *= -1.; }
  mycolv->setElementValue( 2, weight );  
 
  unsigned current=mycolv->current, ider=2*mycolv->getNumberOfDerivatives();
  for(unsigned i=0;i<mycolv->colvar_atoms[current].getNumberActive();++i){
     mycolv->addElementDerivative( ider, mycolv->getCentralAtomDerivative( i, 0, wdf ) ); ider++;
     mycolv->addElementDerivative( ider, mycolv->getCentralAtomDerivative( i, 1, wdf ) ); ider++;
     mycolv->addElementDerivative( ider, mycolv->getCentralAtomDerivative( i, 2, wdf ) ); ider++;
  }

  double ww=mycolv->getElementValue(1);  // Get the standard weight
  bool addval=addValue( 1, ww*weight );
  if(addval){
     double colvar=mycolv->getElementValue(0);
     addValue(0, colvar*weight*ww );
     getAction()->chainRuleForElementDerivatives( 1, 2, ww, this );
     if(diffweight) getAction()->chainRuleForElementDerivatives( 1, 1, weight, this );
     if(!isDensity){
         getAction()->chainRuleForElementDerivatives( 0, 0, ww*weight, this );
         getAction()->chainRuleForElementDerivatives( 0, 2, ww*colvar, this );
         if(diffweight) getAction()->chainRuleForElementDerivatives( 0, 1, colvar*weight, this );
     }
  } 
  return addval; 
}

void Region::finish(){
  if(isDensity){
     setOutputValue( getFinalValue(1) );
     std::vector<double> df(3); df[0]=0.0; df[1]=1.0; df[2]=0.0;
     mergeFinalDerivatives( df );
  } else {
     double denom=getFinalValue(1);
     setOutputValue( getFinalValue(0) / denom );
     std::vector<double> df(3); 
     df[0] = 1 / denom; 
     df[1] = - getFinalValue(0) / (denom*denom);
     df[2] = 0.0;
     mergeFinalDerivatives( df );
  }
}

}
}
