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
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "core/ActionVolume.h"
#include "MultiColvar.h"
#include "tools/HistogramBead.h"

namespace PLMD {

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


class VesselCVDens : public NormedSumVessel {
private:
  bool isDensity;
  Value tmpval, tmpweight;
  ActionVolume* myvol;
  MultiColvar* mycolv;
  std::vector<Value> catom_pos;
  bool not_in;
  HistogramBead bead;
public:
  static void reserveKeyword( Keywords& keys );
  VesselCVDens( const VesselOptions& da );
  void getWeight( const unsigned& i, Value& weight );
  void compute( const unsigned& i, const unsigned& j, Value& theval );
};

PLUMED_REGISTER_VESSEL(VesselCVDens,"REGION")

void VesselCVDens::reserveKeyword( Keywords& keys ){
  keys.reserve("numbered","REGION","calculate the average value of the CV within a portion of the box and store it in a value called label_av. " 
                                    "For more details on how this quantity is calculated see \\ref region.");
}

VesselCVDens::VesselCVDens( const VesselOptions& da ) :
NormedSumVessel(da),
catom_pos(3),
not_in(false)
{
  mycolv=dynamic_cast<MultiColvar*>( getAction() );
  plumed_massert( mycolv, "REGION can only be used with MultiColvars");

  mycolv->useCentralAtom();
  isDensity=mycolv->isDensity();
  if(!isDensity) useNorm(); 

  std::vector<std::string> data=Tools::getWords(da.parameters);
  std::string name; bool found_element=Tools::parse(data,"VOLUME",name);
  if(!found_element){ error("No volume element specified in ccall to REGION"); return; }
  if( name.substr(0,1)=="!" ){ name=name.substr(1); not_in=true; }
  myvol = (mycolv->plumed).getActionSet().selectWithLabel<ActionVolume*>(name);
  if(!myvol){ error( "in REGION " + name + " is not a valid volume element"); return; }
  mycolv->addDependency( myvol ); 

  double sigma; bool found_sigma=Tools::parse(data,"SIGMA",sigma);
  if(!found_sigma){ error("No value for SIGMA specified in call to REGION"); return; }
  myvol->setSigma( sigma );

  std::string kerneltype; bool found_kernel=Tools::parse(data,"KERNEL",kerneltype);
  if(!found_kernel) kerneltype="gaussian";
  bead.isPeriodic( 0, 1.0 ); bead.setKernelType( kerneltype );

  if( isDensity ){
      addOutput( name + "_num" );
      log.printf("  value %s.%s_num contains number of ions in %s \n",(getAction()->getLabel()).c_str(),name.c_str(),name.c_str());
  } else {
      addOutput( name + "_av");
      log.printf("  value %s.%s_av contains average value of cv in %s \n",(getAction()->getLabel()).c_str(),name.c_str(),name.c_str());
  }
}

void VesselCVDens::getWeight( const unsigned& i, Value& weight ){
  mycolv->retrieveCentralAtomPos( catom_pos );
  plumed_assert( catom_pos.size()==3 ); unsigned nder=catom_pos[0].getNumberOfDerivatives();
  plumed_assert( nder==catom_pos[1].getNumberOfDerivatives() && nder==catom_pos[2].getNumberOfDerivatives() );
  if( tmpweight.getNumberOfDerivatives()!=nder ) tmpweight.resizeDerivatives( nder );
  tmpweight.clearDerivatives();
  myvol->calculateNumberInside( catom_pos, bead, tmpweight );
  copy( tmpweight, weight );
  if( not_in ){ weight.set( 1.0 - weight.get() ); weight.chainRule(-1.); }
}

void VesselCVDens::compute( const unsigned& i, const unsigned& j, Value& theval ){
  plumed_assert( j==0 );
  if(isDensity){
     mycolv->retrieveCentralAtomPos( catom_pos );
     plumed_assert( catom_pos.size()==3 ); unsigned nder=catom_pos[0].getNumberOfDerivatives();
     plumed_assert( nder==catom_pos[1].getNumberOfDerivatives() && nder==catom_pos[2].getNumberOfDerivatives() );
     if( theval.getNumberOfDerivatives()!=nder ) theval.resizeDerivatives( nder );
     myvol->calculateNumberInside( catom_pos, bead, theval );
     if( not_in ){ theval.set( 1.0 - theval.get() ); theval.chainRule(-1.); }
  } else {
     tmpval=mycolv->retreiveLastCalculatedValue();
     product( tmpval, tmpweight, theval ); 
  }
}

}
