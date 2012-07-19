/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012 The PLUMED team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

   This file is part of PLUMED, version 2.0.

   PLUMED is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   PLUMED is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with PLUMED.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "ActionPilot.h"
#include "ActionRegister.h"
#include "ActionWithValue.h"
#include "FieldBias.h"
#include "PlumedMain.h"
#include "Atoms.h"

namespace PLMD {

//+PLUMEDOC BIAS FIELD_METAD 
/*
Field metadynamics is a method for enhancing sampling, which, much like conventional
metadynamcis, generates a bias based on the history of visited configurations.  However,
unlike metadynamics, for which the instantaneous state of the system is represented using
a vector of collective variables, the instantaneous state of the system is represented 
using a continueous probability distribution \f$\psi(X,z)\f$ that is calcualted based on
the instantaneous atomic positions.  

This means that the bias at any given time is calculated as an overlap integral \cite field-cvs namely:

\f[
V(X,t) = \int \textrm{d}z \psi(X(t),z) \sum_{t'=0}^t \psi(X(t'),z)
\f] 

\par Examples
The following input is performing field metadynamics using the histogram
of distances between the atoms in the specified group to describe the instantaneous
state fo the system
\verbatim
DISTANCES GROUP=1-7 FIELD=(MIN=0.5 MAX=3.0 NSPLINE=20 SIGMA=0.1) LABEL=f1
FIELD_METAD FIELD=f1 NGRID=400 STRIDE=1 PACE=10 HEIGHT=40.0 LABEL=m1
\endverbatim
(See also \ref DISTANCES)

*/
//+ENDPLUMEDOC

class BiasFieldMetadynamics : public FieldBias {
private:
  unsigned freq;
  double hw;
  double biasf;
  double temp;
  bool welltemp, gnuff;
  unsigned wgridstride;
  unsigned bnumber;
  std::string gridfname;
public:
  BiasFieldMetadynamics(const ActionOptions&);
  static void registerKeywords(Keywords& keys);
  void update();
};

PLUMED_REGISTER_ACTION(BiasFieldMetadynamics,"FIELD_METAD")

void BiasFieldMetadynamics::registerKeywords(Keywords& keys){
  FieldBias::registerKeywords(keys);
  keys.add("compulsory","PACE","the frequency with which to add hills");
  keys.add("compulsory","HEIGHT","the heights of the hills");
  keys.add("optional","BIASFACTOR","use well tempered metadynamics and use this biasfactor.  Please note you must also specify temp");
  keys.add("optional","TEMP","the system temperature - this is only needed if you are doing well-tempered metadynamics");
  keys.add("optional","GRID_WSTRIDE","the frequency with which to output the bias");
  keys.add("optional","GRID_WFILE","an output file to write the bias to");
  keys.addFlag("STORE_BIAS",false,"periodically store the bias - i.e. don't overwrite the bias files");
  keys.addFlag("GNUFF",false,"are we using gnuff");
}

BiasFieldMetadynamics::BiasFieldMetadynamics(const ActionOptions& ao):
Action(ao),
FieldBias(ao),
freq(0),
hw(0),
biasf(1.0),
temp(0.0),
welltemp(false),
gnuff(false),
wgridstride(0)
{
  parse("PACE",freq); 
  parse("HEIGHT",hw);
  log.printf("  adding hills with height %f every %d steps\n",hw,freq);
  parse("BIASFACTOR",biasf); 
  if( biasf<1.0 ) error("Bias factor has not been set properly it must be greater than 1");
  parse("TEMP",temp);
  if( biasf>1.0 && temp<0.0 ) error("You must set the temperature using TEMP when you do well tempered metadynamics");
  parseFlag("GNUFF",gnuff);
  if( biasf>1.0 && gnuff) log.printf("  doing well tempered metadynamics with bias factor %f and gnuff. System temperature is %f \n",biasf, temp );
  else if( biasf>1.0 ) log.printf("  doing well tempered metadynamics with bias factor %f. System temperature is %f \n",biasf, temp );

  parse("GRID_WSTRIDE",wgridstride );
  parse("GRID_WFILE",gridfname ); 
  bool bstore; parseFlag("STORE_BIAS",bstore);
  if( wgridstride!=0 && gridfname.length()==0 ) error("You must set the name of the output file, use GRID_WFILE");
  if( bstore ){
     bnumber=1;
     log.printf("  writing the bias every %d steps to numbered fiels named %s\n", wgridstride, gridfname.c_str() );
  } else {
     bnumber=0;
     log.printf("  writing the bias every %d steps to file %s\n", wgridstride, gridfname.c_str() );
  }
  checkRead();
  if( biasf>1.0 ) welltemp=true;
}

void BiasFieldMetadynamics::update(){
  if( getStep()%freq==0 ){
     double this_ww, deltaT, normali;
     
     Grid* bias=getPntrToBias(); normali=get_normalizer(); 
     std::vector<double> buffer( get_buffer() );
     if(welltemp){
        deltaT=plumed.getAtoms().getKBoltzmann()*temp*(biasf-1.0);
        this_ww = hw*exp(-getPntrToComponent("bias")->get()/deltaT);
        if( gnuff ){
           for(unsigned i=0;i<bias->getSize();++i){ bias->addValue( i, this_ww*exp(-bias->getValue(i)/deltaT)*buffer[i+2]/normali ); }
        } else {
           for(unsigned i=0;i<bias->getSize();++i){ bias->addValue( i, this_ww*buffer[i+2]/normali ); }
        }
     } else {
        for(unsigned i=0;i<bias->getSize();++i){ bias->addValue( i, hw*buffer[i+2]/normali ); }
     }  
  }

  if( wgridstride>0 && getStep()%wgridstride==0 ){
     FILE* gridfile;
     if( bnumber>0 ){
        std::string num, name;
        Tools::convert( bnumber, num );
        name = gridfname + "." + num;
        gridfile=fopen(name.c_str(),"w");
        bnumber++;
     } else {
        gridfile=fopen(gridfname.c_str(),"w"); 
     }
     getPntrToBias()->writeToFile( gridfile );
     fclose( gridfile );
  }
}

}
