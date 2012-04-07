#include "ActionPilot.h"
#include "ActionRegister.h"
#include "ActionWithValue.h"
#include "ActionWithField.h"
#include "PlumedMain.h"
#include "Atoms.h"

namespace PLMD {

//+PLUMEDOC BIAS FIELD_METAD 
/**
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

class BiasFieldMetadynamics : public ActionWithField {
private:
  unsigned freq;
  double hw;
  double biasf;
  double temp;
  bool welltemp;
public:
  BiasFieldMetadynamics(const ActionOptions&);
  static void registerKeywords(Keywords& keys);
  void update();
};

PLUMED_REGISTER_ACTION(BiasFieldMetadynamics,"FIELD_METAD")

void BiasFieldMetadynamics::registerKeywords(Keywords& keys){
  ActionWithField::registerKeywords(keys);
  keys.add("compulsory","PACE","the frequency with which to add hills");
  keys.add("compulsory","HEIGHT","the heights of the hills");
  keys.add("optional","BIASFACTOR","use well tempered metadynamics and use this biasfactor.  Please note you must also specify temp");
  keys.add("optional","TEMP","the system temperature - this is only needed if you are doing well-tempered metadynamics");
}

BiasFieldMetadynamics::BiasFieldMetadynamics(const ActionOptions& ao):
Action(ao),
ActionWithField(ao),
freq(0),
hw(0),
biasf(1.0),
temp(0.0),
welltemp(false)
{
  parse("PACE",freq); 
  parse("HEIGHT",hw);
  parse("BIASFACTOR",biasf); 
  if( biasf<1.0 ) error("Bias factor has not been set properly it must be greater than 1");
  parse("TEMP",temp);
  if( biasf>1.0 && temp<0.0 ) error("You must set the temperature using TEMP when you do well tempered metadynamics");
  checkRead();
  if( biasf>1.0 ) welltemp=true;
}

void BiasFieldMetadynamics::update(){
  if( getStep()%freq==0 ){
     double this_ww;
     if(welltemp){
        this_ww = hw*exp(-getPntrToComponent("bias")->get()/(plumed.getAtoms().getKBoltzmann()*temp*(biasf-1.0)));
     } else {
        this_ww=hw;
     }  
     addFieldToBias( this_ww );
  }
}

}
