/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2019 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed.org for more information.

   This file is part of plumed, version 2.

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
#include "VesselRegister.h"
#include "Vessel.h"
#include "StoreDataVessel.h"
#include "ActionWithVessel.h"

namespace PLMD {
namespace vesselbase {

// This is not the most efficient implementation
// The calculation of all the colvars is parallelized
// but the loops for calculating moments are not
// Feel free to reimplement this if you know how
class Moments : public Vessel {
private:
  unsigned mycomponent;
  StoreDataVessel* mystash;
  std::vector<unsigned> powers;
  std::vector<Value*> value_out;
public:
  static void registerKeywords( Keywords& keys );
  static void reserveKeyword( Keywords& keys );
  explicit Moments( const vesselbase::VesselOptions& da );
  std::string description();
  void resize();
  void calculate( const unsigned& current, MultiValue& myvals, std::vector<double>& buffer, std::vector<unsigned>& der_list ) const {}
  void finish( const std::vector<double>& buffer );
  bool applyForce( std::vector<double>& forces );
};

PLUMED_REGISTER_VESSEL(Moments,"MOMENTS")

void Moments::registerKeywords( Keywords& keys ) {
  Vessel::registerKeywords( keys );
  keys.remove("LABEL");
  keys.add("compulsory","COMPONENT","1","the component of the vector for which to calculate this quantity");
  keys.add("compulsory","MOMENTS","the list of moments that you would like to calculate");
}

void Moments::reserveKeyword( Keywords& keys ) {
  keys.reserve("optional","MOMENTS","calculate the moments of the distribution of collective variables. "
               "The mth moment of a distribution is calculated using \\f$\\frac{1}{N} \\sum_{i=1}^N ( s_i - \\overline{s} )^m \\f$, where \\f$\\overline{s}\\f$ is "
               "the average for the distribution. The moments keyword takes a lists of integers as input or a range. Each integer is a value of \\f$m\\f$. The final "
               "calculated values can be referenced using moment-\\f$m\\f$.  You can use the COMPONENT keyword in this action but the syntax is slightly different. "
               "If you would like the second and third moments of the third component you would use MOMENTS={COMPONENT=3 MOMENTS=2-3}.  The moments would then be referred to "
               "using the labels moment-3-2 and moment-3-3.  This syntax is also required if you are using numbered MOMENT keywords i.e. MOMENTS1, MOMENTS2...");
  keys.reset_style("MOMENTS","vessel");
  keys.addOutputComponent("moment","MOMENTS","the central moments of the distribution of values. The second moment "
                          "would be referenced elsewhere in the input file using "
                          "<em>label</em>.moment-2, the third as <em>label</em>.moment-3, etc.");
}

Moments::Moments( const vesselbase::VesselOptions& da) :
  Vessel(da)
{
  mystash = getAction()->buildDataStashes( NULL );
  ActionWithValue* a=dynamic_cast<ActionWithValue*>( getAction() );
  plumed_massert(a,"cannot create passable values as base action does not inherit from ActionWithValue");

  std::vector<std::string> moments; std::string valstr;
  if( getNumericalLabel()==0 ) {
    valstr = "moment-";
    moments=Tools::getWords(getAllInput(),"\t\n ,");
    Tools::interpretRanges(moments); mycomponent=1;
  } else {
    std::string numstr; parse("COMPONENT",mycomponent);
    Tools::convert( mycomponent, numstr); valstr = "moment-"  + numstr + "-";
    parseVector("MOMENTS",moments); Tools::interpretRanges(moments);
  }
  unsigned nn;
  for(unsigned i=0; i<moments.size(); ++i) {
    a->addComponentWithDerivatives( valstr + moments[i] );
    a->componentIsNotPeriodic( valstr + moments[i] );
    value_out.push_back( a->copyOutput( a->getNumberOfComponents()-1 ) );
    Tools::convert( moments[i], nn );
    if( nn<2 ) error("moments are only possible for m>=2" );
    powers.push_back( nn ); std::string num; Tools::convert(powers[i],num);
  }
}

void Moments::resize() {
  resizeBuffer(0);
  if( getAction()->derivativesAreRequired() ) {
    for(unsigned i=0; i<value_out.size(); ++i) value_out[i]->resizeDerivatives( getAction()->getNumberOfDerivatives() );
  }
}

std::string Moments::description() {
  std::string descri, num;
  Tools::convert(powers[0],num);
  if( getNumericalLabel()==0 ) {
    descri = "value " + getAction()->getLabel() + "." + "moment-" + num + " contains the " + num + "th moment of the distribution";
    for(unsigned i=1; i<powers.size(); ++i) {
      Tools::convert(powers[i],num);
      descri = descri + "\n  value " + getAction()->getLabel() + "." + "moment-" + num + " contains the " + num + "th moment of the distribution";
    }
  } else {
    std::string numlab; Tools::convert( mycomponent, numlab );
    descri = "value " + getAction()->getLabel() + "." + "moment-" + numlab + "-" + num + " contains the " + num + "th moment for the distribution of component " + numlab;
    for(unsigned i=1; i<powers.size(); ++i) {
      Tools::convert(powers[i],num);
      descri = descri + "\n  value " + getAction()->getLabel() + "." + "moment-" + numlab + "-" + num + " contains the " + num + "th moment for the distribution of component " + numlab;
    }
  }
  return descri;
}

void Moments::finish( const std::vector<double>& buffer ) {
  const double pi=3.141592653589793238462643383279502884197169399375105820974944592307;
  unsigned nvals=getAction()->getFullNumberOfTasks();
  std::vector<double>  myvalues( getAction()->getNumberOfQuantities() );

  double mean=0; Value myvalue;
  if( getAction()->isPeriodic() ) {
    std::string str_min, str_max; getAction()->retrieveDomain( str_min, str_max );
    double pfactor, min, max; Tools::convert(str_min,min); Tools::convert(str_max,max);
    pfactor = 2*pi / ( max-min ); myvalue.setDomain( str_min, str_max );
    double sinsum=0, cossum=0, val;
    for(unsigned i=0; i<mystash->getNumberOfStoredValues(); ++i) {
      mystash->retrieveSequentialValue( i, false, myvalues );
      val=pfactor*( myvalues[mycomponent] - min );
      sinsum+=sin(val); cossum+=cos(val);
    }
    mean = 0.5 + atan2( sinsum / static_cast<double>( nvals ), cossum / static_cast<double>( nvals ) ) / (2*pi);
    mean = min + (max-min)*mean;
  } else {
    for(unsigned i=0; i<mystash->getNumberOfStoredValues(); ++i) {
      mystash->retrieveSequentialValue( i, false, myvalues ); mean+=myvalues[mycomponent];
    }
    mean/=static_cast<double>( nvals ); myvalue.setNotPeriodic();
  }

  for(unsigned npow=0; npow<powers.size(); ++npow) {
    double dev1=0;
    if( value_out[0]->getNumberOfDerivatives()>0 ) {
      for(unsigned i=0; i<mystash->getNumberOfStoredValues(); ++i) {
        mystash->retrieveSequentialValue( i, false, myvalues );
        dev1+=pow( myvalue.difference( mean, myvalues[mycomponent] ), powers[npow] - 1 );
      }
      dev1/=static_cast<double>( nvals );
    }

    double moment=0;
    MultiValue myvals( getAction()->getNumberOfQuantities(), getAction()->getNumberOfDerivatives() ); myvals.clearAll();
    for(unsigned i=0; i<mystash->getNumberOfStoredValues(); ++i) {
      mystash->retrieveSequentialValue( i, false, myvalues );
      double tmp=myvalue.difference( mean, myvalues[mycomponent] );
      moment+=pow( tmp, powers[npow] );
      if( value_out[npow]->getNumberOfDerivatives() ) {
        double pref=pow( tmp, powers[npow] - 1 ) - dev1;
        mystash->retrieveDerivatives( i, false, myvals );
        for(unsigned j=0; j<myvals.getNumberActive(); ++j) {
          unsigned jatom=myvals.getActiveIndex(j);
          value_out[npow]->addDerivative(jatom, pref*myvals.getDerivative( mycomponent, jatom ) );
        }
        myvals.clearAll();
      }
    }
    if( value_out[npow]->getNumberOfDerivatives()>0 ) value_out[npow]->chainRule( powers[npow] / static_cast<double>( nvals ) );
    value_out[npow]->set( moment / static_cast<double>( nvals ) );
  }
}

bool Moments::applyForce( std::vector<double>& forces ) {
  std::vector<double> tmpforce( forces.size() );
  forces.assign(forces.size(),0.0); bool wasforced=false;
  for(unsigned i=0; i<value_out.size(); ++i) {
    if( value_out[i]->applyForce( tmpforce ) ) {
      wasforced=true;
      for(unsigned j=0; j<forces.size(); ++j) forces[j]+=tmpforce[j];
    }
  }
  return wasforced;
}

}
}
