/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2017-2020 The plumed team
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
#include "ActionToFetchData.h"
#include "PlumedMain.h"
#include "Atoms.h"

namespace PLMD {

template <class T>
class OutputDataObjectTyped : public OutputDataObject {
private:
/// A pointer to the place we are setting the data
  T* data_out;
public:
  static std::unique_ptr<OutputDataObject> create(unsigned n);
/// Set the pointer to the output
  void setPointer( void* outval ) override;
/// This transfers everything to the output
  void setData( const std::vector<double>& data ) override;
};

std::unique_ptr<OutputDataObject> OutputDataObject::create(unsigned n) {
  if(n==sizeof(double)) {
    return std::unique_ptr<OutputDataObject>(new OutputDataObjectTyped<double>());
  } else  if(n==sizeof(float)) {
    return std::unique_ptr<OutputDataObject>(new OutputDataObjectTyped<float>());
  } 
  std::string pp; Tools::convert(n,pp);
  plumed_merror("cannot create an MD interface with sizeof(real)=="+ pp);
  return NULL;
}

template <class T>
void OutputDataObjectTyped<T>::setPointer( void* outval ) {
   data_out=static_cast<T*>(outval); 
}

template <class T>
void OutputDataObjectTyped<T>::setData( const std::vector<double>& data ) {
   for(unsigned i=0;i<data.size();++i) data_out[i] = static_cast<T>( data[i] );
}

PLUMED_REGISTER_ACTION(ActionToFetchData,"FETCH")

void ActionToFetchData::registerKeywords(Keywords& keys) {
   Action::registerKeywords(keys); ActionPilot::registerKeywords(keys); ActionWithArguments::registerKeywords(keys);
   keys.add("compulsory","STRIDE","1","the frequency with which the quantities of interest should be stored");
   keys.add("compulsory","TYPE","value","what do you want to collect for the value can be derivative/force");
   keys.use("ARG");
} 

ActionToFetchData::ActionToFetchData(const ActionOptions&ao):
Action(ao),
ActionPilot(ao),
ActionWithArguments(ao),
mydata(OutputDataObject::create(plumed.getAtoms().getRealPrecision()))
{
   std::string type; parse("TYPE",type);
   if( type=="value" ) gtype=val;
   else if( type=="derivatives" ) gtype=deriv;
   else if( type=="forces" ) gtype=force;
   else plumed_merror("cannot get " + type + " for value TYPE should be value/derivative/force");

   if( gtype!=val ) error("not implemented functionality to pass derviatives or forces to python.  Email gareth.tribello@gmail.com if you want this.");

   if( getNumberOfArguments()!=1 ) error("python interface works best when you ask for one argument at a time");
   getPntrToArgument(0)->buildDataStore( getLabel() ); data.resize( getPntrToArgument(0)->getNumberOfValues(getLabel()) );
}

void ActionToFetchData::get_rank(long* rank ) {
   rank[0]=getPntrToArgument(0)->getRank();
}

void ActionToFetchData::get_shape(long* dims ) {
   if( getPntrToArgument(0)->getRank()==0 ) { dims[0]=1; return; }
   for(unsigned j=0;j<getPntrToArgument(0)->getRank();++j) dims[j] = getPntrToArgument(0)->getShape()[j];
}

void ActionToFetchData::set_memory(void* val ) {
   mydata->setPointer(val);
}

void ActionToFetchData::calculate() {
   plumed_assert( gtype==val ); Value* val = getPntrToArgument(0);
   for(unsigned i=0;i<data.size();++i) data[i] = val->getRequiredValue( getLabel(), i );
   mydata->setData( data );
}
 
}
