/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2017 The plumed team
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
#include "Value.h"
#include "ActionWithValue.h"
#include "ActionAtomistic.h"
#include "ActionWithArguments.h"
#include "ActionWithVirtualAtom.h"
#include "tools/Exception.h"
#include "tools/OFile.h"
#include "Atoms.h"
#include "PlumedMain.h"

namespace PLMD {

Value::Value():
  action(NULL),
  value_set(false),
  reset(true),
  hasForce(false),
  hasDeriv(true),
  shape(std::vector<unsigned>()),
  storedata(true),
  bufstart(0),
  streampos(0),
  periodicity(unset),
  min(0.0),
  max(0.0),
  max_minus_min(0.0),
  inv_max_minus_min(0.0)
{
  data.resize(1); inputForces.resize(1);
}

Value::Value(ActionWithValue* av, const std::string& name, const bool withderiv,const std::vector<unsigned>&ss):
  action(av),
  value_set(false),
  reset(true),
  hasForce(false),
  name(name),
  hasDeriv(withderiv),
  shape(ss),
  storedata(shape.size()==0),
  bufstart(0),
  streampos(0),
  periodicity(unset),
  min(0.0),
  max(0.0),
  max_minus_min(0.0),
  inv_max_minus_min(0.0)
{
  data.resize(getSize());
  unsigned fsize=1; for(unsigned i=0;i<shape.size();++i) fsize *= shape[i];
  inputForces.resize( fsize );
}

void Value::setupPeriodicity() {
  if( min==0 && max==0 ) {
    periodicity=notperiodic;
  } else {
    periodicity=periodic;
    max_minus_min=max-min;
    plumed_massert(max_minus_min>0, "your function has a very strange domain?");
    inv_max_minus_min=1.0/max_minus_min;
  }
}

void Value::buildDataStore(){
  storedata=true;
}

void Value::interpretDataRequest( const std::string& uselab, const std::string& values ){
  if( userdata.count(uselab) ){
      if( values=="" ){ return; }
      if( userdata[uselab][0]<0 ) plumed_merror("cannot mix use of specific items from value and all items in a single action");
  } else {
      userdata.insert( std::pair<std::string,std::vector<int> >(uselab,std::vector<int>()) );
      if( values=="*" ){ plumed_merror("invalid use of wildcard"); return; }
      else if( values=="" ){ userdata[uselab].push_back(-1); return; }
  }
  // Retrieve the indices of the point from the string requesting the index
  std::vector<unsigned> indices( shape.size() ); std::string indstr=values;
  for(unsigned i=0;i<shape.size()-1;++i){
      std::size_t dot = indstr.find_first_of(".");
      Tools::convert( indstr.substr(0,dot), indices[i] );
      indices[i] -= 1; indstr=indstr.substr(dot+1);
  }
  Tools::convert( indstr, indices[indices.size()-1] );
  userdata[uselab].push_back( getIndex(indices) );
}

// void Value::addStreamIndex( const int& newi ){
//   plumed_dbg_assert( shape.size()>0 );
//   if( indices_in_stream.size()>0 ) plumed_assert( indices_in_stream[0]!=-1 );
//   indices_in_stream.push_back( newi );
// }

bool Value::isPeriodic()const {
  plumed_massert(periodicity!=unset,"periodicity should be set");
  return periodicity==periodic;
}

bool Value::applyForce( std::vector<double>& forces ) const {
  if( !hasForce ) return false;

  if( shape.size()==0 && hasDeriv ){
      for(unsigned i=0;i<forces.size();++i) forces[i] += inputForces[0]*data[1 + i];
  } else if( hasDeriv ){
      const unsigned N=action->getNumberOfDerivatives();
      for(unsigned i=0;i<inputForces.size();++i){
         for(unsigned j=0;j<N;++j) forces[j] += inputForces[i]*data[ i*(1+N) ];
      }
  } else if( shape.size()>0 ) {
      plumed_merror("should not be using apply forces with vectors");
  }
  return true;
}

void Value::setNotPeriodic() {
  min=0; max=0; periodicity=notperiodic;
}

bool Value::hasDerivatives() const {
  if( shape.size()==0 && action->doNotCalculateDerivatives() ) return false;
  return hasDeriv;
}

void Value::setDomain(const std::string& pmin,const std::string& pmax) {
  str_min=pmin;
  if( !Tools::convert(str_min,min) ) action->error("could not convert period string " + str_min + " to real");
  str_max=pmax;
  if( !Tools::convert(str_max,max) ) action->error("could not convert period string " + str_max + " to read");
  setupPeriodicity();
}

void Value::getDomain(std::string&minout,std::string&maxout) const {
  plumed_massert(periodicity==periodic,"function should be periodic");
  minout=str_min;
  maxout=str_max;
}

void Value::getDomain(double&minout,double&maxout) const {
  plumed_massert(periodicity==periodic,"function should be periodic");
  minout=min;
  maxout=max;
}

void Value::setGradients() {
  // Can't do gradients if we don't have derivatives
  if( !hasDeriv ) return;
  plumed_assert( shape.size()==0 );
  gradients.clear();
  ActionAtomistic*aa=dynamic_cast<ActionAtomistic*>(action);
  ActionWithArguments*aw=dynamic_cast<ActionWithArguments*>(action);
  if(aa) {
    Atoms&atoms((aa->plumed).getAtoms());
    for(unsigned j=0; j<aa->getNumberOfAtoms(); ++j) {
      AtomNumber an=aa->getAbsoluteIndex(j);
      if(atoms.isVirtualAtom(an)) {
        ActionAtomistic* a=atoms.getVirtualAtomsAction(an);
        for(const auto & p : a->getVatomGradients(an)) {
// controllare l'ordine del matmul:
          gradients[p.first]+=matmul(Vector(data[1+3*j],data[1+3*j+1],data[1+3*j+2]),p.second);
        }
      } else {
        for(unsigned i=0; i<3; i++) gradients[an][i]+=data[1+3*j+i];
      }
    }
  } else if(aw) {
    std::vector<Value*> values=aw->getArguments();
    for(unsigned j=0; j<data.size()-1; j++) {
      for(const auto & p : values[j]->gradients) {
        AtomNumber iatom=p.first;
        gradients[iatom]+=p.second*data[1+j];
      }
    }
  } else plumed_error();
}

double Value::projection(const Value& v1,const Value&v2) {
  double proj=0.0;
  const std::map<AtomNumber,Vector> & grad1(v1.gradients);
  const std::map<AtomNumber,Vector> & grad2(v2.gradients);
  for(const auto & p1 : grad1) {
    AtomNumber a=p1.first;
    const auto p2=grad2.find(a);
    if(p2!=grad2.end()) {
      proj+=dotProduct(p1.second,(*p2).second);
    }
  }
  return proj;
}

ActionWithValue* Value::getPntrToAction() {
  plumed_assert( action!=NULL );
  return action;
}

void Value::activateTasks( std::vector<unsigned>& taskFlags ) const {
  //plumed_dbg_assert( indices_in_stream.size()>0 );
  //if( indices_in_stream[0]==-1 ){
  //    for(unsigned i=0;i<taskFlags.size();++i) taskFlags[i]=1;
  //}
  //for(unsigned i=0;i<indices_in_stream.size();++i) taskFlags[ indices_in_stream[i] ] = 1;
}

unsigned Value::getSize() const {
  unsigned size=getNumberOfValues(); 
  if( shape.size()>0 && hasDeriv ) return size*( 1 + action->getNumberOfDerivatives() );
  return size;
}

unsigned Value::getNumberOfValues() const {
  unsigned size=1; for(unsigned i=0;i<shape.size();++i) size *= shape[i];
  return size;
}

double Value::get(const unsigned& ival) const {
  if( hasDeriv ) return data[ival*(1+action->getNumberOfDerivatives())];
  return data[ival];
} 

void Value::print( const std::string& uselab, OFile& ofile ) const {
  plumed_dbg_assert( userdata.count(uselab) );
  if( shape.size()==0 ){
      if( isPeriodic() ){ ofile.printField( "min_" + name, str_min ); ofile.printField("max_" + name, str_max ); } 
      ofile.printField( name, data[0] ); 
  } else if( userdata.find(uselab)->second[0]<0 ){
      if( isPeriodic() ){ ofile.printField( "min_" + name, str_min ); ofile.printField("max_" + name, str_max ); }
      std::vector<unsigned> indices( shape.size() ); 
      for(unsigned i=0;i<action->getFullNumberOfTasks();++i){
          convertIndexToindices( i, indices ); std::string num, fname = name; 
          for(unsigned i=0;i<shape.size();++i){ Tools::convert( indices[i]+1, num ); fname += "." + num; }
          ofile.printField( fname, get(i) ); 
      }
  } else {
      if( isPeriodic() ){ ofile.printField( "min_" + name, str_min ); ofile.printField("max_" + name, str_max ); }
      std::vector<unsigned> indices( shape.size() ); 
      for(unsigned i=0;i<userdata.find(uselab)->second.size();++i){ 
         convertIndexToindices( userdata.find(uselab)->second[i], indices ); std::string num, fname = name;
         for(unsigned i=0;i<shape.size();++i){ Tools::convert( indices[i]+1, num ); fname += "." + num; }
         ofile.printField( fname, get( userdata.find(uselab)->second[i] ) );
      }
  }
}

//void Value::setPositionInStream( const unsigned& istream ){
//  streampos = istream; 
//}

unsigned Value::getPositionInStream() const {
  return streampos; 
}

unsigned Value::getPositionInMatrixStash() const {
  return matpos;
}

const std::vector<unsigned>& Value::getShape() const {
  return shape;
}

// void Value::setBufferPosition( const unsigned& ibuf ){
//   bufstart = ibuf;
// }

// unsigned Value::getBufferPosition() const {
//   return bufstart;
// }

}
