/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2019 The plumed team
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
  created_in_plumedmain(false),
  action(NULL),
  value_set(false),
  reset(true),
  norm(1.0),
  hasForce(false),
  hasDeriv(true),
  shape(std::vector<unsigned>()),
  alwaysstore(false),
  storedata(true),
  neverstore(false),
  columnsums(false),
  symmetric(false),
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
  created_in_plumedmain(false),
  action(av),
  value_set(false),
  reset(true),
  norm(1.0),
  hasForce(false),
  name(name),
  hasDeriv(withderiv),
  alwaysstore(false),
  neverstore(false),
  columnsums(false),
  symmetric(false),
  bufstart(0),
  streampos(0),
  periodicity(unset),
  min(0.0),
  max(0.0),
  max_minus_min(0.0),
  inv_max_minus_min(0.0)
{
  setShape( ss );
}

void Value::setShape( const std::vector<unsigned>&ss ) {
  shape.resize( ss.size() );
  for(unsigned i=0; i<shape.size(); ++i) shape[i]=ss[i];

  if( ss.size()>0 && !alwaysstore ) storedata=false;
  else storedata=true;

  if( shape.size()>0 ) data.resize(getSize());
  else if( hasDeriv ) data.resize( 1 +  action->getNumberOfDerivatives() );
  else data.resize(1);
  unsigned fsize=1; for(unsigned i=0; i<shape.size(); ++i) fsize *= shape[i];
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

void Value::buildDataStore( const std::string& actlabel ) {
  if( neverstore ) return ;
  bool found=false;
  for(unsigned i=0; i<store_data_for.size(); ++i) {
    if( actlabel==store_data_for[i].first ) found=true;
  }
  if( !found ) store_data_for.push_back( std::pair<std::string,int>(actlabel,-1) );
  storedata=true;
}

void Value::alwaysStoreValues() {
  plumed_assert( !neverstore); alwaysstore=true; storedata=true;
}

void Value::neverStoreValues() {
  plumed_assert( !alwaysstore ); neverstore=true;
}

void Value::buildColumnSums() {
  columnsums=true; storedata=true;
}

void Value::interpretDataRequest( const std::string& uselab, unsigned& nargs, std::vector<Value*>& args, const std::string& values ) {
  bool found=false;
  if( shape.size()>0 ) {
    for(unsigned i=0; i<args.size(); ++i) {
      if( this==args[i] ) { found=true; break; }
    }
  }
  if( !found ) { args.push_back(this); }

  if( userdata.count(uselab) ) {
    if( values=="" ) { nargs++; return; }
    if( userdata[uselab][0].first<0 ) plumed_merror("cannot mix use of specific items from value and all items in a single action");
  } else {
    userdata.insert( std::pair<std::string,std::vector<std::pair<int,int> > >(uselab,std::vector<std::pair<int,int> >()) );
    if( values=="*" ) { plumed_merror("invalid use of wildcard"); return; }
    else if( values=="" ) { userdata[uselab].push_back(std::pair<int,int>(-1,nargs) ); nargs++; return; }
  }
  // Retrieve the indices of the point from the string requesting the index
  std::vector<unsigned> indices( shape.size() ); std::string indstr=values;
  for(unsigned i=0; i<shape.size()-1; ++i) {
    std::size_t dot = indstr.find_first_of(".");
    if( dot==std::string::npos ) action->error("invalid specification for element of value"); 
    Tools::convert( indstr.substr(0,dot), indices[i] );
    indices[i] -= 1; indstr=indstr.substr(dot+1);
  }
  Tools::convert( indstr, indices[indices.size()-1] ); indices[indices.size()-1] -= 1;
  if( getIndex(indices)>=getNumberOfValues( action->getLabel() ) ) action->error("action does not have this many components");
  userdata[uselab].push_back( std::pair<int,int>(getIndex(indices),nargs) ); nargs++;
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

  if( shape.size()==0 && hasDeriv ) {
    for(unsigned i=0; i<forces.size(); ++i) forces[i] += inputForces[0]*data[1 + i];
  } else if( hasDeriv ) {
    plumed_merror("Cannot apply forces to grids in this way"); return false;
  } else if( shape.size()>0 ) {
    return false;
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
  if( created_in_plumedmain ) return NULL;
  plumed_assert( action!=NULL );
  return action;
}

unsigned Value::getSize() const {
  unsigned size=getNumberOfValues( name );
  if( shape.size()>0 && hasDeriv ) return size*( 1 + action->getNumberOfDerivatives() );
  return size;
}

unsigned Value::getNumberOfValues( const std::string& alab ) const {
  if( usingAllVals(alab) ) {
    unsigned size=1; for(unsigned i=0; i<shape.size(); ++i) size *= shape[i];
    return size;
  } else {
    return userdata.find(alab)->second.size();
  }
}

double Value::get(const unsigned& ival) const {
  if( hasDeriv ) return data[ival*(1+action->getNumberOfDerivatives())] / norm;
#ifdef DNDEBUG 
  if( action ) plumed_dbg_massert( ival<getNumberOfValues( action->getLabel() ), "could not get value from " + name );
#endif
  if( norm>epsilon ) return data[ival] / norm;
  return 0.0;
}

double Value::getGridDerivative(const unsigned& n, const unsigned& j ) const {
  plumed_dbg_assert( hasDeriv && n*(1+action->getNumberOfDerivatives()) + 1 + j < data.size() );
  return data[n*(1+action->getNumberOfDerivatives()) + 1 + j] / norm;
}

void Value::print( const std::string& uselab, OFile& ofile ) const {
  plumed_dbg_assert( userdata.count(uselab) );
  if( shape.size()==0 ) {
    if( isPeriodic() ) { ofile.printField( "min_" + name, str_min ); ofile.printField("max_" + name, str_max ); }
    ofile.printField( name, get(0) );
  } else if( userdata.find(uselab)->second[0].first<0 ) {
    if( isPeriodic() ) { ofile.printField( "min_" + name, str_min ); ofile.printField("max_" + name, str_max ); }
    std::vector<unsigned> indices( shape.size() );
    for(unsigned i=0; i<getNumberOfValues(uselab); ++i) {
      convertIndexToindices( i, indices ); std::string num, fname = name;
      for(unsigned i=0; i<shape.size(); ++i) { Tools::convert( indices[i]+1, num ); fname += "." + num; }
      ofile.printField( fname, get(i) );
    }
  } else {
    if( isPeriodic() ) { ofile.printField( "min_" + name, str_min ); ofile.printField("max_" + name, str_max ); }
    std::vector<unsigned> indices( shape.size() );
    for(unsigned i=0; i<userdata.find(uselab)->second.size(); ++i) {
      convertIndexToindices( userdata.find(uselab)->second[i].first, indices ); std::string num, fname = name;
      for(unsigned i=0; i<shape.size(); ++i) { Tools::convert( indices[i]+1, num ); fname += "." + num; }
      ofile.printField( fname, get( userdata.find(uselab)->second[i].first ) );
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

void Value::set(const unsigned& n, const double& v ) {
  value_set=true;
  if( getRank()==0 ) { plumed_assert( n==0 ); data[n]=v; applyPeriodicity(n); }
  else if( !hasDeriv ) { data[n]=v; applyPeriodicity(n); }
  else { data[n*(1+action->getNumberOfDerivatives())] = v; }
}

bool Value::usingAllVals( const std::string& alabel ) const {
  if( !userdata.count(alabel) ) return true;
  if( userdata.find(alabel)->second.size()==0 ) return true;
  if( userdata.find(alabel)->second[0].first<0 ) return true;
  return false;
}

double Value::getRequiredValue(  const std::string& alabel, const unsigned& num  ) const {
  if( usingAllVals(alabel) ) return get(num);
  return get( userdata.find(alabel)->second[num].first );
}

void Value::getRequiredValue(  const std::string& alabel, const unsigned& num, std::vector<double>& args ) const {
  if( usingAllVals(alabel) ) {
    args[userdata.find(alabel)->second[0].second+num] = getRequiredValue( alabel, num );
  } else {
    args[userdata.find(alabel)->second[num].second] = getRequiredValue( alabel, num );
  }
}

std::string Value::getOutputDescription( const std::string& alabel ) const {
  if( getRank()==0 ) return " " + name;

  if( usingAllVals(alabel) ) {
    if( hasDerivatives() ) return " grid labelled " + name;
    if( getRank()==1 ) return " vector labelled " + name;
    if( getRank()==2 ) return " matrix labelled " + name;
  }
  // N.B. Output for rank 2 values in this case is not very transparent.
  std::string datp;
  for(unsigned i=0; i<userdata.find(alabel)->second.size(); ++i) datp += " " + getOutputDescription( alabel, i ); 
  return datp;
}

std::string Value::getOutputDescription( const std::string& alabel, const unsigned& i ) const {
  if( getRank()==1 ) {
      std::string num; Tools::convert( userdata.find(alabel)->second[i].first+1, num ); return name + "." + num;
  } else if( getRank()==2 ) {
      std::string num; Tools::convert( std::floor(userdata.find(alabel)->second[i].first/shape[1])+1 , num );
      std::string num2; Tools::convert( userdata.find(alabel)->second[i].first%shape[1]+1, num2 );
      return name + "." + num + "." + num2;
  } else {
      action->error("elements of three rank objects have not been implemented");
  }
}

void Value::setSymmetric( const bool& sym ) {
  plumed_assert( shape.size()==2 && !hasDeriv ); 
  if( sym && shape[0]!=shape[1] ) plumed_merror("non-square matrix cannot be symmetric");
  symmetric=sym; 
}

bool Value::isSymmetric() const {
  return symmetric;
}

// void Value::setBufferPosition( const unsigned& ibuf ){
//   bufstart = ibuf;
// }

// unsigned Value::getBufferPosition() const {
//   return bufstart;
// }

}
