/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2020 The plumed team
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
#include "ActionSetup.h"
#include "tools/Exception.h"
#include "tools/OFile.h"
#include "tools/OpenMP.h"
#include "PlumedMain.h"

namespace PLMD {

Value::Value():
  action(NULL),
  value_set(false),
  norm(1.0),
  hasForce(false),
  hasDeriv(true),
  ngrid_der(0),
  historydependent(false),
  istimeseries(false),
  shape(std::vector<unsigned>()),
  ntasks(1),
  reducedTasks(false),
  constant(false),
  alwaysstore(false),
  storedata(true),
  neverstore(false),
  symmetric(false),
  streampos(0),
  periodicity(unset),
  min(0.0),
  max(0.0),
  max_minus_min(0.0),
  inv_max_minus_min(0.0),
  derivativeIsZeroWhenValueIsZero(false)
{
  data.resize(1); inputForces.resize(1);
}

Value::Value(const std::string& name):
  action(NULL),
  value_set(false),
  hasForce(false),
  hasDeriv(true),
  name(name),
  ngrid_der(0),
  historydependent(false),
  istimeseries(false),
  ntasks(1),
  reducedTasks(false),
  constant(false),
  alwaysstore(false),
  storedata(true),
  neverstore(false),
  symmetric(false),
  streampos(0),
  periodicity(unset),
  min(0.0),
  max(0.0),
  max_minus_min(0.0),
  inv_max_minus_min(0.0),
  derivativeIsZeroWhenValueIsZero(false)
{
  data.resize(1); inputForces.resize(1); 
  data[0]=inputForces[0]=0;
}

Value::Value(ActionWithValue* av, const std::string& name, const bool withderiv,const std::vector<unsigned>&ss):
  action(av),
  value_set(false),
  norm(1.0),
  hasForce(false),
  name(name),
  hasDeriv(withderiv),
  ngrid_der(0),
  historydependent(false),
  istimeseries(false),
  ntasks(1),
  reducedTasks(false),
  constant(false),
  alwaysstore(false),
  storedata(true),
  neverstore(false),
  symmetric(false),
  streampos(0),
  periodicity(unset),
  min(0.0),
  max(0.0),
  max_minus_min(0.0),
  inv_max_minus_min(0.0),
  derivativeIsZeroWhenValueIsZero(false)
{
  if( action ) alwaysstore=action->getName()=="PUT";
  if( ss.size()>0 && !alwaysstore ) storedata=false;
  setShape( ss ); 
}

void Value::setShape( const std::vector<unsigned>&ss ) {
  shape.resize( ss.size() );
  for(unsigned i=0; i<shape.size(); ++i) shape[i]=ss[i];
  // Set the number of tasks here that are done in loop
  if( hasDeriv && shape.size()>0 ) ntasks=getNumberOfValues();
  else if( shape.size()>0 ) ntasks=shape[0];
  // Matrices are resized dynamically so we can use the sparsity pattern to reduce the 
  // overhead on memory
  if( getRank()==2 ) {
      if( !hasDeriv && !alwaysstore && !istimeseries ) return;
      if( !hasDeriv && !alwaysstore && istimeseries && getNumberOfColumns()==0 ) return;
  }

  if( (hasDeriv || storedata || istimeseries) && shape.size()>0 ) {
     unsigned ss=getSize(); if( data.size()!=ss ) data.resize(ss); 
     unsigned fsize=1; for(unsigned i=0; i<shape.size(); ++i) fsize *= shape[i];
     if( fsize!=inputForces.size() ) inputForces.resize(fsize); 
     if( hasDeriv ) {
         ngrid_der = shape.size();
         if( action ) ngrid_der = action->getNumberOfDerivatives();
     } 
  } else if( hasDeriv ) {
     data.resize( 1 +  action->getNumberOfDerivatives() ); inputForces.resize(1);
  } else if( shape.size()==0 ) { data.resize(1); inputForces.resize(1); }
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

void Value::use( Action* act, std::vector<Value*>& arg ) {
  arg.push_back( this ); ActionWithArguments* aa=dynamic_cast<ActionWithArguments*>( act );
  if( aa ) userdata.insert(aa->getLabel());
}

void Value::buildDataStore( const std::string& actlabel ) {
  if( neverstore ) return ;
  bool found=false; 
  for(unsigned i=0; i<store_data_for.size(); ++i) {
    if( actlabel==store_data_for[i].first ) found=true;
  }
  if( !found ) store_data_for.push_back( std::pair<std::string,int>(actlabel,-1) );
  storedata=true; setShape( shape );
}

void Value::setConstant() {
  constant=true; storedata=true; setShape( shape );
}

void Value::alwaysStoreValues() {
  plumed_assert( !neverstore); alwaysstore=true; storedata=true; setShape( shape );
}

void Value::makeHistoryDependent() {
  historydependent=true;
  if( shape.size()>0 && !hasDeriv && action->getName()!="AVERAGE" ) { istimeseries=true; setShape( shape ); }
}

void Value::neverStoreValues() {
  plumed_assert( !alwaysstore ); neverstore=true;
}

bool Value::isPeriodic()const {
  plumed_massert(periodicity!=unset,"periodicity should be set");
  return periodicity==periodic;
}

bool Value::applyForce( std::vector<double>& forces ) const {
  if( !hasForce || constant ) return false;

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

void Value::setGradients( ActionAtomistic* aa, unsigned& start ) {
  // Can't do gradients if we don't have derivatives
  if( !hasDeriv ) return;
  plumed_assert( shape.size()==0 );
  for(unsigned j=0; j<aa->getNumberOfAtoms(); ++j) {
      Vector der(data[1+start+3*j],data[1+start+3*j+1],data[1+start+3*j+2]);
      aa->getGradient( j, der, gradients );
  }
  start += aa->getNumberOfAtoms();
}

void Value::passGradients( const double& der, std::map<AtomNumber,Vector>& g ) const {
  for(const auto & p : gradients) { AtomNumber iatom=p.first; g[iatom] += p.second*der; }
}

double Value::projection(const Value& v1,const Value&v2) {
  double proj=0.0; plumed_assert( v1.getRank()==0 && v2.getRank()==0 );
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

unsigned Value::getSize() const {
  unsigned size=getNumberOfValues();
  if( shape.size()>0 && hasDeriv && action ) return size*( 1 + action->getNumberOfDerivatives() ); 
  else if( shape.size()>0 && hasDeriv ) return size*( 1 + shape.size() );
  return size;
}

unsigned Value::getNumberOfValues() const {
  unsigned size=1; for(unsigned i=0; i<shape.size(); ++i) size *= shape[i];
  return size;
}

double Value::get(const unsigned& ival, const bool trueind) const {
  if( hasDeriv ) return data[ival*(1+ngrid_der)] / norm;
#ifdef DNDEBUG 
  if( action ) plumed_dbg_massert( ival<getNumberOfValues(), "could not get value from " + name );
#endif
  if( shape.size()==2 && getNumberOfColumns()<shape[1] && trueind ) {
      unsigned irow = std::floor( ival / shape[0] ), jcol = ival%shape[0];
      // This is a special treatment for the lower triangular matrices that are used when 
      // we do ITRE with COLLECT_FRAMES
      if( getNumberOfColumns()==0 ) {
          if( jcol<=irow ) return data[0.5*irow*(irow+1) + jcol] / norm;
          return 0;
      }
      for(unsigned i=0; i<getRowLength(irow); ++i) {
          if( getRowIndex(irow,i)==jcol ) return data[irow*getNumberOfColumns()+i] / norm;
      }
      return 0.0;
  }
  plumed_massert( ival<data.size(), "cannot get value from " + name );
  if( norm>epsilon ) return data[ival] / norm;
  return 0.0;
}

void Value::addForce(const unsigned& iforce, double f, const bool trueind) {
  if( action->getName()=="COLLECT_FRAMES" ) return; hasForce=true;
  if( shape.size()==2 && !hasDeriv && getNumberOfColumns()<shape[1] && trueind ) { 
      unsigned irow = std::floor( iforce / shape[0] ), jcol = iforce%shape[0];
      for(unsigned i=0; i<getRowLength(irow); ++i) {
          if( getRowIndex(irow,i)==jcol ) { inputForces[irow*getNumberOfColumns()+i]+=f; return; }
      }
      plumed_assert( fabs(f)<epsilon ); return;
  } 
  plumed_massert( iforce<inputForces.size(), "can't add force to " + name );
  inputForces[iforce]+=f;
}

double Value::getGridDerivative(const unsigned& n, const unsigned& j ) const {
  plumed_dbg_assert( hasDeriv && n*(1+ngrid_der) + 1 + j < data.size() );
  return data[n*(1+ngrid_der) + 1 + j] / norm;
}

void Value::setGridDerivative(const unsigned& n, const unsigned& j, const double& val ) {
  plumed_dbg_assert( hasDeriv && n*(1+shape.size()) + 1 + j < data.size() );
  data[n*(1+shape.size()) + 1 + j] = val;
}

void Value::print( const std::string& uselab, OFile& ofile ) const {
  if( isPeriodic() ) { ofile.printField( "min_" + name, str_min ); ofile.printField("max_" + name, str_max ); }
  if( shape.size()==0 ) {
    ofile.printField( name, get(0) );
  } else {
    std::vector<unsigned> indices( shape.size() );
    for(unsigned i=0; i<getNumberOfValues(); ++i) {
      convertIndexToindices( i, indices ); std::string num, fname = name;
      for(unsigned i=0; i<shape.size(); ++i) { Tools::convert( indices[i]+1, num ); fname += "." + num; }
      ofile.printField( fname, get(i) );
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

unsigned Value::getNumberOfColumns() const {
  plumed_massert( shape.size()==2 && !hasDeriv, "failing in " + name );
  if( alwaysstore ) return shape[1];
  return action->getNumberOfColumns();
}

unsigned Value::getRowLength( const unsigned& irow ) const {
  if( alwaysstore || matindexes[(1+getNumberOfColumns())*irow]>shape[1] ) return shape[1];
  return matindexes[(1+getNumberOfColumns())*irow];
}
 
unsigned Value::getRowIndex( const unsigned& irow, const unsigned& jind ) const {
  if( !alwaysstore ) return matindexes[(1+getNumberOfColumns())*irow+1+jind];
  return jind;
}

void Value::reshapeMatrixStore() {
  plumed_assert( shape.size()==2 && !hasDeriv && storedata && action );
  unsigned size = shape[0]*getNumberOfColumns();
  if( matindexes.size()==(size+shape[0]) ) return;
  data.resize( size ); inputForces.resize( size ); 
  matindexes.resize( size + shape[0] );
}

const std::vector<unsigned>& Value::getShape() const {
  return shape;
}

void Value::set(const unsigned& n, const double& v ) {
  value_set=true; 
  if( getRank()==0 ) { plumed_assert( n==0 ); data[n]=v; applyPeriodicity(n); }
  else if( !hasDeriv ) { plumed_dbg_massert( n<data.size(), "failing in " + getName() ); data[n]=v; applyPeriodicity(n); }
  else { data[n*(1+ngrid_der)] = v; }
}

std::string Value::getOutputDescription() const {
  if( getRank()==0 ) return " scalar labelled " + name;
  if( hasDerivatives() ) return " grid labelled " + name;
  if( getRank()==1 ) return " vector labelled " + name;
  return " matrix labelled " + name;
}

void Value::setSymmetric( const bool& sym ) {
  plumed_assert( shape.size()==2 && !hasDeriv ); 
  if( sym && shape[0]!=shape[1] ) plumed_merror("non-square matrix cannot be symmetric");
  symmetric=sym; 
}

bool Value::isSymmetric() const {
  return symmetric;
}

void Value::writeBinary(std::ostream&o) const {
  o.write(reinterpret_cast<const char*>(&data[0]),data.size()*sizeof(double));
}

void Value::readBinary(std::istream&i) {
  i.read(reinterpret_cast<char*>(&data[0]),data.size()*sizeof(double));
}

unsigned Value::getGoodNumThreads( const unsigned& j, const unsigned& k ) const {
  return OpenMP::getGoodNumThreads( &data[j], (k-j) );
}

bool Value::getDerivativeIsZeroWhenValueIsZero() const {
  return derivativeIsZeroWhenValueIsZero;
}

void Value::setNumberOfTasks( const unsigned& nt ) {
  ntasks=nt;
}

unsigned Value::getNumberOfTasks() const {
  return ntasks;
}

}
