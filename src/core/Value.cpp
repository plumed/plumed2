/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2023 The plumed team
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
#include "ActionWithVector.h"
#include "ActionWithVirtualAtom.h"
#include "tools/Exception.h"
#include "tools/OpenMP.h"
#include "tools/OFile.h"
#include "PlumedMain.h"

namespace PLMD {

Value::Value():
  action(NULL),
  value_set(false),
  hasForce(false),
  storedata(false),
  shape(std::vector<unsigned>()),
  hasDeriv(true),
  bufstart(0),
  streampos(0),
  matpos(0),
  ngrid_der(0),
  ncols(0),
  book_start(0),
  symmetric(false),
  periodicity(unset),
  min(0.0),
  max(0.0),
  max_minus_min(0.0),
  inv_max_minus_min(0.0),
  derivativeIsZeroWhenValueIsZero(false)
{
  data.resize(1); inputForce.resize(1);
}

Value::Value(const std::string& name):
  action(NULL),
  value_set(false),
  hasForce(false),
  name(name),
  storedata(false),
  shape(std::vector<unsigned>()),
  hasDeriv(true),
  bufstart(0),
  streampos(0),
  ngrid_der(0),
  matpos(0),
  ncols(0),
  book_start(0),
  symmetric(false),
  periodicity(unset),
  min(0.0),
  max(0.0),
  max_minus_min(0.0),
  inv_max_minus_min(0.0),
  derivativeIsZeroWhenValueIsZero(false)
{
  data.resize(1); inputForce.resize(1);
  data[0]=inputForce[0]=0;
}

Value::Value(ActionWithValue* av, const std::string& name, const bool withderiv, const std::vector<unsigned>&ss):
  action(av),
  value_set(false),
  hasForce(false),
  name(name),
  storedata(false),
  hasDeriv(withderiv),
  bufstart(0),
  streampos(0),
  ngrid_der(0),
  matpos(0),
  ncols(0),
  book_start(0),
  symmetric(false),
  periodicity(unset),
  min(0.0),
  max(0.0),
  max_minus_min(0.0),
  inv_max_minus_min(0.0),
  derivativeIsZeroWhenValueIsZero(false)
{
  if( action ) {
    if( action->getName()=="ACCUMULATE" || action->getName()=="COLLECT" ) valtype=average;
  }
  if( action ) storedata=action->getName()=="PUT" || valtype==average;
  if( ss.size() && withderiv ) storedata=true;
  setShape( ss );
}

void Value::setValType( const std::string& vtype ) {
  if( vtype=="normal" ) valtype=normal;
  else if( vtype=="constant" ) valtype=constant;
  else if( vtype=="average" ) valtype=average;
  else if( vtype=="calcFromAverage" ) valtype=calcFromAverage;
  else plumed_merror("invalid valtype " + vtype );
}

void Value::setShape( const std::vector<unsigned>&ss ) {
  std::size_t tot=1; shape.resize( ss.size() );
  for(unsigned i=0; i<shape.size(); ++i) { tot = tot*ss[i]; shape[i]=ss[i]; }

  if( shape.size()>0 && hasDeriv ) {
    // This is for grids
    ngrid_der = shape.size();
    if( action ) ngrid_der = action->getNumberOfDerivatives();
    std::size_t ndata = tot*(1+ngrid_der);
    data.resize( ndata ); inputForce.resize( tot );
  } else if( shape.size()==0 ) {
    // This is for scalars
    data.resize(1); inputForce.resize(1);
  } else if( storedata && shape.size()<2 ) {
    // This is for vectors (matrices have special version because we have sparse storage)
    data.resize( tot ); inputForce.resize( tot );
  }
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

bool Value::isPeriodic()const {
  plumed_massert(periodicity!=unset,"periodicity should be set");
  return periodicity==periodic;
}

bool Value::applyForce(std::vector<double>& forces ) const {
  if( !hasForce || valtype!=normal ) return false;
  plumed_dbg_massert( data.size()-1==forces.size()," forces array has wrong size" );
  const unsigned N=data.size()-1;
  for(unsigned i=0; i<N; ++i) forces[i]=inputForce[0]*data[1+i];
  return true;
}

void Value::setNotPeriodic() {
  min=0; max=0; periodicity=notperiodic;
}

void Value::setDomain(const std::string& pmin,const std::string& pmax) {
  str_min=pmin;
  if( !Tools::convertNoexcept(str_min,min) ) action->error("could not convert period string " + str_min + " to real");
  str_max=pmax;
  if( !Tools::convertNoexcept(str_max,max) ) action->error("could not convert period string " + str_max + " to read");
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

void Value::set(const std::size_t& n, const double& v ) {
  value_set=true;
  if( getRank()==0 ) { plumed_assert( n==0 ); data[n]=v; applyPeriodicity(n); }
  else if( !hasDeriv ) { plumed_dbg_massert( n<data.size(), "failing in " + getName() ); data[n]=v; applyPeriodicity(n); }
  else { data[n*(1+ngrid_der)] = v; }
}

void Value::push_back( const double& v ) {
  value_set=true;
  if( shape.size()==1 ) {
    data.push_back(v); shape[0]++;
  } else if( shape.size()==2 ) {
    data.push_back(v);
    shape[0] = std::ceil( data.size() / shape[1] );
  }
}

double Value::get(const std::size_t& ival, const bool trueind) const {
  if( hasDeriv ) return data[ival*(1+ngrid_der)];
#ifdef DNDEBUG
  if( action ) plumed_dbg_massert( ival<getNumberOfValues(), "could not get value from " + name );
#endif
  if( shape.size()==2 && ncols<shape[1] && trueind ) {
    unsigned irow = std::floor( ival / shape[1] ), jcol = ival%shape[1];
    // This is a special treatment for the lower triangular matrices that are used when
    // we do ITRE with COLLECT_FRAMES
    if( ncols==0 ) {
      if( jcol<=irow ) return data[0.5*irow*(irow+1) + jcol];
      return 0;
    }
    for(unsigned i=0; i<getRowLength(irow); ++i) {
      if( getRowIndex(irow,i)==jcol ) return data[irow*ncols+i];
    }
    return 0.0;
  }
  plumed_massert( ival<data.size(), "cannot get value from " + name );
  return data[ival];
}

void Value::addForce(const std::size_t& iforce, double f, const bool trueind) {
  hasForce=true;
  if( shape.size()==2 && !hasDeriv && ncols<shape[1] && trueind ) {
    unsigned irow = std::floor( iforce / shape[0] ), jcol = iforce%shape[0];
    for(unsigned i=0; i<getRowLength(irow); ++i) {
      if( getRowIndex(irow,i)==jcol ) { inputForce[irow*ncols+i]+=f; return; }
    }
    plumed_assert( fabs(f)<epsilon ); return;
  }
  plumed_massert( iforce<inputForce.size(), "can't add force to " + name );
  inputForce[iforce]+=f;
}


void Value::buildDataStore( const bool forprint ) {
  if( getRank()==0 ) return;
  storedata=true; setShape( shape );
  if( !forprint ) return ;
  ActionWithVector* av=dynamic_cast<ActionWithVector*>( action );
  if( av ) (av->getFirstActionInChain())->never_reduce_tasks=true;
}

void Value::reshapeMatrixStore( const unsigned& n ) {
  plumed_dbg_assert( shape.size()==2 && !hasDeriv );
  if( !storedata ) return ;
  ncols=n; if( ncols>shape[1] ) ncols=shape[1];
  unsigned size=shape[0]*ncols;
  if( matrix_bookeeping.size()!=(size+shape[0]) ) {
    data.resize( size ); inputForce.resize( size );
    matrix_bookeeping.resize( size + shape[0], 0 );
    if( ncols>=shape[1] ) {
      for(unsigned i=0; i<shape[0]; ++i) {
        matrix_bookeeping[(1+ncols)*i] = shape[1];
        for(unsigned j=0; j<shape[1]; ++j) matrix_bookeeping[(1+ncols)*i+1+j]=j;
      }
    }
  }
  if( ncols<shape[1] ) std::fill(matrix_bookeeping.begin(), matrix_bookeeping.end(), 0);
}

void Value::setPositionInMatrixStash( const unsigned& p ) {
  plumed_dbg_assert( shape.size()==2 && !hasDeriv );
  matpos=p;
}

bool Value::ignoreStoredValue(const std::string& c) const {
  if( !storedata && shape.size()>0 ) return true;
  ActionWithVector* av=dynamic_cast<ActionWithVector*>(action);
  if( av ) return (av->getFirstActionInChain())->getLabel()==c;
  return false;
}

void Value::setConstant() {
  valtype=constant; storedata=true; setShape( shape );
  if( getRank()==2 && !hasDeriv ) reshapeMatrixStore( shape[1] );
}

void Value::writeBinary(std::ostream&o) const {
  o.write(reinterpret_cast<const char*>(&data[0]),data.size()*sizeof(double));
}

void Value::setSymmetric( const bool& sym ) {
  plumed_assert( shape.size()==2 && !hasDeriv );
  if( sym && shape[0]!=shape[1] ) plumed_merror("non-square matrix cannot be symmetric");
  symmetric=sym;
}

void Value::retrieveEdgeList( unsigned& nedge, std::vector<std::pair<unsigned,unsigned> >& active, std::vector<double>& elems ) {
  nedge=0; plumed_dbg_assert( shape.size()==2 && !hasDeriv );
  // Check we have enough space to store the edge list
  if( elems.size()<shape[0]*ncols ) { elems.resize( shape[0]*ncols ); active.resize( shape[0]*ncols ); }

  for(unsigned i=0; i<shape[0]; ++i) {
    unsigned ncol = getRowLength(i);
    for(unsigned j=0; j<ncol; ++j) {
      if( fabs(get(i*ncols+j,false))<epsilon ) continue;
      if( symmetric && getRowIndex(i,j)>i ) continue;
      active[nedge].first = i; active[nedge].second = getRowIndex(i,j);
      elems[nedge] = get(i*ncols+j,false); nedge++;
    }
  }
}

void Value::readBinary(std::istream&i) {
  i.read(reinterpret_cast<char*>(&data[0]),data.size()*sizeof(double));
}

void Value::convertIndexToindices(const std::size_t& index, std::vector<unsigned>& indices ) const {
  if( hasDeriv || getRank()==1 ) {
    std::size_t kk=index; indices[0]=index%shape[0];
    for(unsigned i=1; i<shape.size()-1; ++i) {
      kk=(kk-indices[i-1])/shape[i-1];
      indices[i]=kk%shape[i];
    }
    if(shape.size()>=2) indices[shape.size()-1]=(kk-indices[shape.size()-2])/shape[shape.size()-2];
  } else if( getRank()==2 ) {
    indices[0]=std::floor( index/shape[1] ); indices[1] = index%shape[1];
  }
}

void Value::print( OFile& ofile ) const {
  if( isPeriodic() ) { ofile.printField( "min_" + name, str_min ); ofile.printField("max_" + name, str_max ); }
  if( shape.size()==0 || getNumberOfValues()==1 ) {
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

unsigned Value::getGoodNumThreads( const unsigned& j, const unsigned& k ) const {
  return OpenMP::getGoodNumThreads( &data[j], (k-j) );
}

bool Value::calculateOnUpdate() const {
  return (valtype==average || valtype==calcFromAverage);
}

std::string Value::getValueType() const {
  if( getRank()==0 ) return "scalar";
  if( getRank()>0 && hasDerivatives() ) return "grid";
  if( getRank()==1 ) return "vector";
  if( getRank()==2 ) return "matrix";
  plumed_merror("unknown type for value " + getName() );
  return "";
}

}
