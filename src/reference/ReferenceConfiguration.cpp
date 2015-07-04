/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2015 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

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
#include "ReferenceConfiguration.h"
#include "ReferenceArguments.h"
#include "ReferenceAtoms.h"
#include "core/Value.h"
#include "tools/OFile.h"
#include "tools/PDB.h"

namespace PLMD{

ReferenceConfigurationOptions::ReferenceConfigurationOptions( const std::string& type ):
tt(type)
{
}

bool ReferenceConfigurationOptions::usingFastOption() const {
  return (tt.find("-FAST")!=std::string::npos);
}

std::string ReferenceConfigurationOptions::getMultiRMSDType() const {
  plumed_assert( tt.find("MULTI-")!=std::string::npos );
  std::size_t dot=tt.find_first_of("MULTI-");
  return tt.substr(dot+6); 
}

ReferenceConfiguration::ReferenceConfiguration( const ReferenceConfigurationOptions& ro ):
name(ro.tt)
{
  weight=0.0;
}

ReferenceConfiguration::~ReferenceConfiguration()
{
}

std::string ReferenceConfiguration::getName() const {
  return name;
}

void ReferenceConfiguration::set( const PDB& pdb ){
  line=pdb.getRemark(); 
  std::string ignore; 
  if( parse("TYPE",ignore,true) ){
      if(ignore!=name) error("mismatch for name");
  }
  if( !parse("WEIGHT",weight,true) ) weight=1.0;
  // Read in properties
  parseVector( "PROPERTIES", property_names, true );
  property_values.resize( property_names.size() );
  for(unsigned i=0;i<property_names.size();++i) parse( property_names[i], property_values[i] );
  // And read in rest of pdb
  read( pdb );
}

void ReferenceConfiguration::parseFlag( const std::string&key, bool&t ){
  Tools::parseFlag(line,key,t);
}

void ReferenceConfiguration::error(const std::string& msg){
  plumed_merror("error reading reference configuration of type " + name + " : " + msg );
}

void ReferenceConfiguration::checkRead(){
  if(!line.empty()){
    std::string msg="cannot understand the following words from the input line : ";
    for(unsigned i=0;i<line.size();i++) msg = msg + line[i] + ", ";
    error(msg);
  }
}

bool ReferenceConfiguration::isDirection() const {
  return ( name=="DIRECTION" );
}

double ReferenceConfiguration::calculate( const std::vector<Vector>& pos, const Pbc& pbc, const std::vector<Value*>& vals, 
                                          ReferenceValuePack& myder, const bool& squared ) const {
  std::vector<double> tmparg( vals.size() );
  for(unsigned i=0;i<vals.size();++i) tmparg[i]=vals[i]->get();
  return calc( pos, pbc, vals, tmparg, myder, squared );
}

void ReferenceConfiguration::print( const double& lunits, OFile& ofile, const double& time, const double& weight, const double& old_norm ){
  ofile.printf("REMARK TIME=%f LOG_WEIGHT=%f OLD_NORM=%f\n",time, weight, old_norm );
  print( lunits, ofile, "%f" );  // HARD CODED FORMAT HERE AS THIS IS FOR CHECKPOINT FILE
}

void ReferenceConfiguration::print( const double& lunits, OFile& ofile, const std::string& fmt ){
  ReferenceArguments* args=dynamic_cast<ReferenceArguments*>(this);
  if( property_names.size()>0 ){
      ofile.printf("REMARK PROPERTIES=%s", property_names[0].c_str() ); 
      for(unsigned i=1;i<property_names.size();++i) ofile.printf(",%s", property_names[i].c_str() );
      ofile.printf("\nREMARK "); std::string descr2;
      if(fmt.find("-")!=std::string::npos){
         descr2="%s=" + fmt + " ";
      } else {
         // This ensures numbers are left justified (i.e. next to the equals sign
         std::size_t psign=fmt.find("%");
         plumed_assert( psign!=std::string::npos );
         descr2="%s=%-" + fmt.substr(psign+1) + " ";
      }
      for(unsigned i=0;i<property_names.size();++i) ofile.printf( descr2.c_str(),property_names[i].c_str(), property_values[i] );
      ofile.printf("\n");
  }
  if(args) args->printArguments( ofile, fmt );
  ReferenceAtoms* atoms=dynamic_cast<ReferenceAtoms*>(this);
  if(atoms) atoms->printAtoms( lunits, ofile );
  ofile.printf("END\n");
}

void ReferenceConfiguration::clearAllProperties(){
  property_names.resize(0); property_values.resize(0);
}

double ReferenceConfiguration::getPropertyValue( const std::string& myname ) const {
  bool found=false;
  for(unsigned i=0;i<property_names.size();++i){
      if( myname==property_names[i] ) return property_values[i];
  }
  plumed_assert( false ); return 0.0;
}

void ReferenceConfiguration::attachProperty( const std::string& name, const double& val ){
  bool found=false;
  for(unsigned i=0;i<property_names.size();++i){
      if( property_names[i]==name ){ found=false; property_values[i]=val; break; }
  }
  plumed_dbg_assert( property_names.size()==property_values.size() );
  if( !found ){ property_names.push_back( name ); property_values.push_back( val ); }
}

double distance( const Pbc& pbc, const std::vector<Value*> & vals, ReferenceConfiguration* ref1, ReferenceConfiguration* ref2, const bool& squared ){
  unsigned nder;
  if( ref1->getReferencePositions().size()>0 ) nder=ref1->getReferenceArguments().size() + 3*ref1->getReferencePositions().size() + 9;
  else nder=ref1->getReferenceArguments().size();

  MultiValue myvals( 1, nder ); ReferenceValuePack myder( ref1->getReferenceArguments().size() , ref1->getReferencePositions().size() , myvals ); 
  double dist1=ref1->calc( ref2->getReferencePositions(), pbc, vals, ref2->getReferenceArguments(), myder, squared );
#ifndef NDEBUG
  // Check that A - B = B - A
  double dist2=ref2->calc( ref1->getReferencePositions(), pbc, vals, ref1->getReferenceArguments(), myder, squared );
  plumed_dbg_assert( fabs(dist1-dist2)<epsilon );
#endif 
  return dist1;
}


}
