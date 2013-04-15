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
#include "ReferenceConfiguration.h"
#include "ReferenceArguments.h"
#include "ReferenceAtoms.h"
#include "tools/OFile.h"
#include "tools/PDB.h"

namespace PLMD{

ReferenceConfigurationOptions::ReferenceConfigurationOptions( const std::string& type, Log& log ):
tt(type),
ll(&log)
{
}

bool ReferenceConfigurationOptions::usingFastOption() const {
  return (tt.find("-FAST")!=std::string::npos);
}

ReferenceConfiguration::ReferenceConfiguration( const ReferenceConfigurationOptions& ro ):
name(ro.tt),
log(ro.ll),
arg_ders(0),
atom_ders(0)
{
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
  read( pdb );
}

void ReferenceConfiguration::setNumberOfArguments( const unsigned& n ){
  arg_ders.resize(n);
}

void ReferenceConfiguration::setNumberOfAtoms( const unsigned& n ){
  atom_ders.resize(n);
}

bool ReferenceConfiguration::getVirial( Tensor& virout ) const {
  if(virialWasSet) virout=virial;
  return virialWasSet;
}

void ReferenceConfiguration::parseFlag( const std::string&key, bool&t ){
  Tools::parseFlag(line,key,t);
}

void ReferenceConfiguration::error(const std::string& msg){
  plumed_merror("error reading reference configuration of type " + name + " : " + msg );
}

void ReferenceConfiguration::setNamesAndAtomNumbers( const std::vector<AtomNumber>& numbers, const std::vector<Value*>& arg ){
   ReferenceAtoms* atoms=dynamic_cast<ReferenceAtoms*>( this );
  if(!atoms){
     plumed_massert( numbers.size()==0, "expecting no atomic positions");
     setNumberOfAtoms( 0 );
  } else { 
     atoms->setAtomNumbers( numbers );
     setNumberOfAtoms( numbers.size() );
  }
  // Copy the arguments to the reference
  ReferenceArguments* args=dynamic_cast<ReferenceArguments*>( this );
  if(!args){
     plumed_massert( arg.size()==0, "expecting no arguments");
     setNumberOfArguments(0);
  } else {
     args->setArgumentNames( arg );
     setNumberOfArguments( arg.size() );
  }
}

void ReferenceConfiguration::setReference( const std::vector<Vector>& pos, const std::vector<Value*>& arg, const std::vector<double>& metric ){
  plumed_dbg_assert( pos.size()==atomm_ders.size() && arg.size()==arg_ders.size() );
  // Copy the atomic positions to the reference
  ReferenceAtoms* atoms=dynamic_cast<ReferenceAtoms*>( this );
  if(!atoms){
     plumed_massert( pos.size()==0, "expecting no atomic positions");
  } else { 
     std::vector<double> align_in( pos.size(), 1.0 ), displace_in( pos.size(), 1.0 );
     atoms->setReferenceAtoms( pos, align_in, displace_in );
  }
  // Copy the arguments to the reference
  ReferenceArguments* args=dynamic_cast<ReferenceArguments*>( this );
  if(!args){
     plumed_massert( arg.size()==0 && metric.size()==0, "expecting no arguments");
  } else {
     args->setReferenceArguments( arg, metric );
  }
}

void ReferenceConfiguration::checkRead(){
  if(!line.empty()){
    std::string msg="cannot understand the following words from the input line : ";
    for(unsigned i=0;i<line.size();i++) msg = msg + line[i] + ", ";
    error(msg);
  }
}

void ReferenceConfiguration::clearDerivatives(){
  for(unsigned i=0;i<atom_ders.size();++i) atom_ders[i].zero();
  virial.zero(); virialWasSet=false;
  arg_ders.assign(arg_ders.size(),0.0);
}

double ReferenceConfiguration::calculate( const std::vector<Vector>& pos, const Pbc& pbc, const std::vector<Value*>& vals, const bool& squared ){
  clearDerivatives();
  return calc( pos, pbc, vals, squared );
}

void ReferenceConfiguration::copyDerivatives( const ReferenceConfiguration* ref ){
  plumed_dbg_assert( ref->atom_ders.size()==atom_ders.size() && ref->arg_ders.size()==arg_ders.size() );
  for(unsigned i=0;i<atom_ders.size();++i) atom_ders[i]=ref->atom_ders[i];
  for(unsigned i=0;i<arg_ders.size();++i) arg_ders[i]=ref->arg_ders[i];
  virialWasSet=ref->virialWasSet; virial=ref->virial;
}

void ReferenceConfiguration::print( OFile& ofile, const double& time, const double& weight, const double& old_norm ){
  ofile.printf("REMARK: TIME=%f LOG_WEIGHT=%f OLD_NORM=%f\n",time, weight, old_norm );
  ReferenceArguments* args=dynamic_cast<ReferenceArguments*>(this);
  if(args) args->printArguments( ofile );
  ReferenceAtoms* atoms=dynamic_cast<ReferenceAtoms*>(this);
  if(atoms) atoms->printAtoms( ofile );
  ofile.printf("END\n");
}

}
