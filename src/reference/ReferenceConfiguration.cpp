/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2017 The plumed team
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
#include "ReferenceConfiguration.h"
#include "ReferenceArguments.h"
#include "ReferenceAtoms.h"
#include "Direction.h"
#include "core/Value.h"
#include "tools/OFile.h"
#include "tools/PDB.h"

namespace PLMD {

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
// arg_ders(0),
// atom_ders(0)
{
  weight=0.0;
}

ReferenceConfiguration::~ReferenceConfiguration()
{
}

std::string ReferenceConfiguration::getName() const {
  return name;
}

void ReferenceConfiguration::set( const PDB& pdb ) {
  line=pdb.getRemark();
  std::string ignore;
  if( parse("TYPE",ignore,true) ) {
    if(ignore!=name) error("mismatch for name");
  }
  if( !parse("WEIGHT",weight,true) ) weight=1.0;
  read( pdb );
}

// void ReferenceConfiguration::setNumberOfArguments( const unsigned& n ){
//   arg_ders.resize(n); tmparg.resize(n);
// }

// void ReferenceConfiguration::setNumberOfAtoms( const unsigned& n ){
//   atom_ders.resize(n);
// }

// bool ReferenceConfiguration::getVirial( Tensor& virout ) const {
//   if(virialWasSet) virout=virial;
//   return virialWasSet;
// }

void ReferenceConfiguration::parseFlag( const std::string&key, bool&t ) {
  Tools::parseFlag(line,key,t);
}

void ReferenceConfiguration::error(const std::string& msg) {
  plumed_merror("error reading reference configuration of type " + name + " : " + msg );
}

void ReferenceConfiguration::setNamesAndAtomNumbers( const std::vector<AtomNumber>& numbers, const std::vector<std::string>& arg ) {
  ReferenceAtoms* atoms=dynamic_cast<ReferenceAtoms*>( this );
  if(!atoms) {
    plumed_massert( numbers.size()==0, "expecting no atomic positions");
    //setNumberOfAtoms( 0 );
  } else {
    atoms->setAtomNumbers( numbers );
    // setNumberOfAtoms( numbers.size() );
  }
  // Copy the arguments to the reference
  ReferenceArguments* args=dynamic_cast<ReferenceArguments*>( this );
  if(!args) {
    plumed_massert( arg.size()==0, "expecting no arguments");
    // setNumberOfArguments(0);
  } else {
    args->setArgumentNames( arg );
    // setNumberOfArguments( arg.size() );
  }
}

void ReferenceConfiguration::setReferenceConfig( const std::vector<Vector>& pos, const std::vector<double>& arg, const std::vector<double>& metric ) {
//  plumed_dbg_assert( pos.size()==atom_ders.size() && arg.size()==arg_ders.size() );
  // Copy the atomic positions to the reference
  ReferenceAtoms* atoms=dynamic_cast<ReferenceAtoms*>( this );
  if(!atoms) {
    plumed_massert( pos.size()==0, "expecting no atomic positions");
  } else {
    std::vector<double> align_in( pos.size(), 1.0 ), displace_in( pos.size(), 1.0 );
    atoms->setReferenceAtoms( pos, align_in, displace_in );
  }
  // Copy the arguments to the reference
  ReferenceArguments* args=dynamic_cast<ReferenceArguments*>( this );
  if(!args) {
    plumed_massert( arg.size()==0 && metric.size()==0, "expecting no arguments");
  } else {
    args->setReferenceArguments( arg, metric );
  }
}

void ReferenceConfiguration::checkRead() {
  if(!line.empty()) {
    std::string msg="cannot understand the following words from the input line : ";
    for(unsigned i=0; i<line.size(); i++) {
      if(i>0) msg = msg + ", ";
      msg = msg + line[i];
    }
    error(msg);
  }
}

double ReferenceConfiguration::calculate( const std::vector<Vector>& pos, const Pbc& pbc, const std::vector<Value*>& vals,
    ReferenceValuePack& myder, const bool& squared ) const {
  // clearDerivatives();
  std::vector<double> tmparg( vals.size() );
  for(unsigned i=0; i<vals.size(); ++i) tmparg[i]=vals[i]->get();
  return calc( pos, pbc, vals, tmparg, myder, squared );
}

// void ReferenceConfiguration::copyDerivatives( const ReferenceConfiguration* ref ){
//   plumed_dbg_assert( ref->atom_ders.size()==atom_ders.size() && ref->arg_ders.size()==arg_ders.size() );
//   for(unsigned i=0;i<atom_ders.size();++i) atom_ders[i]=ref->atom_ders[i];
//   for(unsigned i=0;i<arg_ders.size();++i) arg_ders[i]=ref->arg_ders[i];
//   virialWasSet=ref->virialWasSet; virial=ref->virial;
// }

void ReferenceConfiguration::print( OFile& ofile, const double& time, const double& weight, const double& lunits, const double& old_norm ) {
  ofile.printf("REMARK TIME=%f LOG_WEIGHT=%f OLD_NORM=%f\n",time, weight, old_norm );
  print( ofile, "%f", lunits );  // HARD CODED FORMAT HERE AS THIS IS FOR CHECKPOINT FILE
}

void ReferenceConfiguration::print( OFile& ofile, const std::string& fmt, const double& lunits ) {
  ofile.printf("REMARK TYPE=%s\n",getName().c_str() );
  ReferenceArguments* args=dynamic_cast<ReferenceArguments*>(this);
  if(args) args->printArguments( ofile, fmt );
  ReferenceAtoms* atoms=dynamic_cast<ReferenceAtoms*>(this);
  if(atoms) atoms->printAtoms( ofile, lunits );
  ofile.printf("END\n");
}

void ReferenceConfiguration::displaceReferenceConfiguration( const double& weight, Direction& dir ) {
  ReferenceArguments* args=dynamic_cast<ReferenceArguments*>(this);
  if( args ) args->displaceReferenceArguments( weight, dir.getReferenceArguments() );
  ReferenceAtoms* atoms=dynamic_cast<ReferenceAtoms*>(this);
  if( atoms ) atoms->displaceReferenceAtoms( weight, dir.getReferencePositions() );
}

void ReferenceConfiguration::extractDisplacementVector( const std::vector<Vector>& pos, const std::vector<Value*>& vals,
    const std::vector<double>& arg, const bool& nflag,
    Direction& mydir ) const {
  const ReferenceAtoms* atoms=dynamic_cast<const ReferenceAtoms*>( this );
  if( atoms ) atoms->extractAtomicDisplacement( pos, mydir.reference_atoms );
  const ReferenceArguments* args=dynamic_cast<const ReferenceArguments*>( this );
  if( args ) args->extractArgumentDisplacement( vals, arg, mydir.reference_args );

  // Normalize direction if required
  if( nflag ) {
    // Calculate length of vector
    double tmp, norm=0; mydir.normalized = true;
    for(unsigned i=0; i<mydir.getReferencePositions().size(); ++i) {
      for(unsigned k=0; k<3; ++k) { tmp=mydir.getReferencePositions()[i][k]; norm+=tmp*tmp; }
    }
    for(unsigned i=0; i<mydir.getReferenceArguments().size(); ++i) { tmp=mydir.getReferenceArguments()[i]; norm+=tmp*tmp; }
    norm = sqrt( norm );
    // And normalize
    for(unsigned i=0; i<mydir.getReferencePositions().size(); ++i) {
      for(unsigned k=0; k<3; ++k) { mydir.reference_atoms[i][k] /=norm; }
    }
    for(unsigned i=0; i<mydir.getReferenceArguments().size(); ++i) { mydir.reference_args[i] /= norm; }
  }
}

double ReferenceConfiguration::projectDisplacementOnVector( const Direction& mydir,
    const std::vector<Value*>& vals, const std::vector<double>& arg,
    ReferenceValuePack& mypack ) const {
  double proj=0;
  const ReferenceAtoms* atoms=dynamic_cast<const ReferenceAtoms*>( this );
  if( atoms ) proj += atoms->projectAtomicDisplacementOnVector( mydir.normalized, mydir.getReferencePositions(), mypack );
  const ReferenceArguments* args=dynamic_cast<const ReferenceArguments*>( this );
  if( args ) proj += args->projectArgDisplacementOnVector( mydir.getReferenceArguments(), vals, arg, mypack );
  return proj;
}

double distance( const Pbc& pbc, const std::vector<Value*> & vals, ReferenceConfiguration* ref1, ReferenceConfiguration* ref2, const bool& squared ) {
  unsigned nder;
  if( ref1->getReferencePositions().size()>0 ) nder=ref1->getReferenceArguments().size() + 3*ref1->getReferencePositions().size() + 9;
  else nder=ref1->getReferenceArguments().size();

  MultiValue myvals( 1, nder ); ReferenceValuePack myder( ref1->getReferenceArguments().size(), ref1->getReferencePositions().size(), myvals );
  double dist1=ref1->calc( ref2->getReferencePositions(), pbc, vals, ref2->getReferenceArguments(), myder, squared );
#ifndef NDEBUG
  // Check that A - B = B - A
  double dist2=ref2->calc( ref1->getReferencePositions(), pbc, vals, ref1->getReferenceArguments(), myder, squared );
  plumed_dbg_assert( fabs(dist1-dist2)<epsilon );
#endif
  return dist1;
}

}
