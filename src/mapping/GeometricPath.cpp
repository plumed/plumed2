/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016,2017 The plumed team
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
#include "core/ActionWithValue.h"
#include "core/ActionWithArguments.h"
#include "core/ActionRegister.h"
#include "reference/ReferenceValuePack.h"
#include "reference/ReferenceConfiguration.h"
#include "reference/Direction.h"
#include "tools/OpenMP.h"
#include "tools/PDB.h"
#include "Mapping.h"

using namespace std;

namespace PLMD {
namespace mapping {

//+PLUMEDOC FUNCTION GPATH
/*
Calculate path collective variable given a set of distances from a collection of waymarkers.

\par Examples

*/
//+ENDPLUMEDOC


class GeometricPath : public ActionWithValue, public ActionWithArguments {
private:
  Direction projdir;
  std::vector<double> framep;
  std::vector<double> forcesToApply;
  MultiValue mydpack1, mydpack2, mydpack3;
  ReferenceValuePack mypack1, mypack2, mypack3;
  std::vector<double> mypack1_stashd_args;
  std::vector<Vector> mypack1_stashd_atoms;
  ReferenceConfiguration* getReferenceConfiguration( const unsigned& iframe );
  void retrieveReferenceValuePackData( const unsigned& iframe, ReferenceValuePack& mypack );
  void extractDisplacementVector( const unsigned & iframe, const std::vector<Vector>& pos, const std::vector<double>& args, Direction& mydir );
  double calculateDistanceFromPosition( const unsigned & iframe, const std::vector<Vector>& pos, const std::vector<double>& args, ReferenceValuePack& mypack );
  double projectDisplacementOnVector( const unsigned & iframe, const Direction& dir, ReferenceValuePack& mypack );
public:
  static void shortcutKeywords( Keywords& keys );
  static void expandShortcut( const std::string& lab, const std::vector<std::string>& words,
                              const std::map<std::string,std::string>& keys,
                              std::vector<std::vector<std::string> >& actions );
  static void registerKeywords(Keywords& keys);
  explicit GeometricPath(const ActionOptions&);
  void calculate();
  void calculateFunction( const std::vector<double>& args, MultiValue& myvals ) const { plumed_error(); }
  unsigned getNumberOfDerivatives();
  void apply();
};

PLUMED_REGISTER_ACTION(GeometricPath,"GPATH")
PLUMED_REGISTER_SHORTCUT(GeometricPath,"GPATH")

void GeometricPath::shortcutKeywords( Keywords& keys ){
  keys.add("compulsory","REFERENCE","a pdb file containing the set of reference configurations");
  keys.add("compulsory","TYPE","OPTIMAL-FAST","the manner in which distances are calculated. More information on the different "
           "metrics that are available in PLUMED can be found in the section of the manual on "
           "\\ref dists");
}

void GeometricPath::expandShortcut( const std::string& lab, const std::vector<std::string>& words,
                                    const std::map<std::string,std::string>& keys,
                                    std::vector<std::vector<std::string> >& actions ){
  std::vector<std::string> ref_line; ref_line.push_back( lab + "_data:" );
  ref_line.push_back("EUCLIDEAN_DISSIMILARITIES_VECTOR");
  for(const auto & p : keys ) ref_line.push_back( p.first + "=" + p.second );
  ref_line.push_back("SQUARED"); actions.push_back( ref_line );
  std::vector<std::string> path_line; unsigned nfram = 0;
  path_line.push_back( lab + ":" );
  for(unsigned i=0;i<words.size();++i) path_line.push_back(words[i]);
  path_line.push_back("ARG=" + lab + "_data" ); actions.push_back( path_line );
}

void GeometricPath::registerKeywords(Keywords& keys) {
  Action::registerKeywords(keys); ActionWithValue::registerKeywords(keys); ActionWithArguments::registerKeywords(keys); keys.use("ARG"); 
  keys.add("optional","COORDINATES","a vector of coordinates describing the position of each point along the path.  The default "
                                    "is to place these coordinates at 1, 2, 3, ...");
  componentsAreNotOptional(keys);
  keys.addOutputComponent("s","default","the position on the path");
  keys.addOutputComponent("z","default","the distance from the path");
}

GeometricPath::GeometricPath(const ActionOptions&ao):
  Action(ao),
  ActionWithValue(ao),
  ActionWithArguments(ao),
  projdir(ReferenceConfigurationOptions("DIRECTION")),
  mydpack1( 1, 0 ),
  mydpack2( 1, 0 ),
  mydpack3( 1, 0 ),
  mypack1( 0, 0, mydpack1 ),
  mypack2( 0, 0, mydpack2 ),
  mypack3( 0, 0, mydpack3 )
{
  // Ensure that values are stored in base calculation and that PLUMED doesn't try to calculate this in the stream
  done_over_stream = false;  getPntrToArgument(0)->buildDataStore();
  // Some sanity checks
  bool rankOneOutput = getPntrToArgument(0)->getRank()>0;
  if( getPntrToArgument(0)->getRank()>1 ) error("input arguments should be rank 0 or rank 1");
  if( rankOneOutput && getNumberOfArguments()>1 ) error("cannot calculate path for more than one one vector at a time");
  if( arg_ends[1]-arg_ends[0]!=1 ) error("makes no sense to use ARG1, ARG2... with this action use single ARG keyword");

  unsigned maxargs=0, maxatoms=0; std::vector<AtomNumber> myatoms; std::vector<std::string> myargs;
  for(unsigned i=0;i<getNumberOfArguments();++i){
     if( getPntrToArgument(i)->isPeriodic() ) error("cannot use this function on periodic functions");
     // Now get number of atoms
     ActionAtomistic* aa = dynamic_cast<ActionAtomistic*>( getPntrToArgument(i)->getPntrToAction() );
     if( aa->getNumberOfAtoms()>maxatoms ){ maxatoms = aa->getNumberOfAtoms(); myatoms = aa->getAbsoluteIndexes(); }
     // And number of arguments
     ActionWithArguments* aw = dynamic_cast<ActionWithArguments*>( getPntrToArgument(i)->getPntrToAction() );
     if( aw->getNumberOfArguments()>maxargs ){ 
         maxargs = aw->getNumberOfArguments(); 
         myargs.resize( maxargs ); for(unsigned i=0;i<maxargs;++i) myargs[i] = (aw->getPntrToArgument(i))->getName();  
     }
  }
  parseVector("COORDINATES",framep);
  if( framep.size()>0 ){
     if( framep.size()!=getNumberOfArguments() ){ 
        if( framep.size()!=getPntrToArgument(0)->getShape()[0] ) error("wrong number of input coordinates");
     }
  } else {                    
     if( getNumberOfArguments()==1 ) framep.resize( getPntrToArgument(0)->getShape()[0] );
     else framep.resize( getNumberOfArguments() );
     for(unsigned i=0;i<framep.size();++i) framep[i] = static_cast<double>(i+1);
  } 
  log.printf("  coordinates of points on path : ");
  for(unsigned i=0;i<framep.size();++i) log.printf("%f ",framep[i] );
  log.printf("\n");
  projdir.setNamesAndAtomNumbers( myatoms, myargs );
  mypack1_stashd_atoms.resize( maxatoms ); mypack1_stashd_args.resize( maxargs );
  Mapping* am = dynamic_cast<Mapping*>( getPntrToArgument(0)->getPntrToAction() ); 
  unsigned maxderiv=maxargs; if( maxatoms>0 ) maxderiv += 3*maxatoms + 9;
  mydpack1.resize( 1, maxderiv ); mypack1.resize( maxargs, maxatoms ); (am->getReferenceConfiguration(0))->setupPCAStorage( mypack1 );
  mydpack2.resize( 1, maxderiv ); mypack2.resize( maxargs, maxatoms ); (am->getReferenceConfiguration(0))->setupPCAStorage( mypack2 );
  mydpack3.resize( 1, maxderiv ); mypack3.resize( maxargs, maxatoms ); forcesToApply.resize( maxderiv );
  for(unsigned i=0; i<maxatoms; ++i) { mypack1.setAtomIndex(i,i); mypack2.setAtomIndex(i,i); mypack3.setAtomIndex(i,i); }

  checkRead(); 
  addComponentWithDerivatives("s"); componentIsNotPeriodic("s");
  addComponentWithDerivatives("z"); componentIsNotPeriodic("z");
}

unsigned GeometricPath::getNumberOfDerivatives(){
  unsigned nder=0;
  if( mypack1_stashd_atoms.size()>0 ) nder = 3*mypack1_stashd_atoms.size() + 9;
  return nder + mypack1_stashd_args.size();
}

void GeometricPath::retrieveReferenceValuePackData( const unsigned& iframe, ReferenceValuePack& mypack ) {
  mypack.clear(); 
  if( getNumberOfArguments()>1 ){
      Mapping* am = dynamic_cast<Mapping*>( getPntrToArgument(iframe)->getPntrToAction() );
      double ig=am->calculateDistanceFromReference( 0, mypack );
  }
  Mapping* am = dynamic_cast<Mapping*>( getPntrToArgument(0)->getPntrToAction() );
  double ig=am->calculateDistanceFromReference( iframe, mypack );
}

ReferenceConfiguration* GeometricPath::getReferenceConfiguration( const unsigned& iframe ) {
  if( getNumberOfArguments()>1 ){
      Mapping* am = dynamic_cast<Mapping*>( getPntrToArgument(iframe)->getPntrToAction() );      
      return am->getReferenceConfiguration( iframe );
  }
  Mapping* am = dynamic_cast<Mapping*>( getPntrToArgument(0)->getPntrToAction() );
  return am->getReferenceConfiguration( iframe );
}

double GeometricPath::calculateDistanceFromPosition( const unsigned & iframe, const std::vector<Vector>& pos, const std::vector<double>& args, ReferenceValuePack& mypack ){ 
  mypack.clear();
  if( getNumberOfArguments()>1 ){
      Mapping* am = dynamic_cast<Mapping*>( getPntrToArgument(iframe)->getPntrToAction() );
      return am->calculateDistanceBetweenReferenceAndThisPoint( 0, pos, args, mypack );
  }
  Mapping* am = dynamic_cast<Mapping*>( getPntrToArgument(0)->getPntrToAction() );
  return am->calculateDistanceBetweenReferenceAndThisPoint( iframe, pos, args, mypack );
}

double GeometricPath::projectDisplacementOnVector( const unsigned & iframe, const Direction& dir, ReferenceValuePack& mypack ) {
  if( getNumberOfArguments()>1 ){
      Mapping* am = dynamic_cast<Mapping*>( getPntrToArgument(iframe)->getPntrToAction() );
      return am->projectDisplacementOnVector( 0, dir, mypack );
  } 
  Mapping* am = dynamic_cast<Mapping*>( getPntrToArgument(0)->getPntrToAction() );
  return am->projectDisplacementOnVector( iframe, dir, mypack );
}

void GeometricPath::extractDisplacementVector( const unsigned & iframe, const std::vector<Vector>& pos, const std::vector<double>& args, Direction& mydir ){
  if( getNumberOfArguments()>1 ){
      Mapping* am = dynamic_cast<Mapping*>( getPntrToArgument(iframe)->getPntrToAction() );
      am->extractDisplacementVector( 0, pos, args, mydir ); return;
  }
  Mapping* am = dynamic_cast<Mapping*>( getPntrToArgument(0)->getPntrToAction() );
  am->extractDisplacementVector( iframe, pos, args, mydir );
}

void GeometricPath::calculate(){
  // Find two closest points
  unsigned iclose1 = 0, iclose2 = 1;
  double v1v1 = getArgumentScalar(0); 
  double v3v3 = getArgumentScalar(1);
  if( v1v1>v3v3 ){ 
      double tmp=v1v1; v1v1=v3v3; v3v3=tmp; 
      iclose1 = 1; iclose2 = 0;
  }
  for(unsigned i=2; i<getNumberOfScalarArguments(); ++i) {
    double ndist=getArgumentScalar(i);
    if( ndist<v1v1 ) { 
      v3v3=v1v1; iclose2=iclose1;
      v1v1=ndist; iclose1=i;
    } else if( ndist<v3v3 ) {
      v3v3=ndist; iclose2=i;
    }
  }
  // And find third closest point
  int isign = iclose1 - iclose2;
  if( isign>1 ) isign=1; else if( isign<-1 ) isign=-1;
  int iclose3 = iclose1 + isign; double v2v2;

  retrieveReferenceValuePackData( iclose1, mypack1 ); 
  retrieveReferenceValuePackData( iclose2, mypack3 );
  if( iclose3<0 || iclose3>=getNumberOfScalarArguments() ) {
    ReferenceConfiguration* conf2=getReferenceConfiguration( iclose1 );
    v2v2=calculateDistanceFromPosition( iclose2, conf2->getReferencePositions(), conf2->getReferenceArguments(), mypack2 ); 
    extractDisplacementVector( iclose2, conf2->getReferencePositions(), conf2->getReferenceArguments(), projdir );
  } else {
    ReferenceConfiguration* conf2=getReferenceConfiguration( iclose3 );
    v2v2=calculateDistanceFromPosition( iclose1, conf2->getReferencePositions(), conf2->getReferenceArguments(), mypack2 );
    extractDisplacementVector( iclose1, conf2->getReferencePositions(), conf2->getReferenceArguments(), projdir );
  }

  // Stash derivatives of v1v1
  if( !doNotCalculateDerivatives() ){
      for(unsigned i=0; i<mypack1_stashd_args.size(); ++i) mypack1_stashd_args[i]=mypack1.getArgumentDerivative(i);
      for(unsigned i=0; i<mypack1_stashd_atoms.size(); ++i) mypack1_stashd_atoms[i]=mypack1.getAtomDerivative(i);
  }
  if( mypack1_stashd_atoms.size()>0 ) {
      ReferenceAtoms* at = dynamic_cast<ReferenceAtoms*>( getReferenceConfiguration( iclose1 ) );
      const std::vector<double> & displace( at->getDisplace() );
      for(unsigned i=0; i<mypack1_stashd_atoms.size(); ++i) mypack1.getAtomsDisplacementVector()[i] /= displace[i];
  }
  // Calculate the dot product of v1 with v2
  double v1v2 = projectDisplacementOnVector( iclose1, projdir, mypack1 );

  // This computes s value
  double spacing = framep[iclose1] - framep[iclose2];
  double root = sqrt( v1v2*v1v2 - v2v2 * ( v1v1 - v3v3) );
  double dx = 0.5 * ( (root + v1v2) / v2v2 - 1.);
  double path_s = framep[iclose1] + spacing * dx; 
  double fact = 0.25*spacing / v2v2; Value* sp = getPntrToComponent(0); sp->set( path_s );
  if( !doNotCalculateDerivatives() ) {
      // Derivative of s wrt arguments
      for(unsigned i=0; i<mypack1_stashd_args.size(); ++i) {
        sp->setDerivative( i, fact*( mypack2.getArgumentDerivative(i) + (v2v2 * (-mypack1_stashd_args[i] + mypack3.getArgumentDerivative(i))
                                     + v1v2*mypack2.getArgumentDerivative(i) )/root ) );
      }
      // Derivative of s wrt atoms
      if( mypack1_stashd_atoms.size()>0 ) {
        unsigned narg=mypack1_stashd_args.size(); Tensor vir; vir.zero(); fact = 0.5*spacing / v2v2;
        ActionAtomistic* mymap = dynamic_cast<ActionAtomistic*>( getPntrToArgument(0)->getPntrToAction() );
        for(unsigned i=0; i<mypack1_stashd_atoms.size(); ++i) {
          Vector ader = fact*(( v1v2*mypack1.getAtomDerivative(i) + 0.5*v2v2*(-mypack1_stashd_atoms[i] + mypack3.getAtomDerivative(i) ) )/root + mypack1.getAtomDerivative(i) );
          for(unsigned k=0; k<3; ++k) sp->setDerivative( narg+3*i+k, ader[k] );
          vir-=Tensor( mymap->getPosition(i), ader );
        }
        // Set the virial
        unsigned nbase=narg+3*mypack1_stashd_atoms.size();
        for(unsigned i=0; i<3; ++i) for(unsigned j=0; j<3; ++j) sp->setDerivative( nbase+3*i+j, vir(i,j) );
      }
  }

  // Now compute z value
  ReferenceConfiguration* conf2=getReferenceConfiguration( iclose1 );
  double v4v4=calculateDistanceFromPosition( iclose2, conf2->getReferencePositions(), conf2->getReferenceArguments(), mypack2 );
  // Extract vector connecting frames
  extractDisplacementVector( iclose2, conf2->getReferencePositions(), conf2->getReferenceArguments(), projdir ); 
  // Calculate projection of vector on line connnecting frames
  double proj = projectDisplacementOnVector( iclose1, projdir, mypack1 ); double path_z = v1v1 + dx*dx*v4v4 - 2*dx*proj;

  // Derivatives for z path
  path_z = sqrt(path_z); Value* zp = getPntrToComponent(1); zp->set( path_z ); 
  if( !doNotCalculateDerivatives() ) {
      for(unsigned i=0; i<mypack1_stashd_args.size(); ++i) zp->setDerivative( i, (mypack1_stashd_args[i] - 2*dx*mypack1.getArgumentDerivative(i))/(2.0*path_z) );
      // Derivative wrt atoms
      if( mypack1_stashd_atoms.size()>0 ) {
        unsigned narg=mypack1_stashd_args.size(); Tensor vir; vir.zero(); 
        ActionAtomistic* mymap = dynamic_cast<ActionAtomistic*>( getPntrToArgument(0)->getPntrToAction() );
        for(unsigned i=0; i<mypack1_stashd_atoms.size(); ++i) {
          Vector dxder; for(unsigned k=0; k<3; ++k) dxder[k] = ( 2*v4v4*dx - 2*proj )*spacing*sp->getDerivative( narg + 3*i+k );
          Vector ader = ( mypack1_stashd_atoms[i] - 2.*dx*mypack1.getAtomDerivative(i) + dxder )/ (2.0*path_z);
          for(unsigned k=0; k<3; ++k) zp->setDerivative( narg+3*i+k, ader[k] );
          vir-=Tensor( mymap->getPosition(i), ader );
        }
        // Set the virial
        unsigned nbase=narg+3*mypack1_stashd_atoms.size();
        for(unsigned i=0; i<3; ++i) for(unsigned j=0; j<3; ++j) zp->setDerivative( nbase+3*i+j, vir(i,j) );
      }
  }
}

void GeometricPath::apply() {
  if( doNotCalculateDerivatives() ) return;
  std::fill(forcesToApply.begin(),forcesToApply.end(),0);
  if( getForcesFromValues( forcesToApply ) ){
      Mapping* am = dynamic_cast<Mapping*>( getPntrToArgument(0)->getPntrToAction() );
      am->setForcesOnArguments( forcesToApply );
      am->setForcesOnAtoms( forcesToApply, am->getNumberOfArguments()  );
  }
}

}
}


