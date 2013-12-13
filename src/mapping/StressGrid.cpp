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
#include "vesselbase/VesselRegister.h"
#include "vesselbase/FieldGridBase.h"
#include "Mapping.h"

namespace PLMD {
namespace mapping {

class StressGrid : public vesselbase::FieldGridBase {
private:
  bool firsttime;
  std::vector<double> stash;
  Mapping* map;
public:
  static void registerKeywords( Keywords& keys );
  static void reserveKeyword( Keywords& keys );
  StressGrid( const vesselbase::VesselOptions& );
  std::string description();
  void resize();
  void prepare();
  bool calculate();
  void finish();
};

PLUMED_REGISTER_VESSEL(StressGrid,"STRESS_GRID")

void StressGrid::registerKeywords( Keywords& keys ){
  FieldGridBase::registerKeywords( keys );
}

void StressGrid::reserveKeyword( Keywords& keys ){
  keys.reserve("compulsory","STRESS_GRID","the input for a field cv");  
}

StressGrid::StressGrid( const vesselbase::VesselOptions& da ):
FieldGridBase(da),
firsttime(true)
{
  map=dynamic_cast<Mapping*>( getAction() );
  plumed_assert( map );

  unsigned nprop = map->getNumberOfProperties(); 
  unsigned nargs = map->getNumberOfArguments(); 
  if( map->getNumberOfAtoms()>0 ) nargs += 9 + 3*map->getNumberOfAtoms();

  unsigned ntot = nprop + 1 + nprop + nargs + nargs*nprop;
  std::vector<std::string> names( ntot ); unsigned k=0;
  for(unsigned i=0;i<nprop;++i){ names[k]=map->getPropertyName( i ); k++; } 
  names[k]="chi2"; k++;
  for(unsigned i=0;i<nprop;++i){ names[k]="dchi2_" + map->getPropertyName( i ); k++; }

  for(unsigned j=0;j<nargs;++j){
      names[k]="dchi2_" + map->getArgumentName(j); k++;
      for(unsigned i=0;i<nprop;++i){
          names[k]="d2chi2_" + map->getArgumentName(j) + "_" + map->getPropertyName( i ); k++; 
      }
  } 
  std::vector<bool> mypbc( dimension, false );
  finishSetup( names.size()-nprop , mypbc, names );
}

std::string StressGrid::description(){
  return "calculating stress on " + getGridDescription();
}

void StressGrid::prepare(){
  clearForces();

  if(firsttime) return ; 
  unsigned stride=comm.Get_size(); unsigned rank=comm.Get_rank();
  unsigned n=0; for(unsigned i=rank;i<getSizeOfBuffer();i+=stride){ addToBufferElement( i, stash[n]); n++; }
}

void StressGrid::resize(){
  vesselbase::FieldGridBase::resize();
  stash.resize( getSizeOfBuffer() );
}

bool StressGrid::calculate(){
  // Retrieve the high dimensional function
  double highdv = map->getCurrentHighDimFunctionValue( 0 );
  double lowdv = map->getElementValue(0);
  double weight = map->getElementValue(1);

  // And compute the stress
  double tmp = highdv - lowdv;
  double chi2 = weight*tmp*tmp;

  // And add everything to the field stress
  accumulate( chi2, -2.*weight*tmp, 2.*weight*tmp, -2.*weight, 0 ); 
  if(firsttime) return true;

  // Retrieve the old distance from this point
  highdv = map->getCurrentHighDimFunctionValue( 1 );
  // And compute the stress
  tmp = highdv - lowdv;
  double ochi2 = weight*tmp*tmp;

  // And remove everything from the field stress
  accumulate( -ochi2, 2.*weight*tmp, -2.*weight*tmp, 2.*weight, 1 ); 
  return true;
}

void StressGrid::finish(){
  // Store what is currently in the buffers
  unsigned stride=comm.Get_size(); unsigned rank=comm.Get_rank();
  unsigned n=0; for(unsigned i=rank;i<getSizeOfBuffer();i+=stride){ stash[n]=getBufferElement(i); n++; }
  firsttime=false;
}

}
}
