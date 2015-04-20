/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013,2014 The plumed team
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
#include "vesselbase/VesselRegister.h"
#include "vesselbase/ActionWithVessel.h"
#include "tools/DynamicList.h"
#include "MultiColvar.h"
#include "MultiColvarFunction.h"
#include "StoreCentralAtomsVessel.h"

namespace PLMD {
namespace multicolvar {

StoreCentralAtomsVessel::StoreCentralAtomsVessel( const vesselbase::VesselOptions& da ):
StoreDataVessel(da),
tmpdf(3)
{
  mycolv=dynamic_cast<MultiColvarBase*>( getAction() );
  plumed_assert( mycolv ); completeSetup( mycolv->getCentralAtomElementIndex(), 3 );
}

void StoreCentralAtomsVessel::getIndexList( const unsigned& ntotal, const unsigned& jstore, const unsigned& maxder, std::vector<unsigned>& aindexes ){
  plumed_dbg_assert( mycolv->derivativesAreRequired() );

  aindexes[jstore]=3*mycolv->atomsWithCatomDer.getNumberActive();
  if( aindexes[jstore]>maxder ) error("too many derivatives to store. Run with LOWMEM");
  unsigned kder = ntotal + jstore*maxder;
  for(unsigned jder=0;jder<mycolv->atomsWithCatomDer.getNumberActive();++jder){
     unsigned iatom = 3*mycolv->atomsWithCatomDer[jder];
     for(unsigned icomp=0;icomp<3;++icomp){ aindexes[ kder ] = iatom+icomp; kder++; }
  }
}

Vector StoreCentralAtomsVessel::getPosition( const unsigned& iatom ){
  Vector pos; for(unsigned i=0;i<3;++i) pos[i]=getComponent( iatom, i );
  return pos;
}

void StoreCentralAtomsVessel::performTask( const unsigned& itask ){
  mycolv->atomsWithCatomDer.deactivateAll();
  bool check=mycolv->setupCurrentAtomList( mycolv->getCurrentTask() );
  plumed_dbg_assert( check );
  Vector ignore = mycolv->retrieveCentralAtomPos();
}

void StoreCentralAtomsVessel::finishTask( const unsigned& itask ){
  mycolv->atomsWithCatomDer.deactivateAll();
  Vector ignore = mycolv->retrieveCentralAtomPos();
}

void StoreCentralAtomsVessel::addAtomsDerivatives( const unsigned& iatom, const unsigned& jout, const unsigned& base_cv_no, 
                                                   const Vector& df, MultiColvarFunction* funcout ){
  plumed_dbg_assert( mycolv->derivativesAreRequired() );

  for(unsigned ider=0;ider<getNumberOfDerivatives(iatom);ider+=3){
     for(unsigned i=0;i<3;++i) tmpdf[i]=df[0];
     funcout->addStoredDerivative( jout, base_cv_no, getStoredIndex( iatom, ider+0 ), chainRule(iatom, ider+0, tmpdf)  ); 
     for(unsigned i=0;i<3;++i) tmpdf[i]=df[1];
     funcout->addStoredDerivative( jout, base_cv_no, getStoredIndex( iatom, ider+1 ), chainRule(iatom, ider+1, tmpdf)  );
     for(unsigned i=0;i<3;++i) tmpdf[i]=df[2];
     funcout->addStoredDerivative( jout, base_cv_no, getStoredIndex( iatom, ider+2 ), chainRule(iatom, ider+2, tmpdf)  );
  }
}

}
}
