/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2014 The plumed team
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
#include "SingleDomainRMSD.h"
#include "MultiDomainRMSD.h"
#include "MetricRegister.h"
#include "tools/PDB.h"

namespace PLMD {

PLUMED_REGISTER_METRIC(MultiDomainRMSD,"MULTI")

MultiDomainRMSD::MultiDomainRMSD( const ReferenceConfigurationOptions& ro ):
ReferenceConfiguration(ro),
ReferenceAtoms(ro)
{
   ftype=ro.getMultiRMSDType();
}

MultiDomainRMSD::~MultiDomainRMSD(){
  for(unsigned i=0;i<domains.size();++i) delete domains[i];
}

void MultiDomainRMSD::read( const PDB& pdb ){
   unsigned nblocks =  pdb.getNumberOfAtomBlocks();
   if( nblocks<2 ) error("multidomain RMSD only has one block of atoms");
  
   std::vector<AtomNumber> atomnumbers;
   std::vector<Vector> positions; std::vector<double> align, displace;
   std::string num; blocks.resize( nblocks+1 ); blocks[0]=0;
   for(unsigned i=0;i<nblocks;++i) blocks[i+1]=pdb.getAtomBlockEnds()[i]; 

   double lower, upper;
   if( !parse("LOWER_CUTOFF",lower,true) ) lower=0.0;
   if( !parse("UPPER_CUTTOFF",upper,true) ) upper=std::numeric_limits<double>::max( );

   for(unsigned i=1;i<=nblocks;++i){
       Tools::convert(i,num);
       if( ftype=="RMSD" ){
          parse("TYPE"+num, ftype );
          if( !parse("LOWER_CUTOFF"+num,lower,true) ) lower=0.0;
          if( !parse("UPPER_CUTTOFF"+num,upper,true) ) upper=std::numeric_limits<double>::max( );
       }
       domains.push_back( metricRegister().create<SingleDomainRMSD>( ftype ) );
       positions.resize( blocks[i] - blocks[i-1] );
       align.resize( blocks[i] - blocks[i-1] );
       displace.resize( blocks[i] - blocks[i-1] );
       unsigned n=0;
       for(unsigned j=blocks[i-1];j<blocks[i];++j){
           positions[n]=pdb.getPositions()[j];
           align[n]=pdb.getOccupancy()[j];
           displace[n]=pdb.getBeta()[j];
           n++;
       }
       domains[i-1]->setBoundsOnDistances( true, lower, upper );  // Currently no option for nopbc
       domains[i-1]->setReferenceAtoms( positions, align, displace );
       // domains[i-1]->setNumberOfAtoms( positions.size() );
       
       double ww=0; parse("WEIGHT"+num, ww, true );
       if( ww==0 ) weights.push_back( 1.0 );
       else weights.push_back( ww );
   }   
   // And set the atom numbers for this object
   setAtomNumbers( pdb.getAtomNumbers() );
}

void MultiDomainRMSD::setReferenceAtoms( const std::vector<Vector>& conf, const std::vector<double>& align_in, const std::vector<double>& displace_in ){
  plumed_error();
}

double MultiDomainRMSD::calculate( const std::vector<Vector>& pos, const Pbc& pbc, ReferenceValuePack& myder, const bool& squared ) const {
  //clearDerivatives(); 
  double totd=0.; Tensor tvirial; std::vector<Vector> mypos;
  MultiValue tvals( 2, myder.getNumberOfDerivatives() ); ReferenceValuePack tder( 0, 0, tvals );
  for(unsigned i=0;i<domains.size();++i){
     // Must extract appropriate positions here
     mypos.resize( blocks[i+1] - blocks[i] ); 
     unsigned n=0; tder.resize( 0, mypos.size() );
     for(unsigned j=blocks[i];j<blocks[i+1];++j){ tder.setAtomIndex(n,j); mypos[n]=pos[j]; n++; }

     // This actually does the calculation
     totd += weights[i]*domains[i]->calculate( mypos, pbc, tder, true );
     // Now merge the derivative
     myder.copyScaledDerivatives( 1, weights[i], tvals );
     // Make sure virial status is set correctly in output derivative pack
     // This is only done here so I do this by using class friendship
     if( tder.virialWasSet() ) myder.boxWasSet=true;
     // And clear
     tder.clear();
  }
  if( !squared ){
     totd=sqrt(totd); double xx=0.5/totd;
     myder.scaleAllDerivatives( xx );
  }
  if( !myder.updateComplete() ) myder.updateDynamicLists();
  return totd;
}

double MultiDomainRMSD::calc( const std::vector<Vector>& pos, const Pbc& pbc, const std::vector<Value*>& vals, const std::vector<double>& arg, 
                              ReferenceValuePack& myder, const bool& squared ) const {
  plumed_dbg_assert( vals.size()==0 && pos.size()==getNumberOfAtoms() && arg.size()==0 );
  return calculate( pos, pbc, myder, squared );
}

}
