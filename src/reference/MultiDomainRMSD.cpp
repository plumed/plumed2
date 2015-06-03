/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2014,2015 The plumed team
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
ReferenceAtoms(ro),
ftype(ro.getMultiRMSDType())
{
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

   double lower=0.0, upper=std::numeric_limits<double>::max( );
   parse("LOWER_CUTOFF",lower,true); 
   parse("UPPER_CUTOFF",upper,true);

   for(unsigned i=1;i<=nblocks;++i){
       Tools::convert(i,num);
       if( ftype=="RMSD" ){
          parse("TYPE"+num, ftype );
          parse("LOWER_CUTOFF"+num,lower,true); 
          parse("UPPER_CUTOFF"+num,upper,true); 
       }
       domains.push_back( metricRegister().create<SingleDomainRMSD>( ftype ) );
       positions.resize( blocks[i] - blocks[i-1] + 1 );
       align.resize( blocks[i] - blocks[i-1] + 1 );
       displace.resize( blocks[i] - blocks[i-1] + 1 );
       unsigned n=0;
       for(unsigned j=blocks[i-1];j<blocks[i];++j){
           positions[n]=pdb.getPositions()[j];
           align[n]=pdb.getOccupancy()[j];
           displace[n]=pdb.getBeta()[j];
           n++;
       }
       domains[i-1]->setBoundsOnDistances( true, lower, upper );  // Currently no option for nopbc
       domains[i-1]->setReferenceAtoms( positions, align, displace );
       domains[i-1]->setNumberOfAtoms( positions.size() );
       
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

double MultiDomainRMSD::calculate( const std::vector<Vector>& pos, const Pbc& pbc,  const bool& squared ){
  clearDerivatives(); double totd=0.; Tensor tvirial; std::vector<Vector> mypos;
  for(unsigned i=0;i<domains.size();++i){
     // Must extract appropriate positions here
     mypos.resize( blocks[i+1] - blocks[i] + 1 );
     unsigned n=0;
     for(unsigned j=blocks[i];j<blocks[i+1];++j){ mypos[n]=pos[j]; n++; }

     // This actually does the calculation
     totd += weights[i]*domains[i]->calculate( mypos, pbc, true );

     // Must extract derivatives here
     n=0;
     for(unsigned j=blocks[i];j<blocks[i+1];++j){
         addAtomicDerivatives( j, weights[i]*domains[i]->getAtomDerivative(n) ); n++;
     }
     if( domains[i]->getVirial( tvirial ) ){
         addBoxDerivatives( weights[i]*tvirial );
     }
  }
  if( !squared ){
     totd=sqrt(totd); double xx=0.5/totd;
     for(unsigned iat=0;iat<atom_ders.size();iat++) atom_ders[iat]*=xx;
     if( virialWasSet ) virial*=xx;
  }
  return totd;
}

double MultiDomainRMSD::calc( const std::vector<Vector>& pos, const Pbc& pbc, const std::vector<Value*>& vals, const std::vector<double>& arg, const bool& squared ){
  plumed_dbg_assert( vals.size()==0 && pos.size()==getNumberOfAtoms() && arg.size()==0 );
  return calculate( pos, pbc, squared );
}

bool MultiDomainRMSD::pcaIsEnabledForThisReference(){
  bool enabled=true;
  for(unsigned i=0;i<domains.size();++i){
      if( !domains[i]->pcaIsEnabledForThisReference() ) enabled=false;
  }
  return enabled;
}

Vector MultiDomainRMSD::getAtomicDisplacement( const unsigned& iatom ){
  for(unsigned i=0;i<domains.size();++i){
      unsigned n=0;
      for(unsigned j=blocks[i];j<blocks[i+1];++j){
          if( j==iatom ) return weights[i]*domains[i]->getAtomicDisplacement(n);
          n++;
      }
  }
}

double MultiDomainRMSD::projectAtomicDisplacementOnVector( const unsigned& iv, const Matrix<Vector>& vecs, const std::vector<Vector>& pos, std::vector<Vector>& derivatives ){
  double totd=0.; std::vector<Vector> mypos, tderivs; Matrix<Vector> tvecs;
  for(unsigned i=0;i<domains.size();++i){
      // Must extract appropriate positions here 
     mypos.resize( blocks[i+1] - blocks[i] + 1 );
     tderivs.resize( blocks[i+1] - blocks[i] + 1 );
     tvecs.resize( vecs.ncols(), blocks[i+1] - blocks[i] + 1 );
     unsigned n=0;
     for(unsigned j=blocks[i];j<blocks[i+1];++j){
        mypos[n]=pos[j];
        for(unsigned k=0;k<vecs.ncols();++k) tvecs( k, n ) = vecs( k, j );
        n++; 
     }
      
     // Do the calculations
     domains[i]->projectAtomicDisplacementOnVector( iv, tvecs, mypos, tderivs );
      
     // And derivatives
     n=0;
     for(unsigned j=blocks[i];j<blocks[i+1];++j){
         addAtomicDerivatives( j, weights[i]*tderivs[n] ); n++;
     }
  }

  return totd;
}

}
