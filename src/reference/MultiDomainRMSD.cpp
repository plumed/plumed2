/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2020 The plumed team
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
#include "MultiDomainRMSD.h"
#include "SingleDomainRMSD.h"
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

void MultiDomainRMSD::read( const PDB& pdb ) {
  unsigned nblocks =  pdb.getNumberOfAtomBlocks();
  if( nblocks<2 ) error("multidomain RMSD only has one block of atoms");

  std::vector<Vector> positions; std::vector<double> align, displace;
  std::string num; blocks.resize( nblocks+1 ); blocks[0]=0;
  for(unsigned i=0; i<nblocks; ++i) blocks[i+1]=pdb.getAtomBlockEnds()[i];

  double tmp, lower=0.0, upper=std::numeric_limits<double>::max( );
  if( pdb.getArgumentValue("LOWER_CUTOFF",tmp) ) lower=tmp;
  if( pdb.getArgumentValue("UPPER_CUTOFF",tmp) ) upper=tmp;
  bool nopbc=pdb.hasFlag("NOPBC");

  domains.resize(0); weights.resize(0);
  for(unsigned i=1; i<=nblocks; ++i) {
    Tools::convert(i,num);
    if( ftype=="RMSD" ) {
      // parse("TYPE"+num, ftype );
      lower=0.0; upper=std::numeric_limits<double>::max( );
      if( pdb.getArgumentValue("LOWER_CUTOFF"+num,tmp) ) lower=tmp;
      if( pdb.getArgumentValue("UPPER_CUTOFF"+num,tmp) ) upper=tmp;
      nopbc=pdb.hasFlag("NOPBC");
    }
    domains.emplace_back( metricRegister().create<SingleDomainRMSD>( ftype ) );
    positions.resize( blocks[i] - blocks[i-1] );
    align.resize( blocks[i] - blocks[i-1] );
    displace.resize( blocks[i] - blocks[i-1] );
    unsigned n=0;
    for(unsigned j=blocks[i-1]; j<blocks[i]; ++j) {
      positions[n]=pdb.getPositions()[j];
      align[n]=pdb.getOccupancy()[j];
      displace[n]=pdb.getBeta()[j];
      n++;
    }
    domains[i-1]->setBoundsOnDistances( !nopbc, lower, upper );
    domains[i-1]->setReferenceAtoms( positions, align, displace );
    domains[i-1]->setupRMSDObject();

    double ww=0;
    if( !pdb.getArgumentValue("WEIGHT"+num,ww) ) weights.push_back( 1.0 );
    else weights.push_back( ww );
  }
  // And set the atom numbers for this object
  indices.resize(0); atom_der_index.resize(0);
  for(unsigned i=0; i<pdb.size(); ++i) { indices.push_back( pdb.getAtomNumbers()[i] ); atom_der_index.push_back(i); }
  // setAtomNumbers( pdb.getAtomNumbers() );
}

void MultiDomainRMSD::setReferenceAtoms( const std::vector<Vector>& conf, const std::vector<double>& align_in, const std::vector<double>& displace_in ) {
  plumed_error();
}

double MultiDomainRMSD::calculate( const std::vector<Vector>& pos, const Pbc& pbc, ReferenceValuePack& myder, const bool& squared ) const {
  double totd=0.; Tensor tvirial; std::vector<Vector> mypos; MultiValue tvals( 1, 3*pos.size()+9 );
  ReferenceValuePack tder( 0, getNumberOfAtoms(), tvals ); myder.clear();

  for(unsigned i=0; i<domains.size(); ++i) {
    // Must extract appropriate positions here
    mypos.resize( blocks[i+1] - blocks[i] );
    if( myder.calcUsingPCAOption() ) domains[i]->setupPCAStorage( tder );
    unsigned n=0; for(unsigned j=blocks[i]; j<blocks[i+1]; ++j) { tder.setAtomIndex(n,j); mypos[n]=pos[j]; n++; }
    for(unsigned k=n; k<getNumberOfAtoms(); ++k) tder.setAtomIndex(k,3*pos.size()+10);
    // This actually does the calculation
    totd += weights[i]*domains[i]->calculate( mypos, pbc, tder, true );
    // Now merge the derivative
    myder.copyScaledDerivatives( 0, weights[i], tvals );
    // If PCA copy PCA stuff
    if( myder.calcUsingPCAOption() ) {
      unsigned n=0;
      if( tder.centeredpos.size()>0 ) myder.rot[i]=tder.rot[0];
      for(unsigned j=blocks[i]; j<blocks[i+1]; ++j) {
        myder.displacement[j]=weights[i]*tder.displacement[n];  // Multiplication by weights here ensures that normalisation is done correctly
        if( tder.centeredpos.size()>0 ) {
          myder.centeredpos[j]=tder.centeredpos[n];
          for(unsigned p=0; p<3; ++p) for(unsigned q=0; q<3; ++q) myder.DRotDPos(p,q)[j]=tder.DRotDPos(p,q)[n];
        }
        n++;
      }
    }
    // Make sure virial status is set correctly in output derivative pack
    // This is only done here so I do this by using class friendship
    if( tder.virialWasSet() ) myder.boxWasSet=true;
  }
  if( !myder.updateComplete() ) myder.updateDynamicLists();

  if( !squared ) {
    totd=sqrt(totd); double xx=0.5/totd;
    myder.scaleAllDerivatives( xx );
  }
  return totd;
}

double MultiDomainRMSD::calc( const std::vector<Vector>& pos, const Pbc& pbc, const std::vector<Value*>& vals, const std::vector<double>& arg,
                              ReferenceValuePack& myder, const bool& squared ) const {
  plumed_dbg_assert( vals.size()==0 && pos.size()==getNumberOfAtoms() && arg.size()==0 );
  return calculate( pos, pbc, myder, squared );
}

bool MultiDomainRMSD::pcaIsEnabledForThisReference() {
  bool enabled=true;
  for(unsigned i=0; i<domains.size(); ++i) {
    if( !domains[i]->pcaIsEnabledForThisReference() ) enabled=false;
  }
  return enabled;
}

void MultiDomainRMSD::setupPCAStorage( ReferenceValuePack& mypack ) {
  plumed_dbg_assert( pcaIsEnabledForThisReference() );
  mypack.switchOnPCAOption();
  mypack.displacement.resize( getNumberOfAtoms() );
  mypack.centeredpos.resize( getNumberOfAtoms() );
  mypack.DRotDPos.resize(3,3); mypack.rot.resize( domains.size() );
  for(unsigned i=0; i<3; ++i) for(unsigned j=0; j<3; ++j) mypack.DRotDPos(i,j).resize( getNumberOfAtoms() );
}

// Vector MultiDomainRMSD::getAtomicDisplacement( const unsigned& iatom ){
//   for(unsigned i=0;i<domains.size();++i){
//       unsigned n=0;
//       for(unsigned j=blocks[i];j<blocks[i+1];++j){
//           if( j==iatom ) return weights[i]*domains[i]->getAtomicDisplacement(n);
//           n++;
//       }
//   }
// }

void MultiDomainRMSD::extractAtomicDisplacement( const std::vector<Vector>& pos, std::vector<Vector>& direction ) const {
  std::vector<Vector> mypos, mydir;
  for(unsigned i=0; i<domains.size(); ++i) {
    // Must extract appropriate positions here
    mypos.resize( blocks[i+1] - blocks[i] ); mydir.resize( blocks[i+1] - blocks[i] );
    unsigned n=0; for(unsigned j=blocks[i]; j<blocks[i+1]; ++j) { mypos[n]=pos[j]; n++; }
    // Do the calculation
    domains[i]->extractAtomicDisplacement( mypos, mydir );
    // Extract the direction
    n=0; for(unsigned j=blocks[i]; j<blocks[i+1]; ++j) { direction[j]=weights[i]*mydir[n];  n++; }
  }
}

double MultiDomainRMSD::projectAtomicDisplacementOnVector( const bool& normalized, const std::vector<Vector>& vecs, ReferenceValuePack& mypack ) const {
  double totd=0.; std::vector<Vector> tvecs; mypack.clear();
  MultiValue tvals( 1, mypack.getNumberOfDerivatives() ); ReferenceValuePack tder( 0, getNumberOfAtoms(), tvals );
  for(unsigned i=0; i<domains.size(); ++i) {
    // Must extract appropriate positions here
    tvecs.resize( blocks[i+1] - blocks[i] ); domains[i]->setupPCAStorage( tder );
    if( tder.centeredpos.size()>0 ) {
      for(unsigned p=0; p<3; ++p) for(unsigned q=0; q<3; ++q) tder.DRotDPos(p,q).resize( tvecs.size() );
    }
    // Extract information from storage pack and put in local pack
    if( tder.centeredpos.size()>0 ) tder.rot[0]=mypack.rot[i];
    unsigned n=0;
    for(unsigned j=blocks[i]; j<blocks[i+1]; ++j) {
      tder.setAtomIndex(n,j); tvecs[n] = vecs[j]; tder.displacement[n]=mypack.displacement[j] / weights[i];
      if( tder.centeredpos.size()>0 ) {
        tder.centeredpos[n]=mypack.centeredpos[j];
        for(unsigned p=0; p<3; ++p) for(unsigned q=0; q<3; ++q) tder.DRotDPos(p,q)[n]=mypack.DRotDPos(p,q)[j];
      }
      n++;
    }
    for(unsigned k=n; k<getNumberOfAtoms(); ++k) tder.setAtomIndex(k,3*vecs.size()+10);

    // Do the calculations
    totd += weights[i]*domains[i]->projectAtomicDisplacementOnVector( normalized, tvecs, tder );

    // And derivatives
    mypack.copyScaledDerivatives( 0, weights[i], tvals );
  }
  if( !mypack.updateComplete() ) mypack.updateDynamicLists();

  return totd;
}

}
