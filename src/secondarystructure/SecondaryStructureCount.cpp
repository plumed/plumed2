/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2026 The plumed team
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
#include "SecondaryStructureCount.h"
#include "core/ActionRegister.h"
#include "tools/OpenMP.h"

#include <algorithm>
#include <limits>

namespace PLMD {
namespace secondarystructure {

// The documentation of these actions is in ALPHARMSD, ANTIBETARMSD and PARABETARMSD, the shortcuts that create them.

// SecondaryStructureCount gives the same result as the chain of actions created by ALPHARMSD, ANTIBETARMSD and PARABETARMSD
// (SECONDARY_STRUCTURE_(D)RMSD, the DISTANCE and CUSTOM of STRANDS_CUTOFF, LOWEST, LESS_THAN, SUM), but
// - the segments are not stored: they are built from the lengths of the chains (for the sheets, as pairs of strands);
// - the segments whose strands are further apart than STRANDS_CUTOFF cost one distance, and so do the segments
//   whose switching function is zero because their strands are far apart (see setSkipDistance);
// - the derivatives are accumulated on the atoms while the segments are computed, so they are
//   computed once (the chain computes the active segments again in the backward pass), only for the
//   segments where the switching function is not constant, and no vector with an element per segment is created.

PLUMED_REGISTER_ACTION(SecondaryStructureCount,"ALPHARMSD_SCALAR")
PLUMED_REGISTER_ACTION(SecondaryStructureCount,"ANTIBETARMSD_SCALAR")
PLUMED_REGISTER_ACTION(SecondaryStructureCount,"PARABETARMSD_SCALAR")

void SecondaryStructureCount::registerKeywords( Keywords& keys ) {
  Colvar::registerKeywords( keys );
  // the names are ALPHARMSD_SCALAR, ANTIBETARMSD_SCALAR and PARABETARMSD_SCALAR, with _ACC for the GPU version
  const bool helix=keys.getDisplayName().rfind("ALPHARMSD",0)==0;
  if( helix ) {
    keys.setDisplayName("ALPHARMSD");
  } else if( keys.getDisplayName().rfind("ANTIBETARMSD",0)==0 ) {
    keys.setDisplayName("ANTIBETARMSD");
  } else {
    keys.setDisplayName("PARABETARMSD");
  }
  keys.add("atoms","ATOMS","the backbone atoms of the chains (N, CA, CB, C and O of each residue)");
  keys.add("compulsory","CHAIN_SIZES","the number of atoms of each chain in ATOMS");
  if( !helix ) {
    keys.add("compulsory","STYLE","all","the pairs of strands that are used: from the same chain (intra), from different chains (inter) or both (all)");
    keys.add("optional","STRANDS_CUTOFF","the segments whose strands are further apart than this (between the CA atoms of their central residues) are not used");
  }
  keys.add("compulsory","TYPE","DRMSD","the distance from the reference structures: DRMSD, OPTIMAL, OPTIMAL-FAST or SIMPLE");
  keys.add("compulsory","SWITCH","the switching function of the distance from the reference structures that is summed over the segments");
  keys.add("numbered","STRUCTURE","the coordinates (in Angstrom) of the 30 atoms of a reference structure; with more than one, the distance of a segment is the lowest");
  keys.setValueDescription("scalar",helix ? "the number of segments that are similar to an alpha helix" : "the number of segments that are similar to a beta sheet");
}

SecondaryStructureCount::SecondaryStructureCount(const ActionOptions&ao):
  Action(ao),
  Colvar(ao),
  alpha(getName().rfind("ALPHARMSD",0)==0),
  parallel(getName().rfind("PARABETARMSD",0)==0),
  nopbc(false),
  usedrmsd(true),
  strands_cutoff(-1.0),
  skipDistance(-1.0),
  nstructures(0) {
  if( plumed.usingNaturalUnits() ) {
    error("cannot use this collective variable when using natural units");
  }
  std::vector<AtomNumber> atoms;
  parseAtomList("ATOMS",atoms);
  std::vector<unsigned> chains;
  parseVector("CHAIN_SIZES",chains);
  unsigned ntot=0;
  for(const auto c : chains) {
    ntot+=c;
  }
  if( ntot!=atoms.size() ) {
    error("the sum of CHAIN_SIZES is not the number of atoms in ATOMS");
  }
  if( alpha ) {
    // all the segments of six consecutive residues of each chain, as in ALPHARMSD
    unsigned nprevious=0;
    for(const auto c : chains) {
      if( c<maxAtomsInSegment || c%5!=0 ) {
        error("each chain must contain at least 6 residues of 5 backbone atoms");
      }
      for(unsigned ires=0; ires+5<c/5; ++ires) {
        pairs.emplace_back( nprevious+5*ires, nprevious+5*ires+sscount::strandAtoms );
      }
      nprevious+=c;
    }
    log.printf("  %zu segments of six residues that could form an alpha helix\n", pairs.size() );
  } else {
    std::string style;
    parse("STYLE",style);
    SecondaryStructureInput::getSheetPairs( this, chains, parallel, style, pairs );
    log.printf("  %zu pairs of strands that could form a %s sheet\n", pairs.size(), parallel ? "parallel" : "antiparallel" );

    parse("STRANDS_CUTOFF",strands_cutoff);
    if( strands_cutoff>0 ) {
      log.printf("  only the pairs of strands whose central CA atoms are closer than %f are computed\n", strands_cutoff );
    }
  }
  std::string type;
  parse("TYPE",type);
  usedrmsd = (type=="DRMSD");
  parseFlag("NOPBC",nopbc);
  std::string errors;
  parse("SWITCH",switchInput);
  switchingFunction.set(switchInput,errors);
  if( errors.length()>0 ) {
    error("problem reading SWITCH keyword : " + errors );
  }
  log.printf("  sum of %s of the %s distance from the reference structure\n", switchingFunction.description().c_str(), type.c_str() );
  if( switchingFunction.get_dmax()<std::numeric_limits<double>::max() ) {
    log.printf("  the switching function is zero for distances from the reference structure larger than D_MAX=%f\n", switchingFunction.get_dmax() );
  }

  // The reference structures, as in SECONDARY_STRUCTURE_RMSD and SECONDARY_STRUCTURE_DRMSD
  // the bond length of the shortcuts, 0.17 nm
  const double bondlength=0.17/getUnits().getLength();
  std::vector<std::vector<Vector> > references;
  for(unsigned ii=1;; ++ii) {
    std::vector<double> cstruct;
    if( !parseNumberedVector("STRUCTURE",ii,cstruct) ) {
      break ;
    }
    if( cstruct.size()!=3*maxAtomsInSegment ) {
      error("STRUCTURE" + std::to_string(ii) + " must contain the three coordinates of " + std::to_string(maxAtomsInSegment) + " atoms");
    }
    std::vector<Vector> structure( maxAtomsInSegment );
    for(unsigned i=0; i<structure.size(); ++i) {
      for(unsigned j=0; j<3; ++j) {
        structure[i][j] = 0.1*cstruct[3*i+j]/getUnits().getLength();
      }
    }
    if( usedrmsd ) {
      drmsd.setReferenceStructure( type, bondlength, structure );
    } else {
      rmsd.setReferenceStructure( type, 0.0, structure );
    }
    references.push_back( structure );
    nstructures++;
  }
  if( nstructures==0 || nstructures>2 ) {
    error("there must be one or two STRUCTURE keywords");
  }
  setSkipDistance( references, bondlength );
  checkRead();
  log<<"  Bibliography "<<plumed.cite("Pietrucci and Laio, J. Chem. Theory Comput. 5, 2197 (2009)")<<"\n";

  addValueWithDerivatives();
  setNotPeriodic();
  requestAtoms(atoms);
  derivatives.resize( 3*atoms.size()+9 );
}

void SecondaryStructureCount::setSkipDistance( const std::vector<std::vector<Vector> >& references, const double bondlength ) {
  skipDistance=-1.0;
  const double dmax=switchingFunction.get_dmax();
  if( alpha || !(dmax<std::numeric_limits<double>::max()) ) {
    return;
  }
  // The switching function is zero for a distance r>dmax. A segment whose strands have their central CA atoms at a
  // distance d, and radii r1 and r2 (sscount::strandRadius), has all the distances between an atom of the first strand
  // and one of the second larger than or equal to L=d-r1-r2. A lower bound of its distance from a reference is then:
  // - DRMSD: of the npairs distances of the reference, ncross are between the two strands, all at most d0max. If
  //   L>d0max, each of them gives a term at least (L-d0max)^2, so DRMSD >= sqrt(ncross/npairs)*(L-d0max);
  // - RMSD: after the optimal alignment the RMSD is at least the difference of the radii of gyration of the segment
  //   and of the reference, and the radius of gyration of the segment is at least half the distance between the
  //   centers of the strands, which is at least L: RMSD >= L/2-rg.
  // The segment can be skipped if the bound is larger than dmax for all the references, i.e. if L is larger than
  // d0max+dmax*sqrt(npairs/ncross) for DRMSD and 2*(dmax+rg) for RMSD.
  double skip=0.0;
  for(const auto& ref : references) {
    double s;
    if( usedrmsd ) {
      unsigned npairs=0, ncross=0;
      double d0max=0.0;
      for(unsigned i=0; i<ref.size(); ++i) {
        for(unsigned j=i+1; j<ref.size(); ++j) {
          const double d0=delta( ref[i], ref[j] ).modulo();
          if( d0>bondlength ) {
            npairs++;
            if( i<sscount::strandAtoms && j>=sscount::strandAtoms ) {
              ncross++;
              d0max=std::max( d0max, d0 );
            }
          }
        }
      }
      if( ncross==0 ) {
        return;
      }
      s = d0max + dmax*std::sqrt( static_cast<double>(npairs)/ncross );
    } else {
      Vector center;
      center.zero();
      for(const auto& r : ref) {
        center += r;
      }
      center /= ref.size();
      double rg2=0.0;
      for(const auto& r : ref) {
        rg2 += delta( center, r ).modulo2();
      }
      s = 2.0*( dmax + std::sqrt( rg2/ref.size() ) );
    }
    skip=std::max( skip, s );
  }
  // a margin for the rounding
  skip *= 1.0+1e-6;
  if( strands_cutoff>0 && strands_cutoff<=skip ) {
    // STRANDS_CUTOFF already skips these pairs
    return;
  }
  skipDistance=skip;
  for(const auto& p : pairs) {
    strands.push_back( p.first );
    strands.push_back( p.second );
  }
  std::sort( strands.begin(), strands.end() );
  strands.erase( std::unique( strands.begin(), strands.end() ), strands.end() );
  log.printf("  the pairs of strands whose central CA atoms are further apart than %f plus the radii of the strands are skipped, because their switching function is zero\n", skipDistance );
}

void SecondaryStructureCount::calculate() {
  const unsigned npairs=pairs.size();
  const unsigned natoms=getNumberOfAtoms();
  unsigned stride=comm.Get_size();
  unsigned rank=comm.Get_rank();
  // a pair of strands costs much more than reading it: at least 10 pairs per thread (see OpenMP.h)
  const unsigned nt=OpenMP::getGoodNumThreads(npairs/stride,10);
  if( omp_derivatives.size()<nt ) {
    omp_derivatives.resize(nt);
  }
  const unsigned nder=3*natoms+9;
  const double* x=&getPositions()[0][0];
  const Pbc* pbc = nopbc ? nullptr : &getPbc();
  const double* rad=nullptr;
  if( skipDistance>=0 ) {
    radius.resize( natoms );
    rad=radius.data();
  }
  const bool skip = strands_cutoff>0 || rad;
  double value=0.0;
  // The buffers of the threads are cleared and summed with loops over the nt buffers and over the derivatives
  // shared between the threads: all the nt buffers are summed, and fewer threads may run (OpenMP may give fewer
  // than nt, and without OpenMP a single one runs, while nt still comes from PLUMED_NUM_THREADS)
  #pragma omp parallel num_threads(nt) reduction(+:value)
  {
    #pragma omp for
    for(unsigned t=0; t<nt; ++t) {
      omp_derivatives[t].assign( nder, 0.0 );
    }
    if( rad ) {
      #pragma omp for
      for(unsigned i=0; i<strands.size(); ++i) {
        radius[strands[i]]=sscount::strandRadius( x, strands[i], pbc );
      }
    }
    const unsigned t=OpenMP::getThreadNum();
    std::vector<double>& myder( omp_derivatives[t] );
    std::array<Vector,maxAtomsInSegment> pos;
    std::array<double,2*sscount::segmentDerivatives> der;
    #pragma omp for
    for(unsigned i=rank; i<npairs; i+=stride) {
      const unsigned s1=pairs[i].first, s2=pairs[i].second;
      if( skip && !sscount::segmentActive( x, s1, s2, pbc, strands_cutoff, rad, skipDistance ) ) {
        continue;
      }
      sscount::segmentPositions( x, s1, s2, pbc, !alpha, pos );
      unsigned best=0;
      double g=0.0;
      value += usedrmsd ? sscount::segmentSwitch( drmsd, nstructures, switchingFunction, pos, der, best, g )
               : sscount::segmentSwitch( rmsd, nstructures, switchingFunction, pos, der, best, g );
      if( g==0.0 ) {
        continue;
      }
      const double* d=der.data()+best*sscount::segmentDerivatives;
      for(unsigned k=0; k<maxAtomsInSegment; ++k) {
        const unsigned a = (k<sscount::strandAtoms) ? s1+k : s2+k-sscount::strandAtoms;
        for(unsigned j=0; j<3; ++j) {
          myder[3*a+j] += g*d[3*k+j];
        }
      }
      for(unsigned j=0; j<9; ++j) {
        myder[3*natoms+j] += g*d[3*maxAtomsInSegment+j];
      }
    }
    // the sum of the buffers, in the same order for each derivative
    #pragma omp for
    for(unsigned j=0; j<nder; ++j) {
      double sum=0.0;
      for(unsigned k=0; k<nt; ++k) {
        sum += omp_derivatives[k][j];
      }
      derivatives[j]=sum;
    }
  }
  if( stride>1 ) {
    comm.Sum( value );
    comm.Sum( derivatives );
  }
  setValueAndDerivatives( value, derivatives );
}

void SecondaryStructureCount::setValueAndDerivatives( const double value, const std::vector<double>& der ) {
  const unsigned natoms=getNumberOfAtoms();
  setValue( value );
  for(unsigned a=0; a<natoms; ++a) {
    setAtomsDerivatives( a, Vector( der[3*a], der[3*a+1], der[3*a+2] ) );
  }
  setBoxDerivatives( Tensor( der[3*natoms], der[3*natoms+1], der[3*natoms+2],
                             der[3*natoms+3], der[3*natoms+4], der[3*natoms+5],
                             der[3*natoms+6], der[3*natoms+7], der[3*natoms+8] ) );
}

}
}
