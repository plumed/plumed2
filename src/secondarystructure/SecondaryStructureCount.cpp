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

namespace PLMD {
namespace secondarystructure {

// The documentation of these actions is in ALPHARMSD, ANTIBETARMSD and PARABETARMSD, the shortcuts that create them.

// SecondaryStructureCount gives the same result as the chain of actions created by ALPHARMSD, ANTIBETARMSD and PARABETARMSD
// (SECONDARY_STRUCTURE_(D)RMSD, the DISTANCE and CUSTOM of STRANDS_CUTOFF, LOWEST, LESS_THAN, SUM), but
// - the segments are not stored: they are built from the lengths of the chains (for the sheets, as pairs of strands);
// - the segments whose strands are further apart than STRANDS_CUTOFF cost one distance;
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

  // The reference structures, as in SECONDARY_STRUCTURE_RMSD and SECONDARY_STRUCTURE_DRMSD
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
      // the bond length of the shortcuts, 0.17 nm
      drmsd.setReferenceStructure( type, 0.17/getUnits().getLength(), structure );
    } else {
      rmsd.setReferenceStructure( type, 0.0, structure );
    }
    nstructures++;
  }
  if( nstructures==0 || nstructures>2 ) {
    error("there must be one or two STRUCTURE keywords");
  }
  checkRead();
  log<<"  Bibliography "<<plumed.cite("Pietrucci and Laio, J. Chem. Theory Comput. 5, 2197 (2009)")<<"\n";

  addValueWithDerivatives();
  setNotPeriodic();
  requestAtoms(atoms);
  derivatives.resize( 3*atoms.size()+9 );
}

void SecondaryStructureCount::calculate() {
  const unsigned npairs=pairs.size();
  const unsigned natoms=getNumberOfAtoms();
  unsigned stride=comm.Get_size();
  unsigned rank=comm.Get_rank();
  unsigned nt=OpenMP::getNumThreads();
  if( nt*stride*10>npairs ) {
    nt=npairs/stride/10;
  }
  if( nt==0 ) {
    nt=1;
  }
  // The buffers are cleared here and not by the threads, because all the nt buffers are summed below and fewer
  // threads may run: OpenMP may give fewer than nt, and without OpenMP a single one runs (nt still comes from
  // PLUMED_NUM_THREADS)
  if( omp_derivatives.size()<nt ) {
    omp_derivatives.resize(nt);
  }
  for(unsigned t=0; t<nt; ++t) {
    omp_derivatives[t].assign( 3*natoms+9, 0.0 );
  }
  double value=0.0;
  #pragma omp parallel num_threads(nt) reduction(+:value)
  {
    const unsigned t=OpenMP::getThreadNum();
    std::vector<double>& myder( omp_derivatives[t] );
    std::array<Vector,maxAtomsInSegment> pos;
    std::array<double,2*sscount::segmentDerivatives> der;
    const double* x=&getPositions()[0][0];
    const Pbc* pbc = nopbc ? nullptr : &getPbc();
    #pragma omp for
    for(unsigned i=rank; i<npairs; i+=stride) {
      const unsigned s1=pairs[i].first, s2=pairs[i].second;
      if( strands_cutoff>0 && !sscount::withinCutoff( x, s1, s2, pbc, strands_cutoff ) ) {
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
  }
  std::fill( derivatives.begin(), derivatives.end(), 0.0 );
  for(unsigned t=0; t<nt; ++t) {
    for(unsigned j=0; j<derivatives.size(); ++j) {
      derivatives[j] += omp_derivatives[t][j];
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
