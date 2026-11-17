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
#include "SecondaryStructureDRMSD.h"
#include "SecondaryStructureRMSD.h"
#include "core/Colvar.h"
#include "core/ActionRegister.h"
#include "tools/SwitchingFunction.h"
#include "tools/OpenMP.h"

#include <array>

namespace PLMD {
namespace secondarystructure {

// The documentation of these actions is in ANTIBETARMSD and PARABETARMSD, the shortcuts that create them.

/// The number of segments that are similar to a beta sheet, computed in a single action.
/// It gives the same result as the chain of actions created by ANTIBETARMSD and PARABETARMSD
/// (SECONDARY_STRUCTURE_(D)RMSD, the DISTANCE and CUSTOM of STRANDS_CUTOFF, LOWEST, LESS_THAN, SUM), but
/// - the segments are not stored: they are built from the lengths of the chains, as pairs of strands;
/// - the segments whose strands are further apart than STRANDS_CUTOFF cost one distance;
/// - the derivatives are accumulated on the atoms while the segments are computed, so they are
///   computed once (the chain computes the active segments again in the backward pass) and no vector
///   with an element per segment is created.
class SheetCount : public Colvar {
private:
  bool parallel;
  bool nopbc;
  bool usedrmsd;
  double strands_cutoff;
  SwitchingFunction switchingFunction;
  /// The first atom of the two strands of each segment, in the list of atoms
  std::vector<std::pair<unsigned,unsigned> > pairs;
  SecondaryStructureDRMSDInput<double> drmsd;
  SecondaryStructureRMSDInput rmsd;
  unsigned nstructures;
  /// Buffers of the derivatives of each OpenMP thread
  std::vector<std::vector<double> > omp_derivatives;
  std::vector<double> derivatives;
  template <class I>
  double segmentValue( const I& input, unsigned ipair, std::array<double,3*maxAtomsInSegment+9>& der ) const;
public:
  static void registerKeywords( Keywords& keys );
  explicit SheetCount(const ActionOptions&);
  void calculate() override;
};

PLUMED_REGISTER_ACTION(SheetCount,"ANTIBETARMSD_SCALAR")
PLUMED_REGISTER_ACTION(SheetCount,"PARABETARMSD_SCALAR")

void SheetCount::registerKeywords( Keywords& keys ) {
  Colvar::registerKeywords( keys );
  if( keys.getDisplayName()=="ANTIBETARMSD_SCALAR" ) {
    keys.setDisplayName("ANTIBETARMSD");
  } else {
    keys.setDisplayName("PARABETARMSD");
  }
  keys.add("atoms","ATOMS","the backbone atoms of the chains (N, CA, CB, C and O of each residue)");
  keys.add("compulsory","CHAIN_SIZES","the number of atoms of each chain in ATOMS");
  keys.add("compulsory","STYLE","all","the pairs of strands that are used: from the same chain (intra), from different chains (inter) or both (all)");
  keys.add("compulsory","TYPE","DRMSD","the distance from the reference structures: DRMSD, OPTIMAL, OPTIMAL-FAST or SIMPLE");
  keys.add("compulsory","SWITCH","the switching function of the distance from the reference structures that is summed over the segments");
  keys.add("optional","STRANDS_CUTOFF","the segments whose strands are further apart than this (between the CA atoms of their central residues) are not used");
  keys.add("numbered","STRUCTURE","the coordinates (in Angstrom) of the 30 atoms of a reference structure; with more than one, the distance of a segment is the lowest");
  keys.setValueDescription("scalar","the number of segments that are similar to a beta sheet");
}

SheetCount::SheetCount(const ActionOptions&ao):
  Action(ao),
  Colvar(ao),
  parallel(getName()=="PARABETARMSD_SCALAR"),
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
  std::string style;
  parse("STYLE",style);
  SecondaryStructureInput::getSheetPairs( this, chains, parallel, style, pairs );
  log.printf("  %zu pairs of strands that could form a %s sheet\n", pairs.size(), parallel ? "parallel" : "antiparallel" );

  parse("STRANDS_CUTOFF",strands_cutoff);
  if( strands_cutoff>0 ) {
    log.printf("  only the pairs of strands whose central CA atoms are closer than %f are computed\n", strands_cutoff );
  }
  std::string type;
  parse("TYPE",type);
  usedrmsd = (type=="DRMSD");
  parseFlag("NOPBC",nopbc);
  std::string sw, errors;
  parse("SWITCH",sw);
  switchingFunction.set(sw,errors);
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

template <class I>
double SheetCount::segmentValue( const I& input, const unsigned ipair, std::array<double,3*maxAtomsInSegment+9>& der ) const {
  const unsigned s1=pairs[ipair].first, s2=pairs[ipair].second;
  std::array<Vector,maxAtomsInSegment> pos;
  for(unsigned k=0; k<15; ++k) {
    pos[k]=getPosition(s1+k);
    pos[15+k]=getPosition(s2+k);
  }
  if( !nopbc ) {
    // As in SecondaryStructureBase::performTask with ALIGN_STRANDS: the strands are made whole separately,
    // then the second is moved to the image where its central CA atom is closest to that of the first
    for(unsigned start=0; start<maxAtomsInSegment; start+=15) {
      for(unsigned i=start; i+1<start+15; ++i) {
        pos[i+1]=pos[i]+pbcDistance(pos[i],pos[i+1]);
      }
    }
    const Vector shift = pbcDistance(pos[6],pos[21]) + pos[6] - pos[21];
    for(unsigned i=15; i<maxAtomsInSegment; ++i) {
      pos[i]+=shift;
    }
  }
  // the distance from each reference structure and its derivatives
  constexpr unsigned nder=3*maxAtomsInSegment+9;
  std::array<double,2> dist;
  std::array<double,2*nder> d;
  ColvarOutput<double> output( View<double>(dist.data(),nstructures), nder, d.data() );
  for(unsigned n=0; n<nstructures; ++n) {
    I::calculateDistance( n, false, input, View<Vector>(pos.data(),maxAtomsInSegment), output );
  }
  // the lowest distance, as LOWEST in PARABETARMSD
  unsigned best=0;
  for(unsigned n=1; n<nstructures; ++n) {
    if( dist[n]<dist[best] ) {
      best=n;
    }
  }
  std::copy( d.begin()+best*nder, d.begin()+(best+1)*nder, der.begin() );
  return dist[best];
}

void SheetCount::calculate() {
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
  if( omp_derivatives.size()<nt ) {
    omp_derivatives.resize(nt);
  }
  double value=0.0;
  #pragma omp parallel num_threads(nt) reduction(+:value)
  {
    const unsigned t=OpenMP::getThreadNum();
    std::vector<double>& myder( omp_derivatives[t] );
    myder.assign( 3*natoms+9, 0.0 );
    std::array<double,3*maxAtomsInSegment+9> der;
    #pragma omp for
    for(unsigned i=rank; i<npairs; i+=stride) {
      if( strands_cutoff>0 ) {
        const unsigned c1=pairs[i].first+6, c2=pairs[i].second+6;
        const Vector dca = nopbc ? delta( getPosition(c1), getPosition(c2) ) : pbcDistance( getPosition(c1), getPosition(c2) );
        // as step(STRANDS_CUTOFF-x) in the shortcut: the pair is used if the distance is not larger than the cutoff
        if( dca.modulo()>strands_cutoff ) {
          continue;
        }
      }
      const double r = usedrmsd ? segmentValue( drmsd, i, der ) : segmentValue( rmsd, i, der );
      double df=0.0;
      const double s=switchingFunction.calculate( r, df );
      value += s;
      // df is (ds/dr)/r
      const double g=df*r;
      if( g==0.0 ) {
        continue;
      }
      const unsigned s1=pairs[i].first, s2=pairs[i].second;
      for(unsigned k=0; k<maxAtomsInSegment; ++k) {
        const unsigned a = (k<15) ? s1+k : s2+k-15;
        for(unsigned j=0; j<3; ++j) {
          myder[3*a+j] += g*der[3*k+j];
        }
      }
      for(unsigned j=0; j<9; ++j) {
        myder[3*natoms+j] += g*der[3*maxAtomsInSegment+j];
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
  setValue( value );
  for(unsigned a=0; a<natoms; ++a) {
    setAtomsDerivatives( a, Vector( derivatives[3*a], derivatives[3*a+1], derivatives[3*a+2] ) );
  }
  setBoxDerivatives( Tensor( derivatives[3*natoms], derivatives[3*natoms+1], derivatives[3*natoms+2],
                             derivatives[3*natoms+3], derivatives[3*natoms+4], derivatives[3*natoms+5],
                             derivatives[3*natoms+6], derivatives[3*natoms+7], derivatives[3*natoms+8] ) );
}

}
}
