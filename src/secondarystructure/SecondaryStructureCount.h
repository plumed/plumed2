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
#ifndef __PLUMED_secondarystructure_SecondaryStructureCount_h
#define __PLUMED_secondarystructure_SecondaryStructureCount_h
#include "SecondaryStructureDRMSD.h"
#include "SecondaryStructureRMSD.h"
#include "core/Colvar.h"
#include "tools/ColvarOutput.h"
#include "tools/Pbc.h"
#include "tools/SwitchingFunction.h"

#include <array>
#include <cmath>
#include <string>
#include <utility>
#include <vector>

namespace PLMD {
namespace secondarystructure {

/// The calculation of a single segment of ALPHARMSD_SCALAR, ANTIBETARMSD_SCALAR and PARABETARMSD_SCALAR. These functions
/// are used both by SecondaryStructureCount and by its GPU version in the openaccPTM plugin, so that the two give the same result.
/// A segment is given by the first atoms of its two halves of 15 atoms (strands): for a helix they are consecutive.
namespace sscount {

/// The number of atoms of a strand (N, CA, CB, C and O of three residues)
constexpr unsigned strandAtoms=15;
/// The position of the CA of the central residue in a strand, used for STRANDS_CUTOFF and to align the strands
constexpr unsigned centralCA=6;
/// The number of derivatives of the distance of a segment from a reference structure: its atoms and the virial
constexpr unsigned segmentDerivatives=3*maxAtomsInSegment+9;

inline Vector atomPosition( const double* x, const unsigned a ) {
  return Vector( x[3*a], x[3*a+1], x[3*a+2] );
}

/// The largest distance of an atom of the strand that starts at atom s from its central CA, with the strand made
/// whole as in segmentPositions (pbc is null with NOPBC)
inline double strandRadius( const double* x, const unsigned s, const Pbc* pbc ) {
  std::array<Vector,strandAtoms> pos;
  pos[0]=atomPosition( x, s );
  for(unsigned i=0; i+1<strandAtoms; ++i) {
    pos[i+1] = pbc ? pos[i]+pbc->distance( pos[i], atomPosition( x, s+i+1 ) ) : atomPosition( x, s+i+1 );
  }
  double r2=0.0;
  for(unsigned i=0; i<strandAtoms; ++i) {
    const double d2=delta( pos[centralCA], pos[i] ).modulo2();
    r2 = d2>r2 ? d2 : r2;
  }
  return sqrt(r2);
}

/// Is the segment made of the strands starting at atoms s1 and s2 computed? x are the coordinates of the atoms,
/// pbc is null with NOPBC. It is not if the distance between the central CA atoms of the strands is larger than
/// cutoff (if cutoff>0), as step(STRANDS_CUTOFF-x) in the shortcut, or, if radius is not null, if that distance minus
/// the radii of the strands (see strandRadius, radius[s] is the one of the strand that starts at atom s) is larger
/// than skipDistance: then the switching function of the segment is zero (see SecondaryStructureCount::setSkipDistance)
inline bool segmentActive( const double* x, const unsigned s1, const unsigned s2, const Pbc* pbc, const double cutoff,
                           const double* radius, const double skipDistance ) {
  const Vector p1=atomPosition( x, s1+centralCA );
  const Vector p2=atomPosition( x, s2+centralCA );
  const double d = ( pbc ? pbc->distance( p1, p2 ) : delta( p1, p2 ) ).modulo();
  if( cutoff>0 && d>cutoff ) {
    return false;
  }
  return !radius || d-radius[s1]-radius[s2]<=skipDistance;
}

/// The positions of the atoms of the segment made of the strands starting at atoms s1 and s2, as in
/// SecondaryStructureBase::performTask (pbc is null with NOPBC). With alignStrands (the sheets) the strands are made
/// whole separately, then the second is moved to the image where its central CA atom is closest to that of the first;
/// otherwise (the helices, where s2=s1+15) the segment is made whole as a single chain of bonded atoms.
inline void segmentPositions( const double* x, const unsigned s1, const unsigned s2, const Pbc* pbc, const bool alignStrands,
                              std::array<Vector,maxAtomsInSegment>& pos ) {
  for(unsigned k=0; k<strandAtoms; ++k) {
    pos[k]=atomPosition( x, s1+k );
    pos[strandAtoms+k]=atomPosition( x, s2+k );
  }
  if( pbc && !alignStrands ) {
    for(unsigned i=0; i+1<maxAtomsInSegment; ++i) {
      pos[i+1]=pos[i]+pbc->distance(pos[i],pos[i+1]);
    }
  } else if( pbc ) {
    for(unsigned start=0; start<maxAtomsInSegment; start+=strandAtoms) {
      for(unsigned i=start; i+1<start+strandAtoms; ++i) {
        pos[i+1]=pos[i]+pbc->distance(pos[i],pos[i+1]);
      }
    }
    const Vector shift = pbc->distance(pos[centralCA],pos[strandAtoms+centralCA]) + pos[centralCA] - pos[strandAtoms+centralCA];
    for(unsigned i=strandAtoms; i<maxAtomsInSegment; ++i) {
      pos[i]+=shift;
    }
  }
}

/// The switching function of the lowest distance of a segment from the reference structures (as LOWEST and LESS_THAN
/// in the shortcuts). It returns s(r) and sets g=ds/dr. If g is not zero, der[best*segmentDerivatives ...] contains
/// the derivatives of r with respect to the positions of the atoms of the segment and the virial; with several
/// structures the distances are first computed without derivatives, and the derivatives are computed only for the
/// lowest one, and only if g is not zero.
template <class I, class SW>
double segmentSwitch( const I& input, const unsigned nstructures, const SW& switchingFunction,
                      std::array<Vector,maxAtomsInSegment>& pos,
                      std::array<double,2*segmentDerivatives>& der, unsigned& best, double& g ) {
  std::array<double,2> dist{};
  ColvarOutput<double> output( View<double>(dist.data(),nstructures), segmentDerivatives, der.data() );
  const View<Vector> segment( pos.data(), maxAtomsInSegment );
  const bool single = (nstructures==1);
  for(unsigned n=0; n<nstructures; ++n) {
    I::calculateDistance( n, !single, input, segment, output );
  }
  best=0;
  for(unsigned n=1; n<nstructures; ++n) {
    if( dist[n]<dist[best] ) {
      best=n;
    }
  }
  const double r=dist[best];
  double df=0.0;
  const double s=switchingFunction.calculate( r, df );
  // df is (ds/dr)/r
  g=df*r;
  if( !single && g!=0.0 ) {
    I::calculateDistance( best, false, input, segment, output );
  }
  return s;
}

} // namespace sscount

/// The number of segments that are similar to an alpha helix or to a beta sheet, computed in a single action
/// (ALPHARMSD_SCALAR, ANTIBETARMSD_SCALAR and PARABETARMSD_SCALAR, see ALPHARMSD, ANTIBETARMSD and PARABETARMSD).
/// The openaccPTM plugin derives the GPU version from it.
class SecondaryStructureCount : public Colvar {
protected:
  /// ALPHARMSD: the segments are six consecutive residues, made whole as a chain
  bool alpha;
  bool parallel;
  bool nopbc;
  bool usedrmsd;
  double strands_cutoff;
  /// The pairs of strands whose central CA atoms are further apart than this plus the radii of the strands are
  /// skipped (see sscount::segmentActive); negative if no pair can be skipped in this way
  double skipDistance;
  /// The first atom of each strand, if pairs are skipped with skipDistance
  std::vector<unsigned> strands;
  /// The definition of the switching function, as given in SWITCH
  std::string switchInput;
  SwitchingFunction switchingFunction;
  /// The first atom of the two strands of each segment, in the list of atoms (for a helix the second is the first plus 15)
  std::vector<std::pair<unsigned,unsigned> > pairs;
  SecondaryStructureDRMSDInput<double> drmsd;
  SecondaryStructureRMSDInput rmsd;
  unsigned nstructures;
  /// Sets skipDistance from the reference structures and the switching function
  void setSkipDistance( const std::vector<std::vector<Vector> >& references, double bondlength );
  /// Sets the value and the derivatives from the sum over the segments (3*natoms derivatives and the virial)
  void setValueAndDerivatives( double value, const std::vector<double>& derivatives );
private:
  /// Buffers of the derivatives of each OpenMP thread
  std::vector<std::vector<double> > omp_derivatives;
  std::vector<double> derivatives;
  /// The radius of each strand (sscount::strandRadius), at the position of its first atom
  std::vector<double> radius;
public:
  static void registerKeywords( Keywords& keys );
  explicit SecondaryStructureCount(const ActionOptions&);
  void calculate() override;
};

}
}
#endif
