/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2014-2023 The plumed team
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
#include "core/ActionAtomistic.h"
#include "core/ActionPilot.h"
#include "core/ActionRegister.h"
#include "tools/Vector.h"
#include "tools/AtomNumber.h"
#include "tools/Tools.h"
#include "tools/Pbc.h"
#include "tools/PDB.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "core/GenericMolInfo.h"

#include <vector>
#include <array>
#include <algorithm>
#include <cmath>
#include <limits>
#include <tuple>
#include <unordered_map>
#include <unordered_set>

namespace PLMD {
namespace generic {

//+PLUMEDOC GENERIC WRAPAROUND
/*
Rebuild periodic boundary conditions around chosen atoms.

This action modifies the position of the atoms indicated by ATOMS by shifting them by lattice vectors so that they are
as close as possible to the atoms indicated by AROUND. More precisely, for every atom i
in the ATOMS list the following procedure is performed:
- The atom j among those in the AROUND list is searched that is closest to atom i.
- The atom i is replaced with its periodic image that is closest to atom j.

This action works similarly to [WHOLEMOLECULES](WHOLEMOLECULES.md) in that it replaces atoms coordinate. Notice that only
atoms specified with ATOMS are replaced, and that, at variance with [WHOLEMOLECULES](WHOLEMOLECULES.md),
the order in which atoms are specified is irrelevant when grouping and PAIR are not used.

This is often convenient at a post processing stage (using the driver), but sometime
it is required during the simulation if collective variables need atoms to be in a specific periodic image.

!!! caution "modifies stored positions

    This directive modifies the stored position at the precise moment it is executed. This means that only collective variables which are below it in
    the input script will see the corrected positions. As a general rule, put it at the top of the input file. Also, unless you know exactly what you are doing,
    leave the default stride (1), so that this action is performed at every MD step.

This action searches the reference list once per atom or group. Large AROUND
lists use a cell search; smaller lists are searched directly. Reconstructing
residue batches and computing their centers takes linear work in the number
of selected atoms. Large selections can still be expensive.
If you are using it to analyze a trajectory this is usually not a big problem. If you use it to
analyze a simulation on the fly, e.g. with [DUMPATOMS](DUMPATOMS.md) to store a properly wrapped trajectory,
consider using the STRIDE keyword here (with great care).


## Examples

This command instructs plumed to move all the ions to their periodic image that is as close as possible to
the rna group.

```plumed
rna: GROUP ATOMS=1-100
ions: GROUP ATOMS=101-110
# first make the rna molecule whole
WHOLEMOLECULES ENTITY0=rna
WRAPAROUND ATOMS=ions AROUND=rna
DUMPATOMS FILE=dump.xyz ATOMS=rna,ions
```

In case you want to do it during a simulation and you only care about wrapping the ions in
the `dump.xyz` file, you can use the following input:

```plumed
# add some restraint that do not require molecules to be whole:
a: TORSION ATOMS=1,2,10,11
RESTRAINT ARG=a AT=0.0 KAPPA=5


# then do the things that are required for dumping the trajectory
# notice that they are all done every 100 steps, so as not to
# unnecessarily overload the calculation

rna: GROUP ATOMS=1-100
ions: GROUP ATOMS=101-110
# first make the rna molecule whole
WHOLEMOLECULES ENTITY0=rna STRIDE=100
WRAPAROUND ATOMS=ions AROUND=rna STRIDE=100
DUMPATOMS FILE=dump.xyz ATOMS=rna,ions STRIDE=100
```

Notice that if the biased variable requires a molecule to be whole, you might have to put
the [WHOLEMOLECULES](WHOLEMOLECULES.md) command before computing that variable and leave the default STRIDE=1.

This command instructs plumed to center all atoms around the center of mass of a solute molecule.

```plumed
solute: GROUP ATOMS=1-100
all: GROUP ATOMS=1-1000
# center of the solute:
# notice that since plumed 2.2 this also works if the
# solute molecule is broken
com: COM ATOMS=solute
# notice that we wrap around a single atom. this should be fast
WRAPAROUND ATOMS=all AROUND=com
DUMPATOMS FILE=dump.xyz ATOMS=all
```

Notice that whereas [WHOLEMOLECULES](WHOLEMOLECULES.md) is designed to make molecules whole,
WRAPAROUND can easily break molecules. In the last example,
if solvent (atoms 101-1000) is made e.g. of water, then water
molecules could be broken by WRAPAROUND (hydrogen could end up
in an image and oxygen in another one).
One solution is to use [WHOLEMOLECULES](WHOLEMOLECULES.md) on _all_ the water molecules
after WRAPAROUND. This is tedious. A better solution is to use the
GROUPBY option which is going
to consider the atoms listed in ATOMS as a list of groups
each of size GROUPBY. The first atom of the group will be brought
close to the AROUND atoms. The following atoms of the group
will be just brought close to the first atom of the group.
Assuming that oxygen is the first atom of each water molecules,
in the following examples all the water oxygen atoms will be brought
close to the solute, and all the hydrogen atoms will be kept close
to their related oxygen.

```plumed
solute: GROUP ATOMS=1-100
water: GROUP ATOMS=101-1000
com: COM ATOMS=solute
# notice that we wrap around a single atom. this should be fast
WRAPAROUND ATOMS=solute AROUND=com
# notice that we wrap around a single atom. this should be fast
WRAPAROUND ATOMS=water AROUND=com GROUPBY=3
DUMPATOMS FILE=dump.xyz ATOMS=solute,water
```


## Grouping by residue

`GROUP_BY_RESIDUE` handles mixed residue sizes in a single selection, for example
three- and four-site waters and single-atom ions. A preceding MOLINFO supplies
residue membership only: its reference coordinates and box are ignored. Neither
WHOLE nor SMALL_RESIDUES is required.

Duplicate ATOMS entries are removed while retaining their first occurrence.
Each consecutive run of selected atoms of the same residue forms a batch. The
metadata distinguishes chain, residue number, insertion code, TER boundaries and
nonconsecutive occurrences of reused residue identifiers. Partial selections
include only the selected members. If one residue occurs in separate ATOMS
runs, those batches are wrapped independently.

At each active step, the selected atoms in every batch are reconstructed in
ATOMS order, as with WHOLEMOLECULES without a reference tree. WRAPAROUND computes
the geometric center of these current whole coordinates, finds its closest
periodic reference in AROUND, and translates the entire reconstructed batch by
one lattice vector. This preserves its reconstructed internal geometry.

```plumed
#SETTINGS MOLFILE=regtest/basic/rt-wraparound-residue-centers/metadata.pdb
MOLINFO STRUCTURE=regtest/basic/rt-wraparound-residue-centers/metadata.pdb
solute: GROUP ATOMS=1-3
solvent: GROUP ATOMS=4-14
WHOLEMOLECULES ENTITY0=solute
WRAPAROUND ATOMS=solvent AROUND=solute GROUP_BY_RESIDUE
DUMPATOMS ATOMS=solute,solvent FILE=wrapped.xyz
```

Use this option for residues whose selected members can be reconstructed by
minimum-image links in the supplied ATOMS order. Ordinary compact biomolecular
residues normally satisfy this condition. It is an explicit grouping request;
a material encoded as one large residue should use a reconstruction/grouping
appropriate to that material instead. A multi-residue solute can still be split
by independently wrapping its residues around a point. To keep it whole, first
use WHOLEMOLECULES on the solute and wrap solvent around the solute's atoms.

`PAIR` requires exactly one AROUND entry per batch, matched in batch order, and
retains repeated AROUND entries. A count mismatch is an error. The default
`GROUPBY=1`, including an explicit value of 1, is accepted with GROUP_BY_RESIDUE;
other GROUPBY values conflict. Existing fixed-size GROUPBY behavior is unchanged.

*/
//+ENDPLUMEDOC


class WrapAround:
  public ActionPilot,
  public ActionAtomistic {
  // cppcheck-suppress duplInheritedMember
  std::vector<Vector> refatoms;
  // distances from the atom being wrapped to all the reference atoms
  std::vector<Vector> distances;
  std::vector<std::pair<std::size_t,std::size_t> > p_atoms;
  std::vector<std::pair<std::size_t,std::size_t> > p_reference;
  unsigned groupby;
  bool pair_;
  bool byResidue=false;
  std::vector<unsigned> batchStarts, batchParents, batchChildren;
  std::vector<Vector> batchPositions;
  // cells of the reference atoms (orthorhombic boxes): the reference atoms of cell c are
  // cellatoms[cellstart[c]] ... cellatoms[cellstart[c+1]-1], in increasing order
  // (cells in scaled coordinates, so that the boxes need not be orthorhombic)
  std::array<unsigned,3> ncells{};
  // distance between the opposite faces of a cell
  std::array<double,3> cellwidth{};
  std::vector<unsigned> cellstart, cellatoms, cellstamp;
  unsigned stamp=0;
  // distance vectors of the reference atoms of one cell (triclinic boxes)
  std::vector<Vector> celldistances;
  void buildCells(const Pbc& pbc);
  unsigned cellOf(const Pbc& pbc, const Vector& x, std::array<int,3>& c, std::array<double,3>& u) const;
  int closestInCells(const Pbc& pbc, const Vector& first, bool orthorombic, const Vector& diag, const Vector& invdiag);
public:
  explicit WrapAround(const ActionOptions&ao);
  static void registerKeywords( Keywords& keys );
  bool actionHasForces() override {
    return false;
  }
  void calculate() override;
  void apply() override {}
};

PLUMED_REGISTER_ACTION(WrapAround,"WRAPAROUND")

void WrapAround::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys );
  ActionAtomistic::registerKeywords( keys );
  ActionPilot::registerKeywords( keys );
  keys.add("compulsory","STRIDE","1","the frequency with which molecules are reassembled.  Unless you are completely certain about what you are doing leave this set equal to 1!");
  keys.add("atoms","AROUND","the reference atoms: the atoms in ATOMS are moved to the periodic image that is closest to one of these atoms");
  keys.add("atoms","ATOMS","the atoms that should be wrapped");
  keys.add("compulsory","GROUPBY","1","group atoms so as not to break molecules");
  keys.addFlag("GROUP_BY_RESIDUE", false, "Group consecutive selected atoms of each residue, reconstruct them in ATOMS order and wrap their current geometric center. Requires MOLINFO metadata only");
  keys.addFlag("PAIR", false, "Pair each ATOMS group with one AROUND entry in group order");
}

WrapAround::WrapAround(const ActionOptions&ao):
  Action(ao),
  ActionPilot(ao),
  ActionAtomistic(ao),
  groupby(1),
  pair_(false) {
  std::vector<AtomNumber> atoms;
  parseAtomList("ATOMS",atoms);
  std::vector<AtomNumber> reference;
  parseAtomList("AROUND",reference);
  if(reference.empty()) {
    error("no reference atoms given with AROUND");
  }
  parse("GROUPBY",groupby);
  if(groupby<1) {
    error("GROUPBY should be at least 1");
  }
  parseFlag("PAIR", pair_);
  parseFlag("GROUP_BY_RESIDUE",byResidue);
  if(byResidue && groupby!=1) {
    error("GROUP_BY_RESIDUE cannot be combined with GROUPBY other than 1");
  }

  log.printf("  atoms in reference :");
  for(unsigned j=0; j<reference.size(); ++j) {
    log.printf(" %d",reference[j].serial() );
  }
  log.printf("\n");
  log.printf("  atoms to be wrapped :");
  for(unsigned j=0; j<atoms.size(); ++j) {
    log.printf(" %d",atoms[j].serial() );
  }
  log.printf("\n");
  if(groupby>1) {
    log<<"  atoms will be grouped by "<<groupby<<"\n";
  }
  if(pair_) {
    log.printf("  pairing atoms and references\n");
  }

  checkRead();

  if(byResidue) {
    const auto* info=plumed.getActionSet().selectLatest<GenericMolInfo*>(this);
    if(!info) {
      error("GROUP_BY_RESIDUE requires a preceding MOLINFO with structure metadata");
    }
    std::unordered_set<unsigned> seen;
    seen.reserve(atoms.size());
    std::vector<AtomNumber> distinctAtoms;
    for(const auto atom:atoms) {
      if(seen.insert(atom.index()).second) {
        distinctAtoms.push_back(atom);
      }
    }
    atoms.swap(distinctAtoms);
    // Distinguish each reference residue occurrence, including reused PDB numbers.
    // Coordinates and reference box are never consulted.
    const auto& pdb=info->getPDB();
    const auto& referenceAtoms=pdb.getAtomNumbers();
    const auto& blockEnds=pdb.getAtomBlockEnds();
    using ResidueKey=std::tuple<unsigned,std::string,int,std::string>;
    ResidueKey previous;
    unsigned block=0, residue=0;
    std::unordered_map<unsigned,unsigned> residueOfAtom;
    residueOfAtom.reserve(atoms.size());
    for(unsigned j=0; j<referenceAtoms.size(); ++j) {
      while(block<blockEnds.size() && j>=blockEnds[block]) {
        ++block;
      }
      const auto atom=referenceAtoms[j];
      const ResidueKey key(block,pdb.getChainID(atom),pdb.getResidueNumber(atom),pdb.getInsertionCode(atom));
      if(j>0 && key!=previous) {
        ++residue;
      }
      previous=key;
      if(seen.count(atom.index())>0) {
        residueOfAtom.emplace(atom.index(),residue);
      }
    }
    unsigned previousResidue=0;
    for(unsigned j=0; j<atoms.size(); ++j) {
      const auto found=residueOfAtom.find(atoms[j].index());
      if(found==residueOfAtom.end()) {
        error("GROUP_BY_RESIDUE: atom "+std::to_string(atoms[j].serial())+" is missing from the MOLINFO structure");
      }
      if(j==0 || found->second!=previousResidue) {
        batchStarts.push_back(j);
      } else {
        batchParents.push_back(j-1);
        batchChildren.push_back(j);
      }
      previousResidue=found->second;
    }
    if(pair_ && batchStarts.size()!=reference.size()) {
      error("with PAIR and GROUP_BY_RESIDUE you must have one AROUND entry per consecutive residue batch");
    }
    batchStarts.push_back(atoms.size());
    batchPositions.resize(atoms.size());
    log<<"  grouping consecutive atoms by residue and wrapping their current geometric centers\n";
  } else {
    if(atoms.size()%groupby!=0) {
      error("number of atoms should be a multiple of groupby option");
    }
    if(pair_ && atoms.size()!=reference.size()*groupby) {
      error("with PAIR you must have: #ATOMS = #AROUND * #GROUPBY");
    }
    if(!pair_ && groupby<=1) {
      Tools::removeDuplicates(atoms);
    }
  }
  if(!pair_) {
    Tools::removeDuplicates(reference);
  }

  // Group and PAIR lists need not be sorted; deduplicate only the request list.
  std::vector<AtomNumber> merged(atoms);
  merged.insert(merged.end(),reference.begin(),reference.end());
  p_atoms.resize( atoms.size() );
  for(unsigned i=0; i<atoms.size(); ++i) {
    p_atoms[i] = getValueIndices( atoms[i] );
  }
  refatoms.resize( reference.size() );
  p_reference.resize( reference.size() );
  for(unsigned i=0; i<reference.size(); ++i) {
    p_reference[i] = getValueIndices( reference[i] );
  }
  Tools::removeDuplicates(merged);
  requestAtoms(merged);
  doNotRetrieve();
  doNotForce();
}

// Cell of x, and position of x inside the cell along each direction (between 0 and 1)
unsigned WrapAround::cellOf(const Pbc& pbc, const Vector& x, std::array<int,3>& c, std::array<double,3>& u) const {
  const Vector scaled=pbc.realToScaled(x);
  for(unsigned k=0; k<3; ++k) {
    double f=scaled[k];
    f-=std::floor(f);
    c[k]=std::min(static_cast<int>(f*ncells[k]),static_cast<int>(ncells[k])-1);
    u[k]=std::min(std::max(f*ncells[k]-c[k],0.0),1.0);
  }
  return (c[0]*ncells[1]+c[1])*ncells[2]+c[2];
}

void WrapAround::buildCells(const Pbc& pbc) {
  const unsigned nref=refatoms.size();
  const Tensor& box=pbc.getBox();
  const Vector a0(box(0,0),box(0,1),box(0,2)), a1(box(1,0),box(1,1),box(1,2)), a2(box(2,0),box(2,1),box(2,2));
  const double volume=std::fabs(dotProduct(a0,crossProduct(a1,a2)));
  // distance between the faces of the box that are parallel to two lattice vectors
  const std::array<double,3> height= {volume/modulo(crossProduct(a1,a2)),volume/modulo(crossProduct(a2,a0)),volume/modulo(crossProduct(a0,a1))};
  // about one reference atom per cell (the fastest choice in tests with 64 to 2000 reference atoms)
  const double side=std::cbrt(volume/nref);
  for(unsigned k=0; k<3; ++k) {
    ncells[k]=std::max(1,std::min(256,static_cast<int>(height[k]/side)));
    cellwidth[k]=height[k]/ncells[k];
  }
  const unsigned ntot=ncells[0]*ncells[1]*ncells[2];
  cellstart.assign(ntot+1,0);
  std::vector<unsigned> cellof(nref);
  std::array<int,3> c;
  std::array<double,3> u;
  for(unsigned j=0; j<nref; ++j) {
    cellof[j]=cellOf(pbc,refatoms[j],c,u);
    cellstart[cellof[j]+1]++;
  }
  for(unsigned i=0; i<ntot; ++i) {
    cellstart[i+1]+=cellstart[i];
  }
  cellatoms.resize(nref);
  std::vector<unsigned> fill(cellstart.begin(),cellstart.end()-1);
  for(unsigned j=0; j<nref; ++j) {
    cellatoms[fill[cellof[j]]++]=j;
  }
  if(cellstamp.size()!=ntot) {
    cellstamp.assign(ntot,0);
    stamp=0;
  }
}

// The reference atom closest to first, with the same distances and the same choice
// among equal distances (the lowest index) as the loop over all the reference atoms.
// The cells are visited in shells of increasing distance, and the search stops when
// the cells that are left are farther than the closest atom found.
int WrapAround::closestInCells(const Pbc& pbc, const Vector& first, bool orthorombic, const Vector& diag, const Vector& invdiag) {
  std::array<int,3> c0;
  std::array<double,3> u;
  cellOf(pbc,first,c0,u);
  // distance of the atom from the closest face of its cell along each direction
  std::array<double,3> dface;
  int rmax=0;
  for(unsigned k=0; k<3; ++k) {
    dface[k]=std::min(u[k],1.0-u[k])*cellwidth[k];
    rmax=std::max(rmax,static_cast<int>(ncells[k]/2)+1);
  }
  if(++stamp==0) {
    std::fill(cellstamp.begin(),cellstamp.end(),0);
    stamp=1;
  }
  double mindist2=std::numeric_limits<double>::max();
  int closest=-1;
  auto consider=[&](unsigned j, double distance2) {
    if(distance2<mindist2 || (distance2==mindist2 && static_cast<int>(j)<closest)) {
      mindist2=distance2;
      closest=j;
    }
  };
  for(int r=0; r<=rmax; ++r) {
    for(int a=-r; a<=r; ++a) {
      const int ia=((c0[0]+a)%static_cast<int>(ncells[0])+ncells[0])%ncells[0];
      for(int b=-r; b<=r; ++b) {
        const int ib=((c0[1]+b)%static_cast<int>(ncells[1])+ncells[1])%ncells[1];
        // only the surface of the cube of side 2r+1: the inside was visited before
        const bool face=(std::abs(a)==r || std::abs(b)==r);
        for(int e=-r; e<=r; e+=(face || r==0 ? 1 : 2*r)) {
          const int ie=((c0[2]+e)%static_cast<int>(ncells[2])+ncells[2])%ncells[2];
          const unsigned cell=(ia*ncells[1]+ib)*ncells[2]+ie;
          if(cellstamp[cell]==stamp) {
            continue;
          }
          cellstamp[cell]=stamp;
          if(orthorombic) {
            // same operations as in the loop over all the reference atoms
            for(unsigned m=cellstart[cell]; m<cellstart[cell+1]; ++m) {
              const unsigned j=cellatoms[m];
              double distance2=0.0;
              for(unsigned k=0; k<3; ++k) {
                const double d=Tools::pbc((refatoms[j][k]-first[k])*invdiag[k])*diag[k];
                distance2+=d*d;
              }
              consider(j,distance2);
            }
          } else {
            // Pbc::apply, as in the loop over all the reference atoms
            const unsigned n=cellstart[cell+1]-cellstart[cell];
            if(n==0) {
              continue;
            }
            celldistances.resize(n);
            for(unsigned m=0; m<n; ++m) {
              celldistances[m]=delta(first,refatoms[cellatoms[cellstart[cell]+m]]);
            }
            pbc.apply(VectorView(&celldistances[0][0],n));
            for(unsigned m=0; m<n; ++m) {
              consider(cellatoms[cellstart[cell]+m],modulo2(celldistances[m]));
            }
          }
        }
      }
    }
    // A cell that is not visited yet differs by more than r cells along a direction k in which
    // the visited cells do not cover the whole box: it is at least r*cellwidth[k]+dface[k] away
    // (also for its periodic images, which differ by whole boxes along the same directions)
    double bound=std::numeric_limits<double>::max();
    for(unsigned k=0; k<3; ++k) {
      if(2*r+1<static_cast<int>(ncells[k])) {
        bound=std::min(bound,r*cellwidth[k]+dface[k]);
      }
    }
    if(bound==std::numeric_limits<double>::max()) {
      break; // all the cells were visited
    }
    // (with a margin for rounding, so that atoms at the same distance are never split)
    if(closest>=0 && mindist2<bound*bound*(1.0-1e-9)) {
      break;
    }
  }
  return closest;
}

void WrapAround::calculate() {
  const unsigned nref=p_reference.size();
  for(unsigned j=0; j<nref; ++j) {
    refatoms[j] = getGlobalPosition(p_reference[j]);
  }
  if(!pair_) {
    distances.resize(nref);
  }
  const Pbc& pbc=getPbc();
  if(byResidue) {
    for(unsigned j=0; j<p_atoms.size(); ++j) {
      batchPositions[j]=getGlobalPosition(p_atoms[j]);
    }
    // One forest reconstruction for all batches: reuse the optimized integer
    // image shifts and box-dependent setup, without per-residue allocations.
    pbc.makeWhole(batchPositions,batchParents.size(),batchParents.data(),batchChildren.data());
  }
  const bool orthorombic=pbc.isOrthorombic();
  const Vector diag(pbc.getBox()(0,0),pbc.getBox()(1,1),pbc.getBox()(2,2));
  const Vector invdiag(pbc.getInvBox()(0,0),pbc.getInvBox()(1,1),pbc.getInvBox()(2,2));
  // with few reference atoms the loop over all of them is faster than the cells
  const unsigned cellsMinReference=128;
  for(unsigned i=0, group=0; i<p_atoms.size(); ++group) {
    const unsigned end=byResidue ? batchStarts[group+1] : i+groupby;
    Vector second, first;
    if(byResidue) {
      const Vector origin=batchPositions[i];
      Vector sum;
      for(unsigned j=i+1; j<end; ++j) {
        sum+=batchPositions[j]-origin;
      }
      first=origin+sum/(end-i);
    } else {
      first=getGlobalPosition(p_atoms[i]);
    }
    int closest=-1;
    if(pair_) {
      closest = group;
    } else {
      double mindist2=std::numeric_limits<double>::max();
      if(pbc.isSet() && nref>=cellsMinReference) {
        if(i==0) {
          buildCells(pbc);
        }
        closest=closestInCells(pbc,first,orthorombic,diag,invdiag);
      } else if(orthorombic) {
        // same operations as Pbc::distance, written out so that the loop is inlined
        for(unsigned j=0; j<nref; ++j) {
          double distance2=0.0;
          for(unsigned k=0; k<3; ++k) {
            const double d=Tools::pbc((refatoms[j][k]-first[k])*invdiag[k])*diag[k];
            distance2+=d*d;
          }
          if(distance2<mindist2) {
            mindist2=distance2;
            closest=j;
          }
        }
      } else {
        // all the distances to the reference atoms are computed in one call,
        // the result is the same as calling pbcDistance(first,refatoms[j]) for each j
        for(unsigned j=0; j<nref; ++j) {
          distances[j]=delta(first,refatoms[j]);
        }
        pbc.apply(VectorView(&distances[0][0],nref));
        for(unsigned j=0; j<nref; ++j) {
          const double distance2=modulo2(distances[j]);
          if(distance2<mindist2) {
            mindist2=distance2;
            closest=j;
          }
        }
      }
      plumed_massert(closest>=0,"closest not found");
    }
    second=refatoms[closest];
    if(byResidue) {
      // All atoms receive the same final lattice translation. Do not perform
      // another minimum-image reconstruction around the translated center.
      Vector shift;
      if(pbc.isSet()) {
        const Vector difference=first-second;
        const Vector correction=pbcDistance(second,first)-difference;
        if(orthorombic) {
          for(unsigned k=0; k<3; ++k) {
            shift[k]=std::nearbyint(correction[k]*invdiag[k])*diag[k];
          }
        } else {
          Vector lattice=pbc.realToScaled(correction);
          for(unsigned k=0; k<3; ++k) {
            lattice[k]=std::nearbyint(lattice[k]);
          }
          shift=pbc.scaledToReal(lattice);
        }
      }
      for(unsigned j=i; j<end; ++j) {
        const Vector moved=batchPositions[j]+shift;
        const Vector original=getGlobalPosition(p_atoms[j]);
        if(moved[0]!=original[0] || moved[1]!=original[1] || moved[2]!=original[2]) {
          setGlobalPosition(p_atoms[j],moved);
        }
      }
    } else {
      first=second+pbcDistance(second,first);
      setGlobalPosition(p_atoms[i],first);
      for(unsigned j=i+1; j<end; ++j) {
        second=getGlobalPosition(p_atoms[j]);
        setGlobalPosition(p_atoms[j],first+pbcDistance(first,second));
      }
    }
    i=end;
  }
}



}

}
