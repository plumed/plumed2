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
#include "PDB.h"
#include "Exception.h"
#ifdef __PLUMED_HAS_LIBCHEMFILES
#include <chemfiles.hpp>
#include <algorithm>
#include <limits>
#endif

// Reader of GROMACS TPR files for the PDB class, with chemfiles (the documentation is in PDB.cpp).

namespace PLMD {

#ifdef __PLUMED_HAS_LIBCHEMFILES
namespace {
/// chain ids A, B, ... Z, AA, AB, ...
std::string chainName(std::size_t n) {
  std::string s;
  for(n++; n>0; n/=26) {
    n--;
    s.insert(s.begin(),static_cast<char>('A'+n%26));
  }
  return s;
}
}
#endif

void PDB::readTpr(const std::string& file,double scale) {
#ifdef __PLUMED_HAS_LIBCHEMFILES
  // the TPR format was added in chemfiles 0.10.4
  bool hasTpr=false;
  for(const auto & metadata : chemfiles::formats_list()) {
    hasTpr=hasTpr || metadata.get().name==std::string("TPR");
  }
  if(!hasTpr) {
    plumed_merror("cannot read the TPR file " + file + ": this needs chemfiles 0.10.4 or later");
  }
  chemfiles::Frame frame;
  try {
    chemfiles::Trajectory trajectory(file,'r',"TPR");
    frame=trajectory.read();
  } catch(const std::exception & e) {
    plumed_merror("cannot read the TPR file " + file + " with chemfiles: " + e.what());
  }
  const std::size_t natoms=frame.size();
  if(natoms==0) {
    plumed_merror("the TPR file " + file + " has no atoms");
  }
  if(natoms>std::size_t(std::numeric_limits<int>::max())) {
    plumed_merror("the TPR file " + file + " has too many atoms");
  }
  const auto & topology=frame.topology();
  const auto & residues=topology.residues();

  // the residues in the order of their atoms
  std::vector<std::size_t> residueOfAtom(natoms,residues.size());
  for(std::size_t r=0; r<residues.size(); r++) {
    for(const auto a : residues[r]) {
      residueOfAtom[a]=r;
    }
  }
  std::vector<std::size_t> sequence;
  std::vector<std::size_t> position(residues.size(),0);
  for(std::size_t i=0; i<natoms; i++) {
    const std::size_t r=residueOfAtom[i];
    if(r==residues.size()) {
      plumed_merror("the atom " + std::to_string(i+1) + " of the TPR file " + file + " is not in a residue");
    }
    if(sequence.empty() || sequence.back()!=r) {
      position[r]=sequence.size();
      sequence.push_back(r);
    }
  }

  // the chains are the sequences of residues bonded to the next one (protein and nucleic acid chains...):
  // as gmx editconf does, they are named A, B, ... and keep their residue numbers, while the other
  // residues (water, ions...) have no chain and are numbered after the largest residue number of the chains
  std::vector<bool> bondedToNext(sequence.size(),false);
  for(const auto & bond : topology.bonds()) {
    const std::size_t p1=std::min(position[residueOfAtom[bond[0]]],position[residueOfAtom[bond[1]]]);
    const std::size_t p2=std::max(position[residueOfAtom[bond[0]]],position[residueOfAtom[bond[1]]]);
    if(p2==p1+1) {
      bondedToNext[p1]=true;
    }
  }
  std::vector<int> number(sequence.size(),0);
  std::vector<bool> inChain(sequence.size(),false);
  int lastNumber=0;
  for(std::size_t k=0; k<sequence.size(); k++) {
    const auto & res=residues[sequence[k]];
    if(res.id()) {
      const auto id=*res.id();
      if(id<std::numeric_limits<int>::min() || id>std::numeric_limits<int>::max()) {
        plumed_merror("the TPR file " + file + " has the residue number " + std::to_string(id) + ", which is too large");
      }
      number[k]=static_cast<int>(id);
    }
    inChain[k]=bondedToNext[k] || (k>0 && bondedToNext[k-1]);
    if(inChain[k]) {
      lastNumber=std::max(lastNumber,number[k]);
    }
  }
  std::vector<std::string> chainOf(sequence.size()," ");
  std::size_t nchains=0;
  for(std::size_t k=0; k<sequence.size(); k++) {
    if(inChain[k]) {
      chainOf[k]=(k>0 && bondedToNext[k-1] ? chainOf[k-1] : chainName(nchains++));
    } else {
      if(lastNumber==std::numeric_limits<int>::max()) {
        plumed_merror("the TPR file " + file + " has too many residues");
      }
      number[k]=++lastNumber;
    }
  }

  const auto coordinates=frame.positions();
  for(std::size_t i=0; i<natoms; i++) {
    const std::size_t k=position[residueOfAtom[i]];
    const auto & res=residues[sequence[k]];
    const auto & atom=frame[i];
    std::string resname=res.name();
    // as in the columns 18-20 of a PDB file, where short names are right-justified
    if(resname.length()<3) {
      resname=std::string(3-resname.length(),' ')+resname;
    }
    std::string icode;
    const auto ic=res.get<chemfiles::Property::STRING>("insertion_code");
    if(ic && *ic!=" ") {
      icode=*ic;
    }
    AtomNumber a;
    a.setIndex(i);
    const Vector p(coordinates[i][0],coordinates[i][1],coordinates[i][2]);
    addAtom(a,atom.name(),resname,number[k],chainOf[k],icode,(atom.atomic_number() ? atom.type() : ""),p*scale,1.0,1.0);
  }
  block_ends.push_back( size() );

  // the simulation box, as the CRYST1 record written by gmx editconf
  const auto & cell=frame.cell();
  if(cell.shape()!=chemfiles::UnitCell::INFINITE) {
    const auto lengths=cell.lengths();
    const auto angles=cell.angles();
    BoxXYZ=Vector(lengths[0],lengths[1],lengths[2]);
    BoxXYZ*=scale;
    BoxABG=Vector(angles[0],angles[1],angles[2]);
    setBoxFromLengthsAndAngles();
    spaceGroup="P 1";
  }
#else
  (void) scale;
  plumed_merror("cannot read the TPR file " + file + ": PLUMED must be configured with --enable-libchemfiles to read TPR files");
#endif
}

}
