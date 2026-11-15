/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2023 The plumed team
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
#include "Colvar.h"
#include "tools/Communicator.h"
#include "tools/PDB.h"
#include "core/ActionRegister.h"
#include "tools/SwitchingFunction.h"
#include "tools/OpenMP.h"

#include <algorithm>
#include <array>
#include <cstdlib>
#include <limits>
#include <sstream>

namespace PLMD {
namespace colvar {

//+PLUMEDOC COLVAR CONTACTMAP
/*
Calculate the distances between a number of pairs of atoms and transform each distance by a switching function.

This CV calculates a series of distances between pairs of atoms and transforms them by a switching function.
The transformed distance can be compared with a reference value in order to calculate the squared distance
between two contact maps. Each distance can also be weighted for a given value. CONTACTMAP can be used together
with [FUNCPATHMSD](FUNCPATHMSD.md) to define a path in contactmap space as was done in the 2008 paper by Bonomi that
is cited below.

If $s_i$ is the switching function of the $i$-th contact, $r_i$ its distance, $w_i$ its weight and $s_i^0$ its reference
value, CONTACTMAP computes:

- without SUM and CMDIST, one component for each contact, $w_i (s_i(r_i) - s_i^0)$;
- with SUM, $\sum_i w_i (s_i(r_i) - s_i^0)$;
- with CMDIST, $\sum_i w_i^2 (s_i(r_i) - s_i^0)^2$.

The weights are 1 and the reference values 0 unless they are given with WEIGHT and REFERENCE.
With NORMALIZE_WEIGHTS the weights (the default ones or those given with WEIGHT) are divided by their sum, so that
they sum to one: with SUM the CV is then the weighted average of the contacts, and with CMDIST it is
$\sum_i w_i^2 (s_i(r_i) - s_i^0)^2 / (\sum_j w_j)^2$.
The individual contacts can be accessed as components
named `cm.contact-1`, `cm.contact-2`, etc, assuming that the label of the CONTACTMAP is `cm`.

## Examples

The following example calculates switching functions based on the distances between atoms
1 and 2, 3 and 4, 4 and 5 and 5 and 6. The values of these four switching functions are then output
to a file named colvar.

```plumed
f1: CONTACTMAP ATOMS1=1,2 ATOMS2=3,4 ATOMS3=4,5 ATOMS4=5,6 SWITCH={RATIONAL R_0=1.5}
PRINT ARG=f1.* FILE=colvar
```

The following example calculates the difference of the current contact map with respect
to a reference provided. In this case REFERENCE is the fraction of contact that is formed
(i.e. the distance between two atoms transformed with the SWITCH), while R_0 is the contact
distance. WEIGHT gives the relative weight of each contact to the final distance measure.

```plumed
cmap: CONTACTMAP ...
  ATOMS1=1,2 REFERENCE1=0.1 WEIGHT1=0.5
  ATOMS2=3,4 REFERENCE2=0.5 WEIGHT2=1.0
  ATOMS3=4,5 REFERENCE3=0.25 WEIGHT3=1.0
  ATOMS4=5,6 REFERENCE4=0.0 WEIGHT4=0.5
  SWITCH={RATIONAL R_0=1.5}
  CMDIST
...

PRINT ARG=cmap FILE=colvar
```

The reference contact map can also be taken from a reference structure, by giving a PDB file
(or any other structure file that [MOLINFO](MOLINFO.md) can read, e.g. mmCIF) to REFERENCE.
The reference value of each contact is then its switching function evaluated at the distance between the two atoms
in the reference structure. The atoms are identified by their serial numbers in the file,
the distances are computed without periodic boundary conditions (so the molecule should be whole in the file)
and the coordinates are read in Angstrom (or in the natural units if those are used).

```plumed
#SETTINGS INPUTFILES=regtest/basic/rt-contactmap-pdb/native.pdb
cmap: CONTACTMAP ...
  ATOMS1=1,67 ATOMS2=1,68 ATOMS3=5,61 ATOMS4=15,57
  SWITCH={RATIONAL R_0=0.6}
  REFERENCE=regtest/basic/rt-contactmap-pdb/native.pdb
  CMDIST
...

PRINT ARG=cmap FILE=colvar
```

The next example calculates the fraction of native contacts (Q)
for Trp-cage mini-protein. R_0 is the distance at which the switch function is guaranteed to
be 1.0 – it doesn't really matter for Q and  should be something very small, like 1 A.
REF is the reference distance for the contact, e.g. the distance from a crystal structure.
LAMBDA is the tolerance for the distance – if set to 1.0, the contact would have to have exactly
the reference value to be formed; instead for lambda values of 1.5–1.8 are usually used to allow some slack.
BETA is the softness of the switch function, default is 50 nm$^{-1}$.
WEIGHT is the 1/(number of contacts) giving equal weight to each contact.

When using native contact Q switch function, please cite the 2013 paper from Best from the bibliography below

```plumed
# The full (much-longer) example available in regtest/basic/rt72/

cmap: CONTACTMAP ...
   ATOMS1=1,67 SWITCH1={Q R_0=0.01 BETA=50.0 LAMBDA=1.5 REF=0.4059} WEIGHT1=0.003597
   ATOMS2=1,68 SWITCH2={Q R_0=0.01 BETA=50.0 LAMBDA=1.5 REF=0.4039} WEIGHT2=0.003597
   ATOMS3=1,69 SWITCH3={Q R_0=0.01 BETA=50.0 LAMBDA=1.5 REF=0.3215} WEIGHT3=0.003597
   ATOMS4=5,61 SWITCH4={Q R_0=0.01 BETA=50.0 LAMBDA=1.5 REF=0.4277} WEIGHT4=0.003597
   ATOMS5=5,67 SWITCH5={Q R_0=0.01 BETA=50.0 LAMBDA=1.5 REF=0.3851} WEIGHT5=0.003597
   ATOMS6=5,68 SWITCH6={Q R_0=0.01 BETA=50.0 LAMBDA=1.5 REF=0.3811} WEIGHT6=0.003597
   ATOMS7=5,69 SWITCH7={Q R_0=0.01 BETA=50.0 LAMBDA=1.5 REF=0.3133} WEIGHT7=0.003597
   SUM
...

PRINT ARG=cmap FILE=colvar
```

With a reference structure given to REFERENCE, the Q switching functions that do not have REF
take it from the distance in the reference structure, so the native contacts need no per-contact numbers.
Without CMDIST the reference structure only provides these distances: no reference value is subtracted,
so that the CV below, where NORMALIZE_WEIGHTS gives each contact the weight 1/8, is the fraction of native contacts, between 0 and 1.

```plumed
#SETTINGS INPUTFILES=regtest/basic/rt-contactmap-pdb/native.pdb
q: CONTACTMAP ...
   ATOMS1=1,67 ATOMS2=1,68 ATOMS3=1,69 ATOMS4=5,61
   ATOMS5=5,67 ATOMS6=5,68 ATOMS7=5,69 ATOMS8=11,68
   SWITCH={Q R_0=0.01 BETA=50.0 LAMBDA=1.5}
   REFERENCE=regtest/basic/rt-contactmap-pdb/native.pdb
   NORMALIZE_WEIGHTS
   SUM
...

PRINT ARG=q FILE=colvar
```

The contacts themselves can be taken from the reference structure: with CUTOFF instead of the ATOMS keywords,
the contacts are all the pairs of atoms of the reference structure that are closer than CUTOFF.
To use only some of the atoms (e.g. the heavy atoms, the $C_\alpha$ atoms or a single domain), give a reference structure
that contains only those atoms; the atoms are identified by their serial numbers, that must be those of the simulation.
CHAINS chooses the pairs of atoms in the same chain (INTRA), in different chains (INTER) or both (ALL, the default),
with the chains given by the chain identifiers of the reference structure.
Two atoms in the same chain form a contact only if their residue numbers differ by at least MIN_SEQ_SEPARATION (4 by default,
as in the 2013 paper by Best cited below); this does not apply to atoms in different chains.
The contacts that are found are listed in the log.

The input below calculates the fraction of native contacts between the heavy atoms of the reference structure
(that has two chains) that are closer than 0.45 nm, with residues at least 4 apart when in the same chain.

```plumed
#SETTINGS INPUTFILES=regtest/basic/rt-contactmap-cutoff/native-heavy.pdb
q: CONTACTMAP ...
   REFERENCE=regtest/basic/rt-contactmap-cutoff/native-heavy.pdb
   CUTOFF=0.45
   SWITCH={Q R_0=0.01 BETA=50.0 LAMBDA=1.5}
   NORMALIZE_WEIGHTS
   SUM
...

# only the contacts between the two chains
qinter: CONTACTMAP ...
   REFERENCE=regtest/basic/rt-contactmap-cutoff/native-heavy.pdb
   CUTOFF=0.45
   CHAINS=INTER
   SWITCH={Q R_0=0.01 BETA=50.0 LAMBDA=1.5}
   NORMALIZE_WEIGHTS
   SUM
...

PRINT ARG=q,qinter FILE=colvar
```

*/
//+ENDPLUMEDOC

class ContactMap : public Colvar {
private:
  bool pbc, serial, docomp, dosum, docmdist;
  std::vector<SwitchingFunction> sfs;
  std::vector<double> reference, weight;
/// derivative of each contact with respect to its second atom (the first one has the opposite)
  std::vector<Vector> deriv;
/// minimum number of contacts per OpenMP thread (measured: below these the threads are slower);
/// a component costs more than a term of the sum, since it has its own array of derivatives
  static constexpr unsigned minContactsPerThread=300;
  static constexpr unsigned minComponentsPerThread=100;
public:
  static void registerKeywords( Keywords& keys );
  explicit ContactMap(const ActionOptions&);
// active methods:
  void calculate() override;
  void checkFieldsAllowed() override {}
};

PLUMED_REGISTER_ACTION(ContactMap,"CONTACTMAP")

namespace {
/// With a reference structure, a native-contact (Q) switching function without REF
/// takes the distance in the reference structure. Returns true if REF was added.
bool addNativeDistance(std::string& sw, const double r) {
  const auto words=Tools::getWords(sw);
  if(words.empty() || words[0]!="Q") {
    return false;
  }
  for(const auto & w : words) {
    if(w.compare(0,4,"REF=")==0) {
      return false;
    }
  }
  std::ostringstream ostr;
  ostr.precision(std::numeric_limits<double>::max_digits10);
  ostr<<r;
  sw+=" REF="+ostr.str();
  return true;
}
/// The pairs of atoms of a reference structure closer than cutoff, sorted by atom numbers.
/// chains is ALL, INTRA (pairs in the same chain) or INTER (pairs in different chains);
/// the atoms of a pair in the same chain must have residue numbers at least minsep apart.
std::vector<std::pair<AtomNumber,AtomNumber>> contactsFromStructure(const PDB& pdb, const double cutoff,
const std::string& chains, const int minsep) {
  const auto & atoms=pdb.getAtomNumbers();
  const auto & pos=pdb.getPositions();
  const unsigned natoms=atoms.size();
  std::vector<std::string> chain(natoms);
  std::vector<int> resnum(natoms);
  for(unsigned i=0; i<natoms; ++i) {
    chain[i]=pdb.getChainID(atoms[i]);
    resnum[i]=pdb.getResidueNumber(atoms[i]);
  }
  // with the atoms sorted along x, the search of the partners of an atom stops at the first one farther than cutoff along x
  std::vector<unsigned> order(natoms);
  for(unsigned i=0; i<natoms; ++i) {
    order[i]=i;
  }
  std::sort(order.begin(),order.end(),[&](unsigned i,unsigned j) {
    return pos[i][0]<pos[j][0];
  });
  const double cutoff2=cutoff*cutoff;
  std::vector<std::pair<AtomNumber,AtomNumber>> pairs;
  for(unsigned a=0; a<natoms; ++a) {
    const unsigned i=order[a];
    for(unsigned b=a+1; b<natoms; ++b) {
      const unsigned j=order[b];
      if( pos[j][0]-pos[i][0]>=cutoff ) {
        break;
      }
      if( delta(pos[i],pos[j]).modulo2()>=cutoff2 ) {
        continue;
      }
      const bool samechain=( chain[i]==chain[j] );
      if( (samechain && chains=="INTER") || (!samechain && chains=="INTRA") ) {
        continue;
      }
      if( samechain && std::abs(resnum[i]-resnum[j])<minsep ) {
        continue;
      }
      pairs.emplace_back(std::min(atoms[i],atoms[j]),std::max(atoms[i],atoms[j]));
    }
  }
  std::sort(pairs.begin(),pairs.end());
  return pairs;
}
}

void ContactMap::registerKeywords( Keywords& keys ) {
  Colvar::registerKeywords( keys );
  keys.add("numbered","ATOMS","the atoms involved in each of the contacts you wish to calculate. "
           "Keywords like ATOMS1, ATOMS2, ATOMS3,... should be listed and one contact will be "
           "calculated for each ATOM keyword you specify.");
  keys.reset_style("ATOMS","atoms");
  keys.add("numbered","SWITCH","The switching functions to use for each of the contacts in your map. "
           "You can either specify a global switching function using SWITCH or one "
           "switching function for each contact. Details of the various switching "
           "functions you can use are provided in the documentation for LESS_THAN.");
  keys.add("numbered","REFERENCE","A reference value for a given contact, by default is 0.0 "
           "You can either specify a global reference value using REFERENCE or one "
           "reference value for each contact. REFERENCE can also be a PDB file with a reference structure: "
           "the reference values are then computed from the distances in it (used with CMDIST), and "
           "the Q switching functions without REF take REF from them. With CUTOFF the contacts are also taken from the reference structure.");
  keys.add("optional","CUTOFF","instead of the ATOMS keywords, use as contacts all the pairs of atoms that are closer than this distance "
           "in the reference structure given with REFERENCE");
  keys.add("optional","CHAINS","with CUTOFF, the pairs of atoms to use: ALL (default), INTRA (only pairs in the same chain) "
           "or INTER (only pairs in different chains)");
  keys.add("optional","MIN_SEQ_SEPARATION","with CUTOFF, the minimum difference between the residue numbers of the two atoms "
           "of a pair in the same chain (default 4); it does not apply to pairs in different chains");
  keys.add("numbered","WEIGHT","A weight value for a given contact, by default is 1.0 "
           "You can either specify a global weight value using WEIGHT or one "
           "weight value for each contact.");
  keys.addFlag("NORMALIZE_WEIGHTS",false,"divide the weights by their sum, so that they sum to one "
               "(with SUM and no other weight, the CV is the average of the contacts)");
  keys.reset_style("SWITCH","compulsory");
  keys.linkActionInDocs("SWITCH","LESS_THAN");
  keys.addFlag("SUM",false,"calculate the sum of all the contacts in the input");
  keys.addFlag("CMDIST",false,"calculate the distance with respect to the provided reference contact map");
  keys.addFlag("SERIAL",false,"Perform the calculation in serial - for debug purpose");
  keys.addOutputComponent("contact","default","scalar","By not using SUM or CMDIST each contact will be stored in a component");
  keys.setValueDescription("scalar","the sum of the contacts (SUM) or the squared distance from the reference contact map (CMDIST)");
  keys.addDOI("10.1021/ja803652f");
  keys.addDOI("10.1073/pnas.1311599110");
}

ContactMap::ContactMap(const ActionOptions&ao):
  PLUMED_COLVAR_INIT(ao),
  pbc(true),
  serial(false),
  docomp(true),
  dosum(false),
  docmdist(false) {
  parseFlag("SERIAL",serial);
  parseFlag("SUM",dosum);
  parseFlag("CMDIST",docmdist);
  if(docmdist==true&&dosum==true) {
    error("You cannot use SUM and CMDIST together");
  }
  bool nopbc=!pbc;
  parseFlag("NOPBC",nopbc);
  pbc=!nopbc;

  // Read in the atoms
  std::vector<AtomNumber> t, ga_lista, gb_lista;
  for(int i=1;; ++i ) {
    parseAtomList("ATOMS", i, t );
    if( t.empty() ) {
      break;
    }

    std::string num;
    Tools::convert(i,num);
    if( t.size()!=2 ) {
      error("ATOMS" + num + " keyword has the wrong number of atoms");
    }
    if( t[0]==t[1] ) {
      error("ATOMS" + num + " keyword contains the same atom twice");
    }
    ga_lista.push_back(t[0]);
    gb_lista.push_back(t[1]);
    t.resize(0);
  }

  // The reference: a single value, or a reference structure
  std::string refstr;
  parse("REFERENCE",refstr);
  double ref0=0.0;
  const bool globalref=( refstr.length()>0 && Tools::convertNoexcept(refstr,ref0) );
  PDB pdb;
  const bool refstructure=( refstr.length()>0 && !globalref );
  if( refstructure ) {
    if( !pdb.read(refstr,usingNaturalUnits(),0.1/getUnits().getLength()) ) {
      error("the value " + refstr + " given to REFERENCE is neither a number nor a readable PDB file");
    }
    log.printf("  reference structure read from file %s\n", refstr.c_str());
  }

  // The contacts can also be the pairs of atoms closer than CUTOFF in the reference structure
  std::string cutstr;
  parse("CUTOFF",cutstr);
  const bool usecutoff=( cutstr.length()>0 );
  double cutoff=0.0;
  if( usecutoff && ( !Tools::convertNoexcept(cutstr,cutoff) || !(cutoff>0.0) ) ) {
    error("CUTOFF should be a positive distance");
  }
  std::string chains;
  parse("CHAINS",chains);
  std::string minsepstr;
  parse("MIN_SEQ_SEPARATION",minsepstr);
  int minsep=4;
  if( minsepstr.length()>0 && ( !Tools::convertNoexcept(minsepstr,minsep) || minsep<0 ) ) {
    error("MIN_SEQ_SEPARATION should be a non-negative integer");
  }
  if( usecutoff ) {
    if( !ga_lista.empty() ) {
      error("use either ATOMS1, ATOMS2, ... or CUTOFF to choose the contacts");
    }
    if( !refstructure ) {
      error("CUTOFF needs a reference structure given with REFERENCE");
    }
    if( chains.length()==0 ) {
      chains="ALL";
    }
    if( chains!="ALL" && chains!="INTRA" && chains!="INTER" ) {
      error("CHAINS should be ALL, INTRA or INTER");
    }
    log.printf("  contacts are the pairs of atoms of the reference structure closer than %f", cutoff);
    if( chains!="INTER" ) {
      log.printf(", in the same chain with residue numbers at least %d apart", minsep);
    }
    log.printf(" (CHAINS=%s)\n", chains.c_str());
    for(const auto & c : contactsFromStructure(pdb,cutoff,chains,minsep)) {
      ga_lista.push_back(c.first);
      gb_lista.push_back(c.second);
    }
    if( ga_lista.empty() ) {
      error("no pairs of atoms in the reference structure " + refstr + " satisfy CUTOFF, CHAINS and MIN_SEQ_SEPARATION");
    }
  } else if( chains.length()>0 || minsepstr.length()>0 ) {
    error("CHAINS and MIN_SEQ_SEPARATION can only be used with CUTOFF");
  }
  if( ga_lista.empty() ) {
    error("no contacts have been specified, use the ATOMS1, ATOMS2, ... keywords");
  }
  const unsigned ncontacts=ga_lista.size();
  std::string snc;
  Tools::convert(ncontacts,snc);
  if( usecutoff ) {
    log.printf("  %s contacts found in the reference structure\n", snc.c_str());
  }

  // Add a value for each contact
  if(!dosum&&!docmdist) {
    for(unsigned i=0; i<ncontacts; ++i) {
      std::string num;
      Tools::convert(i+1,num);
      addComponentWithDerivatives("contact-"+num);
      componentIsNotPeriodic("contact-"+num);
    }
  }

  // Read numbered keywords KEY1, KEY2, ...: either none of them or one for each contact
  auto parseAllNumbered=[&](const std::string& key, std::vector<double>& v) {
    unsigned nread=0;
    for(unsigned i=0; i<ncontacts; ++i) {
      if( parseNumbered( key, i+1, v[i] ) ) {
        nread++;
      }
    }
    if( nread>0 && nread<ncontacts ) {
      std::string sn;
      Tools::convert(nread,sn);
      error("found " + sn + " " + key + " keywords for " + snc + " contacts: use either " + key +
            " or one numbered " + key + " keyword for each contact");
    }
    return nread>0;
  };

  // Read in the reference values: one value for each contact, a single value or from the reference structure
  reference.assign(ncontacts, 0.);
  const bool numberedref=parseAllNumbered("REFERENCE",reference);
  if( numberedref && refstr.length()>0 ) {
    error("use either REFERENCE or REFERENCE1, REFERENCE2, ...");
  }
  const bool hasref=( numberedref || refstr.length()>0 );
  if( globalref ) {
    reference.assign(ncontacts, ref0);
  }
  std::vector<double> refdist; // the distances in the reference structure
  if( refstructure ) {
    refdist.resize(ncontacts);
    for(unsigned i=0; i<ncontacts; ++i) {
      for(const AtomNumber a : std::array<AtomNumber,2> {ga_lista[i], gb_lista[i]}) {
        if( !pdb.checkForAtom(a) ) {
          std::string sa;
          Tools::convert(a.serial(),sa);
          error("atom " + sa + " is not in the reference structure " + refstr);
        }
      }
      refdist[i]=delta(pdb.getPosition(ga_lista[i]),pdb.getPosition(gb_lista[i])).modulo();
    }
  }
  if( !hasref && docmdist ) {
    error("with CMDIST one must use REFERENCE to setup the reference contact map");
  }

  // Read in switching functions
  std::vector<std::string> swdef(ncontacts);
  unsigned nswitch=0;
  for(unsigned i=0; i<ncontacts; ++i) {
    if( !parseNumbered( "SWITCH", i+1, swdef[i] ) ) {
      break;
    }
    nswitch++;
  }
  if( nswitch==0 ) {
    std::string sw;
    parse("SWITCH",sw);
    if(sw.length()==0) {
      error("no switching function specified use SWITCH keyword");
    }
    swdef.assign(ncontacts,sw);
  } else if( nswitch!=ncontacts ) {
    std::string num;
    Tools::convert(nswitch+1, num);
    error("missing SWITCH" + num + " keyword");
  }
  unsigned nnative=0;
  sfs.resize(ncontacts);
  for(unsigned i=0; i<ncontacts; ++i) {
    if( !refdist.empty() && addNativeDistance(swdef[i],refdist[i]) ) {
      nnative++;
    }
    std::string errors;
    sfs[i].set(swdef[i],errors);
    if( errors.length()!=0 ) {
      std::string num;
      if( nswitch>0 ) {
        Tools::convert(i+1, num);
      }
      error("problem reading SWITCH" + num + " keyword : " + errors );
    }
  }

  // Reference values from the reference structure
  if( !refdist.empty() ) {
    if( nnative>0 ) {
      std::string sn;
      Tools::convert(nnative,sn);
      log.printf("  REF of %s Q switching functions taken from the reference structure\n", sn.c_str());
    }
    if( docmdist ) {
      for(unsigned i=0; i<ncontacts; ++i) {
        double dfunc;
        reference[i]=sfs[i].calculate(refdist[i],dfunc);
      }
    } else if( nnative==0 && !usecutoff ) {
      error("without CMDIST the reference structure given to REFERENCE is only used to set REF in Q switching functions, and there are none without REF");
    } else {
      log.printf("  no reference value is subtracted from the contacts (the reference structure is used for the reference contact map only with CMDIST)\n");
    }
  }

  // Read in weight values
  weight.assign(ncontacts, 1.0);
  const bool numberedweight=parseAllNumbered("WEIGHT",weight);
  std::string wstr;
  parse("WEIGHT",wstr);
  if( wstr.length()>0 ) {
    if( numberedweight ) {
      error("use either WEIGHT or WEIGHT1, WEIGHT2, ...");
    }
    double w0=1.0;
    if( !Tools::convertNoexcept(wstr,w0) ) {
      error("cannot understand the value " + wstr + " given to the WEIGHT keyword");
    }
    weight.assign(ncontacts, w0);
  }
  bool normalize=false;
  parseFlag("NORMALIZE_WEIGHTS",normalize);
  if(normalize) {
    double wsum=0.;
    for(const auto w : weight) {
      wsum+=w;
    }
    if( wsum==0. ) {
      error("NORMALIZE_WEIGHTS cannot be used when the weights sum to zero");
    }
    for(auto & w : weight) {
      w/=wsum;
    }
    log.printf("  weights normalized to sum to one\n");
  }

  // Output details of all contacts
  for(unsigned i=0; i<ncontacts; ++i) {
    log.printf("  contact %u between atoms %u and %u: inflection point of the switching function at %s. Reference value %f, weight %f\n",
               i+1, ga_lista[i].serial(), gb_lista[i].serial(), ( sfs[i].description() ).c_str(), reference[i], weight[i] );
  }

  if(dosum) {
    addValueWithDerivatives();
    setNotPeriodic();
    log.printf("  colvar is sum of all contacts in contact map\n");
  }
  if(docmdist) {
    addValueWithDerivatives();
    setNotPeriodic();
    log.printf("  colvar is distance between the contact map matrix and the provided reference matrix\n");
  }

  if(dosum || docmdist) {
    docomp=false;
    deriv.resize(ncontacts);
  } else {
    // with components each contact has its own value: the calculation is done in serial
    serial=true;
    docomp=true;
  }

  // contact i is between atoms i and ncontacts+i
  std::vector<AtomNumber> atoms(ga_lista);
  atoms.insert(atoms.end(),gb_lista.begin(),gb_lista.end());
  requestAtoms(atoms);
  checkRead();
}

void ContactMap::calculate() {
  const unsigned ncontacts=sfs.size();

  unsigned stride=1;
  unsigned rank=0;
  if(!serial) {
    stride=comm.Get_size();
    rank=comm.Get_rank();
  }
  // OpenMP threads only when each of them has enough contacts
  unsigned nt=OpenMP::getNumThreads();
  if(nt*stride*(docomp ? minComponentsPerThread : minContactsPerThread)>ncontacts) {
    nt=1;
  }

  if(docomp) {
    // each contact writes only its own component
    #pragma omp parallel for num_threads(nt)
    for(unsigned i=0; i<ncontacts; ++i) {
      const Vector distance=(pbc ? pbcDistance(getPosition(i),getPosition(ncontacts+i))
                             : delta(getPosition(i),getPosition(ncontacts+i)));
      double dfunc=0.;
      const double coord = weight[i]*(sfs[i].calculateSqr(distance.modulo2(), dfunc) - reference[i]);
      const Vector der = weight[i]*dfunc*distance;
      Value* val=getPntrToComponent( i );
      setAtomsDerivatives( val, i, -der );
      setAtomsDerivatives( val, ncontacts+i, der );
      setBoxDerivatives( val, -Tensor(distance,der) );
      val->set(coord);
    }
    return;
  }

  double cv=0.;
  Tensor virial;
  if(!serial) {
    // each rank sets only its own contacts
    std::fill(deriv.begin(),deriv.end(),Vector(0.,0.,0.));
  }
  #pragma omp parallel num_threads(nt)
  {
    Tensor omp_virial;
    #pragma omp for reduction(+:cv) nowait
    for(unsigned i=rank; i<ncontacts; i+=stride) {
      const Vector distance=(pbc ? pbcDistance(getPosition(i),getPosition(ncontacts+i))
                             : delta(getPosition(i),getPosition(ncontacts+i)));
      double dfunc=0.;
      const double coord = weight[i]*(sfs[i].calculateSqr(distance.modulo2(), dfunc) - reference[i]);
      Vector der = weight[i]*dfunc*distance;
      if(docmdist) {
        der *= 2.*coord;
        cv  += coord*coord;
      } else {
        cv  += coord;
      }
      deriv[i] = der;
      omp_virial -= Tensor(distance,der);
    }
    #pragma omp critical
    virial += omp_virial;
  }

  if(!serial) {
    comm.Sum(cv);
    comm.Sum(&deriv[0][0],3*deriv.size());
    comm.Sum(virial);
  }

  #pragma omp parallel for num_threads(nt)
  for(unsigned i=0; i<ncontacts; ++i) {
    setAtomsDerivatives(i,-deriv[i]);
    setAtomsDerivatives(ncontacts+i,deriv[i]);
  }
  setValue           (cv);
  setBoxDerivatives  (virial);
}

}
}
