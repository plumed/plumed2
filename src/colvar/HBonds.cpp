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
#include "Colvar.h"
#include "core/ActionRegister.h"
#include "tools/NeighborList.h"
#include "tools/SwitchingFunction.h"
#include "tools/Communicator.h"
#include "tools/OpenMP.h"
#include "tools/Tools.h"

#include <algorithm>
#include <cmath>
#include <iterator>
#include <limits>

namespace PLMD {
namespace colvar {

//+PLUMEDOC COLVAR HBONDS
/*
Calculate the number of hydrogen bonds using the Wernet-Nilsson criterion.

This action counts the hydrogen bonds between a set of donors (with their hydrogens) and a set of acceptors.
A donor D, one of its hydrogens H and an acceptor A form a hydrogen bond according to the criterion of Wernet and Nilsson
(the one used by the `wernet_nilsson` function of [mdtraj](https://mdtraj.org)) if

$$
r_{DA} < 0.33\,\textrm{nm} - 0.000044\,\textrm{nm}\,\delta_{HDA}^2
$$

where $r_{DA}$ is the donor-acceptor distance and $\delta_{HDA}$ is the angle, in degrees, between the D-H and the D-A vectors.
The criterion is a condition on a single variable,

$$
s_{DHA} = r_{DA} + k\, \delta_{HDA}^2 \qquad k = 0.000044\,\textrm{nm}/\textrm{deg}^2 \simeq 0.1444\,\textrm{nm}/\textrm{rad}^2
$$

which has to be smaller than 0.33 nm. The step function is replaced by a switching function of $s_{DHA}$, so that the CV is

$$
N_{HB} = \sum_{D} \sum_{H \in D} \sum_{A \neq D} \sigma(s_{DHA})
$$

and it has continuous derivatives. The default switching function is `{CUBIC D_0=0.31 D_MAX=0.35}` (in nm),
which is one below 0.31 nm, zero above 0.35 nm and 1/2 at 0.33 nm. You can change it with the SWITCH keyword
(see [LESS_THAN](LESS_THAN.md) for the available switching functions). For instance a very sharp switching function
such as `{CUBIC D_0=0.3299 D_MAX=0.3301}` reproduces the number of hydrogen bonds computed by mdtraj, but should not be used to bias a simulation.
The angle cutoff of 45 degrees that mdtraj also applies is not needed here: at that angle the criterion
already requires a donor-acceptor distance shorter than 0.24 nm.

The donors, the hydrogens and the acceptors are given as three lists of atoms, in any order.
Each hydrogen is assigned to the closest donor, which must be closer than 0.13 nm (a covalent D-H bond).
This assignment is done once, on the first configuration.
Hydrogens that are not bonded to any of the donors are ignored, as are the donors without hydrogens,
and the number of each is written in the log. A donor that is also an acceptor does not accept from itself.
An atom cannot be both in DONORS and in HYDROGENS, and an atom cannot be repeated in the same list.

The [MOLINFO](MOLINFO.md) shortcuts `@hbdonors`, `@polarhydrogens` and `@hbacceptors` select the N and O atoms bonded to a hydrogen,
the hydrogens bonded to an N or an O atom and all the N and O atoms of the reference structure respectively.
As in mdtraj, all N and O atoms are acceptors.
In the example below all the hydrogen bonds formed within a protein are counted:

```plumed
#SETTINGS MOLFILE=regtest/basic/rt-hbonds/ref.pdb
MOLINFO STRUCTURE=regtest/basic/rt-hbonds/ref.pdb
hb: HBONDS DONORS=@hbdonors HYDROGENS=@polarhydrogens ACCEPTORS=@hbacceptors
PRINT ARG=hb FILE=colvar
```

The atoms can also be given explicitly. This input counts the hydrogen bonds in a box of 64 water molecules,
in which the oxygens are atoms 1, 4, 7... and the hydrogens are the atoms that follow each oxygen.
The oxygens are both donors and acceptors, and a neighbor list on the donor-acceptor distance is used
(see [COORDINATION](COORDINATION.md) for the choice between NLIST and NLISTCELLS, and of NL_CUTOFF and NL_STRIDE; NL_CUTOFF should be larger than D_MAX of the switching function):

```plumed
hb: HBONDS ...
  DONORS=1-192:3 HYDROGENS=2-192:3,3-192:3 ACCEPTORS=1-192:3
  SWITCH={CUBIC D_0=0.31 D_MAX=0.35}
  NLISTCELLS NL_CUTOFF=0.6 NL_STRIDE=10
...
PRINT ARG=hb FILE=colvar
```

Only N and O atoms (or other electronegative atoms) should be given as donors and acceptors:
the hydrogens bonded to any atom in DONORS are used, so if a carbon was in the list its hydrogens would be donated.

*/
//+ENDPLUMEDOC

class HBonds : public Colvar {
  bool pbc=true;
  bool serial=false;
  unsigned nD=0;
  unsigned nA=0;
  unsigned nH=0;
  SwitchingFunction switchingFunction;
/// coefficient of the squared angle, in length units per rad^2
  double kangle=0.0;
/// largest value of s for which the switching function is not zero
  double dmax=0.0;
  double dmax2=0.0;
/// largest D-H bond length
  double bondCutoff=0.0;
/// neighbor list of the donor-acceptor pairs
  std::unique_ptr<NeighborList> nl;
  unsigned nlStride=0;
  bool firsttime=true;
  bool invalidateList=true;
/// hydrogens of donor i are hydrogens[hfirst[i]..hfirst[i+1]), as indexes in the positions array
  std::vector<unsigned> hfirst;
  std::vector<unsigned> hydrogens;
  bool paired=false;
/// scratch positions for the neighbor list
  std::vector<Vector> nlpos;
/// the D-H vectors and their inverse lengths (zero for a degenerate bond), in the order of hydrogens
  std::vector<Vector> dh;
  std::vector<double> invdh;
/// Assign each hydrogen to the closest donor, on the first configuration
  void pairHydrogens();
/// Errors for empty lists, atoms repeated in a list and atoms both in DONORS and HYDROGENS
  void checkAtomLists(const std::vector<AtomNumber>& donors,
                      const std::vector<AtomNumber>& hlist,
                      const std::vector<AtomNumber>& acceptors) const;
/// Read SWITCH and set the constants of the criterion in the length units
  void parseSwitchingFunction();
/// Read the neighbor-list keywords and create the list of donor-acceptor pairs
  void setupNeighborList(const std::vector<AtomNumber>& donors,const std::vector<AtomNumber>& acceptors);
  Vector distance(const Vector&a,const Vector&b) const {
    return pbc ? pbcDistance(a,b) : delta(a,b);
  }
public:
  explicit HBonds(const ActionOptions&);
  ~HBonds();
  static void registerKeywords( Keywords& keys );
  void prepare() override;
  void calculate() override;
};

PLUMED_REGISTER_ACTION(HBonds,"HBONDS")

void HBonds::registerKeywords( Keywords& keys ) {
  Colvar::registerKeywords(keys);
  keys.add("atoms","DONORS","The atoms that can donate a hydrogen bond (usually the N and O atoms bonded to a hydrogen, e.g. the @hbdonors MOLINFO shortcut)");
  keys.add("atoms","HYDROGENS","The hydrogens that can be donated (e.g. the @polarhydrogens MOLINFO shortcut). Each is assigned to the closest donor, which must be closer than 0.13 nm");
  keys.add("atoms","ACCEPTORS","The atoms that can accept a hydrogen bond (usually the N and O atoms, e.g. the @hbacceptors MOLINFO shortcut)");
  keys.add("optional","SWITCH","The switching function applied to s = r_DA + k delta^2. The default is {CUBIC D_0=0.31 D_MAX=0.35} (in nm)");
  keys.linkActionInDocs("SWITCH","LESS_THAN");
  keys.addFlag("SERIAL",false,"Perform the calculation in serial - for debug purpose");
  keys.addFlag("NLIST",false,"Use a neighbor list of the donor-acceptor pairs to speed up the calculation");
  keys.addFlag("NLISTCELLS",false,"Use a neighbor list of the donor-acceptor pairs to speed up the calculation - use the cell list implementation instead of the classical one");
  keys.add("optional","NL_CUTOFF","The cutoff on the donor-acceptor distance for the neighbor list (it should be larger than D_MAX of the switching function)");
  keys.add("optional","NL_STRIDE","The frequency with which the neighbor list is updated");
  keys.setValueDescription("scalar","the number of hydrogen bonds");
}

HBonds::HBonds(const ActionOptions&ao):
  PLUMED_COLVAR_INIT(ao) {
  parseFlag("SERIAL",serial);
  bool nopbc=!pbc;
  parseFlag("NOPBC",nopbc);
  pbc=!nopbc;

  std::vector<AtomNumber> donors,acceptors,hlist;
  parseAtomList("DONORS",donors);
  parseAtomList("HYDROGENS",hlist);
  parseAtomList("ACCEPTORS",acceptors);
  checkAtomLists(donors,hlist,acceptors);
  parseSwitchingFunction();
  setupNeighborList(donors,acceptors);
  checkRead();

  nD=donors.size();
  nA=acceptors.size();
  nH=hlist.size();
  // positions: donors, acceptors, hydrogens
  std::vector<AtomNumber> all(donors);
  all.insert(all.end(),acceptors.begin(),acceptors.end());
  all.insert(all.end(),hlist.begin(),hlist.end());
  addValueWithDerivatives();
  setNotPeriodic();
  requestAtoms(all);

  log.printf("  %u donors, %u hydrogens and %u acceptors\n",nD,nH,nA);
  log.printf("  switching function on s = r_DA + %f delta_HDA^2 (delta in radians): %s\n",kangle,switchingFunction.description().c_str());
  if(pbc) {
    log.printf("  using periodic boundary conditions\n");
  } else {
    log.printf("  without periodic boundary conditions\n");
  }
  log<<"  Bibliography "<<cite("Wernet et al, Science 304, 995 (2004)")<<"\n";
}

void HBonds::checkAtomLists(const std::vector<AtomNumber>& donors,
                            const std::vector<AtomNumber>& hlist,
                            const std::vector<AtomNumber>& acceptors) const {
  if(donors.empty()) {
    error("DONORS should contain at least one atom");
  }
  if(hlist.empty()) {
    error("HYDROGENS should contain at least one atom");
  }
  if(acceptors.empty()) {
    error("ACCEPTORS should contain at least one atom");
  }
  auto sorted=[](std::vector<AtomNumber> v) {
    std::sort(v.begin(),v.end());
    return v;
  };
  const auto sd=sorted(donors);
  const auto sh=sorted(hlist);
  const auto sa=sorted(acceptors);
  // an atom repeated in a list would be counted twice
  for(const auto & [name,v] : {
        std::pair{"DONORS",&sd},std::pair{"HYDROGENS",&sh},std::pair{"ACCEPTORS",&sa}
      }) {
    if(std::adjacent_find(v->begin(),v->end())!=v->end()) {
      error(std::string(name)+" contains the same atom more than once");
    }
  }
  // a donor in HYDROGENS would be assigned to another donor (e.g. an O to the H of another O-H)
  std::vector<AtomNumber> common;
  std::set_intersection(sd.begin(),sd.end(),sh.begin(),sh.end(),std::back_inserter(common));
  if(!common.empty()) {
    error("atom "+std::to_string(common[0].serial())+" is both in DONORS and in HYDROGENS");
  }
}

void HBonds::parseSwitchingFunction() {
  // the default and the constants of the criterion are in nm, converted to the length units
  const double nm=1.0/getUnits().getLength();
  std::string sw;
  parse("SWITCH",sw);
  if(sw.empty()) {
    sw="CUBIC D_0="+std::to_string(0.31*nm)+" D_MAX="+std::to_string(0.35*nm);
  }
  std::string errors;
  switchingFunction.set(sw,errors);
  if(errors.length()>0) {
    error("problem reading SWITCH keyword : " + errors );
  }
  dmax=switchingFunction.get_dmax();
  dmax2=switchingFunction.get_dmax2();
  // 0.00044 Angstrom/deg^2, converted to nm/rad^2
  const double rad=180.0/pi;
  kangle=0.000044*rad*rad*nm;
  bondCutoff=0.13*nm;
}

void HBonds::setupNeighborList(const std::vector<AtomNumber>& donors,const std::vector<AtomNumber>& acceptors) {
  bool doneighClassic=false;
  bool doneighCells=false;
  parseFlag("NLIST",doneighClassic);
  parseFlag("NLISTCELLS",doneighCells);
  if(doneighClassic && doneighCells) {
    error("use only one of NLIST and NLISTCELLS");
  }
  if(!doneighClassic && !doneighCells) {
    // all the donor-acceptor pairs
    nl=Tools::make_unique<NeighborList>(donors,acceptors,serial,false,pbc,getPbc(),comm);
    return;
  }
  double nlCutoff=0.0;
  parse("NL_CUTOFF",nlCutoff);
  if(nlCutoff<=0.0) {
    error("NL_CUTOFF should be explicitly specified and positive");
  }
  int st=0;
  parse("NL_STRIDE",st);
  if(st<=0) {
    error("NL_STRIDE should be explicitly specified and positive");
  }
  nlStride=st;
  // s >= r_DA: the pairs beyond D_MAX do not contribute
  if(nlCutoff<dmax) {
    log<<"  WARNING: NL_CUTOFF is smaller than D_MAX of the switching function: some hydrogen bonds may be missed\n";
  }
  nl=Tools::make_unique<NeighborList>(donors,acceptors,serial,false,pbc,getPbc(),comm,nlCutoff,nlStride,doneighCells);
  log.printf("  using neighbor lists of the donor-acceptor pairs with cutoff %f, updated every %u steps\n",nlCutoff,nlStride);
}

HBonds::~HBonds()=default;

void HBonds::pairHydrogens() {
  // each hydrogen goes with the closest donor within bondCutoff (the search uses a cell list)
  std::vector<AtomNumber> hatoms(nH),datoms(nD);
  std::vector<Vector> pos(nH+nD);
  for(unsigned i=0; i<nH; ++i) {
    hatoms[i]=getAbsoluteIndex(nD+nA+i);
    pos[i]=getPosition(nD+nA+i);
  }
  for(unsigned i=0; i<nD; ++i) {
    datoms[i]=getAbsoluteIndex(i);
    pos[nH+i]=getPosition(i);
  }
  NeighborList bonds(hatoms,datoms,true,false,pbc,getPbc(),comm,bondCutoff,1,true);
  bonds.update(pos);
  constexpr unsigned none=std::numeric_limits<unsigned>::max();
  std::vector<unsigned> donorOf(nH,none);
  std::vector<double> closest(nH,bondCutoff*bondCutoff);
  for(unsigned p=0; p<bonds.size(); ++p) {
    const auto [h,d]=bonds.getClosePair(p);
    if(hatoms[h]==datoms[d-nH]) {
      continue;
    }
    const double d2=distance(pos[h],pos[d]).modulo2();
    // the list is built in serial, so the order of the pairs (and the choice in a tie) is reproducible
    if(d2<closest[h]) {
      closest[h]=d2;
      donorOf[h]=d-nH;
    }
  }
  hfirst.assign(nD+1,0);
  unsigned unpaired=0;
  for(unsigned h=0; h<nH; ++h) {
    if(donorOf[h]==none) {
      ++unpaired;
    } else {
      ++hfirst[donorOf[h]+1];
    }
  }
  for(unsigned i=0; i<nD; ++i) {
    hfirst[i+1]+=hfirst[i];
  }
  hydrogens.assign(hfirst[nD],0);
  std::vector<unsigned> fill(hfirst.begin(),hfirst.end()-1);
  for(unsigned h=0; h<nH; ++h) {
    if(donorOf[h]!=none) {
      hydrogens[fill[donorOf[h]]++]=nD+nA+h;
    }
  }
  unsigned nobonds=0;
  for(unsigned i=0; i<nD; ++i) {
    if(hfirst[i+1]==hfirst[i]) {
      ++nobonds;
    }
  }
  if(hydrogens.empty()) {
    error("none of the HYDROGENS is closer than 0.13 nm to one of the DONORS");
  }
  log.printf("  HBONDS %s: %u D-H bonds found; %u hydrogens are not bonded to a donor and %u donors have no hydrogen: they are ignored\n",
             getLabel().c_str(),static_cast<unsigned>(hydrogens.size()),unpaired,nobonds);
  paired=true;
}

void HBonds::prepare() {
  if(nlStride>0) {
    // the list is also updated on exchange steps, where the coordinates are swapped
    invalidateList=firsttime || getStep()%nlStride==0 || getExchangeStep();
    firsttime=false;
  }
}

void HBonds::calculate() {
  if(!paired) {
    pairHydrogens();
  }
  if(nlStride>0 && invalidateList) {
    nlpos.assign(getPositions().begin(),getPositions().begin()+nD+nA);
    nl->update(nlpos);
  }

  const unsigned stride=serial ? 1 : comm.Get_size();
  const unsigned rank=serial ? 0 : comm.Get_rank();
  const unsigned nn=nl->size();
  unsigned nt=OpenMP::getNumThreads();
  if(nt*stride*10>nn) {
    nt=1;
  }
  const unsigned perRank=(nn+stride-1)/stride;
  const unsigned start=std::min(rank*perRank,nn);
  const unsigned end=std::min(start+perRank,nn);

  // each D-H vector is used with all the acceptors close to the donor: computed once
  dh.resize(hydrogens.size());
  invdh.resize(hydrogens.size());
  for(unsigned d=0; d<nD; ++d) {
    for(unsigned k=hfirst[d]; k<hfirst[d+1]; ++k) {
      dh[k]=distance(getPosition(d),getPosition(hydrogens[k]));
      const double u2=dh[k].modulo2();
      invdh[k]=u2<epsilon ? 0.0 : 1.0/std::sqrt(u2);
    }
  }

  const unsigned natoms=getNumberOfAtoms();
  std::vector<Vector> deriv(natoms);
  Tensor virial;
  double nhb=0.0;

  #pragma omp parallel num_threads(nt)
  {
    std::vector<Vector> omp_deriv(nt>1 ? natoms : 0);
    Tensor omp_virial;
    std::vector<Vector> & myderiv(nt>1 ? omp_deriv : deriv);
    Tensor & myvirial(nt>1 ? omp_virial : virial);
    #pragma omp for reduction(+:nhb) nowait
    for(unsigned p=start; p<end; ++p) {
      const auto [d,a]=nl->getClosePair(p);
      if(hfirst[d]==hfirst[d+1] || getAbsoluteIndex(d)==getAbsoluteIndex(a)) {
        continue;
      }
      const Vector v=distance(getPosition(d),getPosition(a));
      const double r2=v.modulo2();
      // s >= r_DA, so the pair does not contribute beyond dmax
      if(r2>=dmax2 || r2<epsilon) {
        continue;
      }
      const double r=std::sqrt(r2);
      const double invr=1.0/r;
      for(unsigned k=hfirst[d]; k<hfirst[d+1]; ++k) {
        const double invu=invdh[k];
        if(invu==0.0) {
          continue;
        }
        const Vector & u=dh[k];
        const double c=std::clamp(dotProduct(u,v)*invu*invr,-1.0,1.0);
        // delta^2 >= 2(1-cos(delta)): most hydrogens point away from the acceptor and are discarded without acos
        if(r+2.0*kangle*(1.0-c)>=dmax) {
          continue;
        }
        const unsigned h=hydrogens[k];
        const double delta=std::acos(c);
        const double s=r+kangle*delta*delta;
        if(s>=dmax) {
          continue;
        }
        double df=0.0;
        nhb+=switchingFunction.calculate(s,df);
        // df is the derivative divided by s
        const double g=df*s;
        if(g==0.0) {
          continue;
        }
        // d(delta^2)/dc = -2 delta/sin(delta), which is finite for delta -> 0
        const double sind=std::sqrt(std::max(0.0,1.0-c*c));
        double f;
        if(delta<1e-4) {
          f=1.0+delta*delta/6.0;
        } else if(sind<1e-12) {
          f=0.0;
        } else {
          f=delta/sind;
        }
        const double gc=-2.0*kangle*f*g;
        const Vector dcdu=invu*invr*v-c*invu*invu*u;
        const Vector dcdv=invu*invr*u-c*invr*invr*v;
        const Vector gu=gc*dcdu;
        const Vector gv=g*invr*v+gc*dcdv;
        myderiv[h]+=gu;
        myderiv[a]+=gv;
        myderiv[d]-=gu+gv;
        myvirial-=Tensor(gu,u)+Tensor(gv,v);
      }
    }
    if(nt>1) {
      #pragma omp critical
      {
        for(unsigned i=0; i<natoms; ++i) {
          deriv[i]+=omp_deriv[i];
        }
        virial+=omp_virial;
      }
    }
  }

  if(!serial) {
    comm.Sum(nhb);
    comm.Sum(&deriv[0][0],3*deriv.size());
    comm.Sum(virial);
  }
  for(unsigned i=0; i<natoms; ++i) {
    setAtomsDerivatives(i,deriv[i]);
  }
  setValue(nhb);
  setBoxDerivatives(virial);
}

}
}
