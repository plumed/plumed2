/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2017-2021 Omar Valsson

   This file is part of S2 contact model module

   The S2 contact model module is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   The S2 contact model module is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with the S2 contact model module.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

#include "colvar/Colvar.h"
#include "tools/NeighborList.h"
#include "tools/Communicator.h"
#include "core/ActionRegister.h"
#include "core/PlumedMain.h"


#include <string>
#include <cmath>

namespace PLMD {
namespace s2cm {


//

//+PLUMEDOC S2CMMOD_COLVAR S2CM
/*
S2 contact model CV.

This CV was used in \cite Palazzesi_s2_2017, based on NH order parameter from \cite Zhang_s2_2002 and methyl order parameter from \cite Ming_s2_2004. Input parameters can be found in the relevant papers.

\par Examples


*/
//+ENDPLUMEDOC







class S2ContactModel : public Colvar {
  bool pbc_;
  bool serial_;
  std::unique_ptr<NeighborList> nl;
  bool invalidateList;
  bool firsttime;
  //
  double r_eff_;
  double inv_r_eff_;
  double prefactor_a_;
  double exp_b_;
  double offset_c_;
  double n_i_;
  double total_prefactor_;
  double r_globalshift_;

  enum ModelType {methyl,nh} modeltype_;

  //
public:
  explicit S2ContactModel(const ActionOptions&);
  ~S2ContactModel();
  virtual void calculate();
  virtual void prepare();
  static void registerKeywords( Keywords& keys );
};

PLUMED_REGISTER_ACTION(S2ContactModel,"S2CM")

void S2ContactModel::registerKeywords( Keywords& keys ) {
  Colvar::registerKeywords(keys);
  keys.addFlag("SERIAL",false,"Perform the calculation in serial - for debug purpose");
  keys.addFlag("NLIST",false,"Use a neighbour list to speed up the calculation");
  keys.add("optional","NL_CUTOFF","The cutoff for the neighbour list");
  keys.add("optional","NL_STRIDE","The frequency with which we are updating the atoms in the neighbour list");
  keys.add("atoms","METHYL_ATOM","the methyl carbon atom of the residue (i)");
  keys.add("atoms","NH_ATOMS","the hydrogen atom of the NH group of the residue (i) and carbonyl oxygen of the preceding residue (i-1)");
  keys.add("atoms","HEAVY_ATOMS","the heavy atoms to be included in the calculation");
  //
  keys.add("compulsory","R_EFF","the effective distance, r_eff in the equation, given in nm.");
  keys.add("compulsory","PREFACTOR_A","the prefactor, a in the equation");
  keys.add("compulsory","EXPONENT_B","the exponent, b in the equation");
  keys.add("compulsory","OFFSET_C","the offset, c in the equation");
  keys.add("compulsory","N_I"," n_i in the equation");
  keys.add("optional","R_SHIFT","shift all distances by given amount");
}

S2ContactModel::S2ContactModel(const ActionOptions&ao):
  PLUMED_COLVAR_INIT(ao),
  pbc_(true),
  serial_(false),
  invalidateList(true),
  firsttime(true),
  r_eff_(0.0),
  inv_r_eff_(0.0),
  prefactor_a_(0.0),
  exp_b_(0.0),
  offset_c_(0.0),
  n_i_(0.0),
  total_prefactor_(0.0),
  r_globalshift_(0.0),
  modeltype_(methyl)
{

  parseFlag("SERIAL",serial_);

  std::vector<AtomNumber> methyl_atom;
  parseAtomList("METHYL_ATOM",methyl_atom);
  std::vector<AtomNumber> nh_atoms;
  parseAtomList("NH_ATOMS",nh_atoms);

  if(methyl_atom.size()==0 && nh_atoms.size()==0) {
    error("you have to give either METHYL_ATOM or NH_ATOMS");
  }
  if(methyl_atom.size()>0 && nh_atoms.size()>0) {
    error("you cannot give both METHYL_ATOM or NH_ATOMS");
  }
  if(methyl_atom.size()>0 && methyl_atom.size()!=1) {
    error("you should give one atom in METHYL_ATOM, the methyl carbon atom of the residue");
  }
  if(nh_atoms.size()>0 && nh_atoms.size()!=2) {
    error("you should give two atoms in NH_ATOMS, the hydrogen atom of the NH group of the residue (i) and carbonyl oxygen of the preceding residue (i-1)");
  }

  std::vector<AtomNumber> heavy_atoms;
  parseAtomList("HEAVY_ATOMS",heavy_atoms);
  if(heavy_atoms.size()==0) {
    error("HEAVY_ATOMS is not given");
  }

  std::vector<AtomNumber> main_atoms;
  if(methyl_atom.size()>0) {
    modeltype_= methyl;
    main_atoms = methyl_atom;
  } else if(nh_atoms.size()>0) {
    modeltype_= nh;
    main_atoms = nh_atoms;
  }

  bool nopbc=!pbc_;
  parseFlag("NOPBC",nopbc);
  pbc_=!nopbc;

  // neighbor list stuff
  bool doneigh=false;
  double nl_cut=0.0;
  int nl_st=0;
  parseFlag("NLIST",doneigh);
  if(doneigh) {
    parse("NL_CUTOFF",nl_cut);
    if(nl_cut<=0.0) error("NL_CUTOFF should be explicitly specified and positive");
    parse("NL_STRIDE",nl_st);
    if(nl_st<=0) error("NL_STRIDE should be explicitly specified and positive");
  }

  parse("R_EFF",r_eff_);
  inv_r_eff_ = 1.0/r_eff_;
  parse("PREFACTOR_A",prefactor_a_);
  parse("EXPONENT_B",exp_b_);
  parse("OFFSET_C",offset_c_);
  unsigned int n_i_int;
  parse("N_I",n_i_int);
  n_i_ = static_cast<double>(n_i_int);
  total_prefactor_ = prefactor_a_/pow(n_i_,exp_b_);
  //
  parse("R_SHIFT",r_globalshift_);

  checkRead();

  addValueWithDerivatives();
  setNotPeriodic();

  bool dopair=false;
  if(doneigh) {
    nl=Tools::make_unique<NeighborList>(main_atoms,heavy_atoms,serial_,dopair,pbc_,getPbc(),comm,nl_cut,nl_st);
  }
  else {
    nl=Tools::make_unique<NeighborList>(main_atoms,heavy_atoms,serial_,dopair,pbc_,getPbc(),comm);
  }

  requestAtoms(nl->getFullAtomList());


  log.printf("  NMR S2 contact model CV, please read and cite ");
  log << plumed.cite("Palazzesi, Valsson, and Parrinello, J. Phys. Chem. Lett. 8, 4752 (2017) - DOI:10.1021/acs.jpclett.7b01770");

  if(modeltype_==methyl) {
    log << plumed.cite("Ming and Bruschweiler, J. Biomol. NMR, 29, 363 (2004) - DOI:10.1023/B:JNMR.0000032612.70767.35");
    log.printf("\n");
    log.printf("  calculation of methyl order parameter using atom %d \n",methyl_atom[0].serial());
  }
  else if(modeltype_==nh) {
    log << plumed.cite("Zhang and Bruschweiler, J. Am. Chem. Soc. 124, 12654 (2002) - DOI:10.1021/ja027847a");
    log.printf("\n");
    log.printf("  calculation of NH order parameter using atoms %d and %d\n",nh_atoms[0].serial(),nh_atoms[1].serial());
  }
  log.printf("  heavy atoms used in the calculation (%u in total):\n",static_cast<unsigned int>(heavy_atoms.size()));
  for(unsigned int i=0; i<heavy_atoms.size(); ++i) {
    if( (i+1) % 25 == 0 ) {log.printf("  \n");}
    log.printf("  %d", heavy_atoms[i].serial());
  }
  log.printf("  \n");
  log.printf("  total number of distances: %u\n",nl->size());
  //
  log.printf("  using parameters");
  log.printf(" r_eff=%f,",r_eff_);
  log.printf(" a=%f,",prefactor_a_);
  log.printf(" b=%f,",exp_b_);
  log.printf(" c=%f,",offset_c_);
  log.printf(" n_i=%u",n_i_int);
  if(r_globalshift_!=0.0) {
    log.printf(", r_shift=%f",r_globalshift_);
  }
  log.printf("\n");
  if(pbc_) {
    log.printf("  using periodic boundary conditions\n");
  } else {
    log.printf("  without periodic boundary conditions\n");
  }
  if(doneigh) {
    log.printf("  using neighbor lists with\n");
    log.printf("  update every %d steps and cutoff %f\n",nl_st,nl_cut);
  }
  if(serial_) {
    log.printf("  calculation done in serial\n");
  }
}

S2ContactModel::~S2ContactModel() {
}

void S2ContactModel::prepare() {
  if(nl->getStride()>0) {
    if(firsttime || (getStep()%nl->getStride()==0)) {
      requestAtoms(nl->getFullAtomList());
      invalidateList=true;
      firsttime=false;
    } else {
      requestAtoms(nl->getReducedAtomList());
      invalidateList=false;
      if(getExchangeStep()) error("Neighbor lists should be updated on exchange steps - choose a NL_STRIDE which divides the exchange stride!");
    }
    if(getExchangeStep()) firsttime=true;
  }
}

// calculator
void S2ContactModel::calculate() {

  Tensor virial;
  std::vector<Vector> deriv(getNumberOfAtoms());

  if(nl->getStride()>0 && invalidateList) {
    nl->update(getPositions());
  }

  unsigned int stride=comm.Get_size();
  unsigned int rank=comm.Get_rank();
  if(serial_) {
    stride=1;
    rank=0;
  }

  double contact_sum = 0.0;

  const unsigned int nn=nl->size();

  for(unsigned int i=rank; i<nn; i+=stride) {
    Vector distance;
    unsigned int i0=nl->getClosePair(i).first;
    unsigned int i1=nl->getClosePair(i).second;
    if(getAbsoluteIndex(i0)==getAbsoluteIndex(i1)) {continue;}

    if(pbc_) {
      distance=pbcDistance(getPosition(i0),getPosition(i1));
    } else {
      distance=delta(getPosition(i0),getPosition(i1));
    }

    double exp_arg = exp(-(distance.modulo()-r_globalshift_)*inv_r_eff_);
    contact_sum += exp_arg;

    exp_arg /= distance.modulo();
    Vector dd(exp_arg*distance);
    Tensor vv(dd,distance);
    deriv[i0]-=dd;
    deriv[i1]+=dd;
    virial-=vv;
  }

  if(!serial_) {
    comm.Sum(contact_sum);
    if(!deriv.empty()) {
      comm.Sum(&deriv[0][0],3*deriv.size());
    }
    comm.Sum(virial);
  }

  double value = tanh(total_prefactor_*contact_sum);
  // using that d/dx[tanh(x)]= 1-[tanh(x)]^2
  double deriv_f = -inv_r_eff_*total_prefactor_*(1.0-value*value);
  value -= offset_c_;

  for(unsigned int i=0; i<deriv.size(); ++i) {
    setAtomsDerivatives(i,deriv_f*deriv[i]);
  }
  setValue(value);
  setBoxDerivatives(deriv_f*virial);

}
}
}
