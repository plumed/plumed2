/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2018 The plumed team
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
/*
 This class was originally written by Alexander Jussupow with some contribution
 by Carlo Camilloni
*/

#include "colvar/Colvar.h"
#include "colvar/ActionRegister.h"
#include "core/ActionSet.h"
#include "core/PlumedMain.h"
#include "core/SetupMolInfo.h"
#include "tools/Communicator.h"
#include "tools/OpenMP.h"
#include "tools/Pbc.h"

#include <string>
#include <cmath>
#include <map>

#ifdef __PLUMED_HAS_ARRAYFIRE
#include <arrayfire.h>
#include <af/util.h>
#endif

using namespace std;

namespace PLMD {
namespace isdb {

//+PLUMEDOC ISDB_COLVAR SAXSGPU
/*
Calculate SAXS scattered intensity on GPU  using the Debye equation.

Intensities are calculated for a set of scattering lenght set using QVALUES numbered keywords, QVALUE cannot be 0.
Structure factors can be either assigned using a polynomial expansion to any order using the PARAMETERS keywords;
automatically assigned to atoms using the ATOMISTIC flag reading a PDB file, a correction for the water density is automatically added;
automatically assigned to Martini pseudoatoms usign the MARTINI flag.
The calculated intensities can be scaled using the SCEXP keywords. This is applied by rescaling the structure factors.
Experimental reference intensities can be added using the ADDEXP and EXPINT flag and keywords.

\par Examples
in the following example the saxs intensities for a martini model are calculated. structure factors
are obtained from the pdb file indicated in the MOLINFO.

\plumedfile
MOLINFO STRUCTURE=template.pdb

SAXSGPU ...
LABEL=saxs
ATOMS=1-355
ADDEXP
SCEXP=3920000
MARTINI
QVALUE1=0.02 EXPINT1=1.0902
QVALUE2=0.05 EXPINT2=0.790632
QVALUE3=0.08 EXPINT3=0.453808
QVALUE4=0.11 EXPINT4=0.254737
QVALUE5=0.14 EXPINT5=0.154928
QVALUE6=0.17 EXPINT6=0.0921503
QVALUE7=0.2 EXPINT7=0.052633
QVALUE8=0.23 EXPINT8=0.0276557
QVALUE9=0.26 EXPINT9=0.0122775
QVALUE10=0.29 EXPINT10=0.00880634
QVALUE11=0.32 EXPINT11=0.0137301
QVALUE12=0.35 EXPINT12=0.0180036
QVALUE13=0.38 EXPINT13=0.0193374
QVALUE14=0.41 EXPINT14=0.0210131
QVALUE15=0.44 EXPINT15=0.0220506
...

PRINT ARG=(saxs\.q_.*),(saxs\.exp_.*) FILE=colvar STRIDE=1

\endplumedfile

*/
//+ENDPLUMEDOC

class SAXSGPU : public Colvar {
private:
  bool                pbc;
  bool                serial;
  int                 deviceid;
  vector<double>      q_list;
#ifdef  __PLUMED_HAS_ARRAYFIRE
  af::array           FF_value;
  af::array           q_value;
#endif
  void getMartiniSFparam(const vector<AtomNumber> &atoms, vector<vector<long double> > &parameter);
  void calculateASF(const vector<AtomNumber> &atoms, vector<vector<long double> > &FF_tmp, const double rho);

public:
  static void registerKeywords( Keywords& keys );
  explicit SAXSGPU(const ActionOptions&);
  virtual void calculate();
};

PLUMED_REGISTER_ACTION(SAXSGPU,"SAXSGPU")

void SAXSGPU::registerKeywords(Keywords& keys) {
  Colvar::registerKeywords( keys );
  componentsAreNotOptional(keys);
  useCustomisableComponents(keys);
  keys.addFlag("SERIAL",false,"Perform the calculation in serial - for debug purpose");
  keys.add("compulsory","DEVICEID","0","Identifier of the GPU to be used");
  keys.addFlag("ATOMISTIC",false,"calculate SAXS for an atomistic model");
  keys.addFlag("MARTINI",false,"calculate SAXS for a Martini model");
  keys.add("atoms","ATOMS","The atoms to be included in the calculation, e.g. the whole protein.");
  keys.add("numbered","QVALUE","Selected scattering lenghts in Angstrom are given as QVALUE1, QVALUE2, ... .");
  keys.add("numbered","PARAMETERS","Used parameter Keywords like PARAMETERS1, PARAMETERS2. These are used to calculate the structure factor for the i-th atom/bead.");
  keys.add("compulsory","WATERDENS","0.334","Density of the water to be used for the correction of atomistic structure factors.");
  keys.addFlag("ADDEXP",false,"Set to TRUE if you want to have fixed components with the experimental values.");
  keys.add("numbered","EXPINT","Add an experimental value for each q value.");
  keys.add("compulsory","SCEXP","1.0","SCALING value of the experimental data. Usefull to simplify the comparison.");
  keys.addOutputComponent("q","default","the # SAXS of q");
  keys.addOutputComponent("exp","ADDEXP","the # experimental intensity");
}

SAXSGPU::SAXSGPU(const ActionOptions&ao):
  PLUMED_COLVAR_INIT(ao),
  pbc(true),
  serial(false),
  deviceid(0)
{
#ifndef __PLUMED_HAS_ARRAYFIRE
  error("SAXSGPU can only be used if ARRAYFIRE is installed");
#else
  std::vector<AtomNumber> atoms;
  parseAtomList("ATOMS",atoms);
  const unsigned size = atoms.size();

  parseFlag("SERIAL",serial);

  bool nopbc=!pbc;
  parseFlag("NOPBC",nopbc);
  pbc=!nopbc;

  parse("DEVICEID",deviceid);
  af::setDevice(deviceid);
  af::info();

  long double scexp = 0;
  parse("SCEXP",scexp);
  if(scexp==0) scexp=1.0;

  unsigned ntarget=0;
  for(unsigned i=0;; ++i) {
    double t_list;
    if( !parseNumbered( "QVALUE", i+1, t_list) ) break;
    q_list.push_back(t_list);
    ntarget++;
  }
  const unsigned numq = ntarget;

  bool atomistic=false;
  parseFlag("ATOMISTIC",atomistic);
  bool martini=false;
  parseFlag("MARTINI",martini);

  if(martini&&atomistic) error("You cannot use martini and atomistic at the same time");

  double rho = 0.334;
  parse("WATERDENS", rho);

  vector<vector<long double> >  FF_tmp;
  FF_tmp.resize(numq,vector<long double>(size));
  if(!atomistic&&!martini) {
    //read in parameter vector
    vector<vector<long double> > parameter;
    parameter.resize(size);
    ntarget=0;
    for(unsigned i=0; i<size; ++i) {
      if( !parseNumberedVector( "PARAMETERS", i+1, parameter[i]) ) break;
      ntarget++;
    }
    if( ntarget!=size ) error("found wrong number of parameter vectors");
    for(unsigned i=0; i<size; ++i) {
      for(unsigned k=0; k<numq; ++k) {
        for(unsigned j=0; j<parameter[i].size(); ++j) {
          FF_tmp[k][i]+= parameter[i][j]*powl(static_cast<long double>(q_list[k]),j);
        }
      }
    }
  } else if(martini) {
    //read in parameter vector
    vector<vector<long double> > parameter;
    parameter.resize(size);
    getMartiniSFparam(atoms, parameter);
    for(unsigned i=0; i<size; ++i) {
      for(unsigned k=0; k<numq; ++k) {
        for(unsigned j=0; j<parameter[i].size(); ++j) {
          FF_tmp[k][i]+= parameter[i][j]*powl(static_cast<long double>(q_list[k]),j);
        }
      }
    }
  } else if(atomistic) {
    calculateASF(atoms, FF_tmp, rho);
  }

  bool exp=false;
  parseFlag("ADDEXP",exp);

  vector<double> expint;
  expint.resize( numq );
  ntarget=0;
  for(unsigned i=0; i<numq; ++i) {
    if( !parseNumbered( "EXPINT", i+1, expint[i] ) ) break;
    ntarget++;
  }
  if( ntarget!=numq && exp==true) error("found wrong number of EXPINT values");

  if(pbc)      log.printf("  using periodic boundary conditions\n");
  else         log.printf("  without periodic boundary conditions\n");
  for(unsigned i=0; i<numq; i++) {
    if(q_list[i]==0.) error("it is not possible to set q=0\n");
    log.printf("  my q: %lf \n",q_list[i]);
  }

  for(unsigned i=0; i<numq; i++) {
    std::string num; Tools::convert(i,num);
    addComponentWithDerivatives("q_"+num);
    componentIsNotPeriodic("q_"+num);
  }
  if(exp) {
    for(unsigned i=0; i<numq; i++) {
      std::string num; Tools::convert(i,num);
      addComponent("exp_"+num);
      componentIsNotPeriodic("exp_"+num);
      Value* comp=getPntrToComponent("exp_"+num);
      comp->set(expint[i]);
    }
  }
  // convert units to nm^-1
  for(unsigned i=0; i<numq; ++i) {
    q_list[i]=q_list[i]*10.0;    //factor 10 to convert from A^-1 to nm^-1
  }

  log<<"  Bibliography ";
  log<<plumed.cite("Jussupow, et al. (in preparation)");
  if(martini)   log<<plumed.cite("Niebling, Björling, Westenhoff, J Appl Crystallogr 47, 1190–1198 (2014).");
  if(atomistic) {
    log<<plumed.cite("Fraser, MacRae, Suzuki, J. Appl. Crystallogr., 11, 693–694 (1978).");
    log<<plumed.cite("Brown, Fox, Maslen, O'Keefe, Willis, International Tables for Crystallography C, 554–595 (International Union of Crystallography, 2006).");
  }
  log<< plumed.cite("Bonomi, Camilloni, Bioinformatics, 33, 3999 (2017)");
  log<<"\n";

  requestAtoms(atoms);
  checkRead();

  // move structure factor to the GPU
  vector<float> FF_new;
  FF_new.resize(numq*size);
  for(unsigned k=0; k<numq; ++k) {
    for(unsigned i=0; i<size; i++) {
      FF_new[k+i*numq] = static_cast<float>(FF_tmp[k][i]/sqrt(scexp));
    }
  }
  af::array allFFa = af::array(numq, size, &FF_new.front());
  FF_value = allFFa;
#endif
}

void SAXSGPU::calculate() {
#ifdef __PLUMED_HAS_ARRAYFIRE
  if(pbc) makeWhole();

  const unsigned size = getNumberOfAtoms();
  const unsigned numq = q_list.size();

  vector<float> posi;
  posi.resize(3*size);
  #pragma omp parallel for num_threads(OpenMP::getNumThreads())
  for (unsigned i=0; i<size; i++) {
    const Vector tmp = getPosition(i);
    posi[3*i]   = static_cast<float>(tmp[0]);
    posi[3*i+1] = static_cast<float>(tmp[1]);
    posi[3*i+2] = static_cast<float>(tmp[2]);
  }

  // create array a and b containing atomic coordinates
  af::setDevice(deviceid);
  // 3,size,1,1
  af::array pos_a = af::array(3, size, &posi.front());
  // size,3,1,1
  pos_a = af::moddims(pos_a.T(), size, 1, 3);
  // size,3,1,1
  af::array pos_b = pos_a(af::span, af::span);
  // size,1,3,1
  pos_a = af::moddims(pos_a, size, 1, 3);
  // 1,size,3,1
  pos_b = af::moddims(pos_b, 1, size, 3);

  // size,size,3,1
  af::array xyz_dist = af::tile(pos_a, 1, size, 1) - af::tile(pos_b, size, 1, 1);
  // size,size,1,1
  af::array square = af::sum(xyz_dist*xyz_dist,2);
  // size,size,1,1
  af::array dist_sqrt = af::sqrt(square);
  // replace the zero of square with one to avoid nan in the derivatives (the number does not matter becasue this are multiplied by zero)
  af::replace(square,!(af::iszero(square)),1.);
  // size,size,3,1
  xyz_dist = xyz_dist / af::tile(square, 1, 1, 3);
  // numq,1,1,1
  af::array sum_device   = af::constant(0, numq, f32);
  // numq,size,3,1
  af::array deriv_device = af::constant(0, numq, size, 3, f32);

  for (unsigned k=0; k<numq; k++) {
    // calculate FF matrix
    // size,size,1,1
    af::array FFdist_mod = af::tile(af::moddims(FF_value(k, af::span), size, 1), 1, size)*af::tile(FF_value(k, af::span), size, 1);

    // get q
    const float qvalue = static_cast<float>(q_list[k]);
    // size,size,1,1
    af::array dist_q = qvalue*dist_sqrt;
    // size,size,1
    af::array dist_sin = af::sin(dist_q)/dist_q;
    af::replace(dist_sin,!(af::isNaN(dist_sin)),1.);
    // 1,1,1,1
    sum_device(k) = af::sum(af::flat(dist_sin)*af::flat(FFdist_mod));

    // size,size,1,1
    af::array tmp = FFdist_mod*(dist_sin - af::cos(dist_q));
    // size,size,3,1
    af::array dd_all = af::tile(tmp, 1, 1, 3)*xyz_dist;
    // it should become 1,size,3
    deriv_device(k, af::span, af::span) = af::sum(dd_all,0);
  }

  // read out results
  std::vector<float> inten; 
  inten.resize(numq);
  sum_device.host(&inten.front());

  std::vector<float> deriv;
  deriv.resize(size*3*numq);
  deriv_device = af::reorder(deriv_device, 2, 1, 0);
  deriv_device = af::flat(deriv_device);
  deriv_device.host(&deriv.front());

  for(unsigned k=0; k<numq; k++) {
    Value* val=getPntrToComponent(k);
    val->set(inten[k]);
    Tensor deriv_box;
    for(unsigned i=0; i<size; i++) {
      const unsigned di = k*size*3+i*3;
      const Vector dd(deriv[di+0],deriv[di+1],deriv[di+2]);
      setAtomsDerivatives(val, i, 2*dd);
      const Vector posi=getPosition(i);
      deriv_box += Tensor(posi,2*dd);
    }
    setBoxDerivatives(val, -deriv_box);
  }
#endif
}

void SAXSGPU::getMartiniSFparam(const vector<AtomNumber> &atoms, vector<vector<long double> > &parameter)
{
  vector<SetupMolInfo*> moldat=plumed.getActionSet().select<SetupMolInfo*>();
  if( moldat.size()==1 ) {
    log<<"  MOLINFO DATA found, using proper atom names\n";
    for(unsigned i=0; i<atoms.size(); ++i) {
      string Aname = moldat[0]->getAtomName(atoms[i]);
      string Rname = moldat[0]->getResidueName(atoms[i]);
      if(Rname=="ALA") {
        if(Aname=="BB") {
          parameter[i].push_back(9.045);
          parameter[i].push_back(-0.098114);
          parameter[i].push_back(7.54281);
          parameter[i].push_back(-1.97438);
          parameter[i].push_back(-8.32689);
          parameter[i].push_back(6.09318);
          parameter[i].push_back(-1.18913);
        } else error("Atom name not known");
      } else if(Rname=="ARG") {
        if(Aname=="BB") {
          parameter[i].push_back(10.729);
          parameter[i].push_back(-0.0392574);
          parameter[i].push_back(1.15382);
          parameter[i].push_back(-0.155999);
          parameter[i].push_back(-2.43619);
          parameter[i].push_back(1.72922);
          parameter[i].push_back(-0.33799);
        } else if(Aname=="SC1") {
          parameter[i].push_back(-2.796);
          parameter[i].push_back(0.472403);
          parameter[i].push_back(8.07424);
          parameter[i].push_back(4.37299);
          parameter[i].push_back(-10.7398);
          parameter[i].push_back(4.95677);
          parameter[i].push_back(-0.725797);
        } else if(Aname=="SC2") {
          parameter[i].push_back(15.396);
          parameter[i].push_back(0.0636736);
          parameter[i].push_back(-1.258);
          parameter[i].push_back(1.93135);
          parameter[i].push_back(-4.45031);
          parameter[i].push_back(2.49356);
          parameter[i].push_back(-0.410721);
        } else error("Atom name not known");
      } else if(Rname=="ASN") {
        if(Aname=="BB") {
          parameter[i].push_back(10.738);
          parameter[i].push_back(-0.0402162);
          parameter[i].push_back(1.03007);
          parameter[i].push_back(-0.254174);
          parameter[i].push_back(-2.12015);
          parameter[i].push_back(1.55535);
          parameter[i].push_back(-0.30963);
        } else if(Aname=="SC1") {
          parameter[i].push_back(9.249);
          parameter[i].push_back(-0.0148678);
          parameter[i].push_back(5.52169);
          parameter[i].push_back(0.00853212);
          parameter[i].push_back(-6.71992);
          parameter[i].push_back(3.93622);
          parameter[i].push_back(-0.64973);
        } else error("Atom name not known");
      } else if(Rname=="ASP") {
        if(Aname=="BB") {
          parameter[i].push_back(10.695);
          parameter[i].push_back(-0.0410247);
          parameter[i].push_back(1.03656);
          parameter[i].push_back(-0.298558);
          parameter[i].push_back(-2.06064);
          parameter[i].push_back(1.53495);
          parameter[i].push_back(-0.308365);
        } else if(Aname=="SC1") {
          parameter[i].push_back(9.476);
          parameter[i].push_back(-0.0254664);
          parameter[i].push_back(5.57899);
          parameter[i].push_back(-0.395027);
          parameter[i].push_back(-5.9407);
          parameter[i].push_back(3.48836);
          parameter[i].push_back(-0.569402);
        } else error("Atom name not known");
      } else if(Rname=="CYS") {
        if(Aname=="BB") {
          parameter[i].push_back(10.698);
          parameter[i].push_back(-0.0233493);
          parameter[i].push_back(1.18257);
          parameter[i].push_back(0.0684464);
          parameter[i].push_back(-2.792);
          parameter[i].push_back(1.88995);
          parameter[i].push_back(-0.360229);
        } else if(Aname=="SC1") {
          parameter[i].push_back(8.199);
          parameter[i].push_back(-0.0261569);
          parameter[i].push_back(6.79677);
          parameter[i].push_back(-0.343845);
          parameter[i].push_back(-5.03578);
          parameter[i].push_back(2.7076);
          parameter[i].push_back(-0.420714);
        } else error("Atom name not known");
      } else if(Rname=="GLN") {
        if(Aname=="BB") {
          parameter[i].push_back(10.728);
          parameter[i].push_back(-0.0391984);
          parameter[i].push_back(1.09264);
          parameter[i].push_back(-0.261555);
          parameter[i].push_back(-2.21245);
          parameter[i].push_back(1.62071);
          parameter[i].push_back(-0.322325);
        } else if(Aname=="SC1") {
          parameter[i].push_back(8.317);
          parameter[i].push_back(-0.229045);
          parameter[i].push_back(12.6338);
          parameter[i].push_back(-7.6719);
          parameter[i].push_back(-5.8376);
          parameter[i].push_back(5.53784);
          parameter[i].push_back(-1.12604);
        } else error("Atom name not known");
      } else if(Rname=="GLU") {
        if(Aname=="BB") {
          parameter[i].push_back(10.694);
          parameter[i].push_back(-0.0521961);
          parameter[i].push_back(1.11153);
          parameter[i].push_back(-0.491995);
          parameter[i].push_back(-1.86236);
          parameter[i].push_back(1.45332);
          parameter[i].push_back(-0.29708);
        } else if(Aname=="SC1") {
          parameter[i].push_back(8.544);
          parameter[i].push_back(-0.249555);
          parameter[i].push_back(12.8031);
          parameter[i].push_back(-8.42696);
          parameter[i].push_back(-4.66486);
          parameter[i].push_back(4.90004);
          parameter[i].push_back(-1.01204);
        } else error("Atom name not known");
      } else if(Rname=="GLY") {
        if(Aname=="BB") {
          parameter[i].push_back(9.977);
          parameter[i].push_back(-0.0285799);
          parameter[i].push_back(1.84236);
          parameter[i].push_back(-0.0315192);
          parameter[i].push_back(-2.88326);
          parameter[i].push_back(1.87323);
          parameter[i].push_back(-0.345773);
        } else error("Atom name not known");
      } else if(Rname=="HIS") {
        if(Aname=="BB") {
          parameter[i].push_back(10.721);
          parameter[i].push_back(-0.0379337);
          parameter[i].push_back(1.06028);
          parameter[i].push_back(-0.236143);
          parameter[i].push_back(-2.17819);
          parameter[i].push_back(1.58357);
          parameter[i].push_back(-0.31345);
        } else if(Aname=="SC1") {
          parameter[i].push_back(-0.424);
          parameter[i].push_back(0.665176);
          parameter[i].push_back(3.4369);
          parameter[i].push_back(2.93795);
          parameter[i].push_back(-5.18288);
          parameter[i].push_back(2.12381);
          parameter[i].push_back(-0.284224);
        } else if(Aname=="SC2") {
          parameter[i].push_back(5.363);
          parameter[i].push_back(-0.0176945);
          parameter[i].push_back(2.9506);
          parameter[i].push_back(-0.387018);
          parameter[i].push_back(-1.83951);
          parameter[i].push_back(0.9703);
          parameter[i].push_back(-0.1458);
        } else if(Aname=="SC3") {
          parameter[i].push_back(5.784);
          parameter[i].push_back(-0.0293129);
          parameter[i].push_back(2.74167);
          parameter[i].push_back(-0.520875);
          parameter[i].push_back(-1.62949);
          parameter[i].push_back(0.902379);
          parameter[i].push_back(-0.139957);
        } else error("Atom name not known");
      } else if(Rname=="ILE") {
        if(Aname=="BB") {
          parameter[i].push_back(10.699);
          parameter[i].push_back(-0.0188962);
          parameter[i].push_back(1.217);
          parameter[i].push_back(0.242481);
          parameter[i].push_back(-3.13898);
          parameter[i].push_back(2.07916);
          parameter[i].push_back(-0.392574);
        } else if(Aname=="SC1") {
          parameter[i].push_back(-4.448);
          parameter[i].push_back(1.20996);
          parameter[i].push_back(11.5141);
          parameter[i].push_back(6.98895);
          parameter[i].push_back(-19.1948);
          parameter[i].push_back(9.89207);
          parameter[i].push_back(-1.60877);
        } else error("Atom name not known");
      } else if(Rname=="LEU") {
        if(Aname=="BB") {
          parameter[i].push_back(10.692);
          parameter[i].push_back(-0.0414917);
          parameter[i].push_back(1.1077);
          parameter[i].push_back(-0.288062);
          parameter[i].push_back(-2.17187);
          parameter[i].push_back(1.59879);
          parameter[i].push_back(-0.318545);
        } else if(Aname=="SC1") {
          parameter[i].push_back(-4.448);
          parameter[i].push_back(2.1063);
          parameter[i].push_back(6.72381);
          parameter[i].push_back(14.6954);
          parameter[i].push_back(-23.7197);
          parameter[i].push_back(10.7247);
          parameter[i].push_back(-1.59146);
        } else error("Atom name not known");
      } else if(Rname=="LYS") {
        if(Aname=="BB") {
          parameter[i].push_back(10.706);
          parameter[i].push_back(-0.0468629);
          parameter[i].push_back(1.09477);
          parameter[i].push_back(-0.432751);
          parameter[i].push_back(-1.94335);
          parameter[i].push_back(1.49109);
          parameter[i].push_back(-0.302589);
        } else if(Aname=="SC1") {
          parameter[i].push_back(-2.796);
          parameter[i].push_back(0.508044);
          parameter[i].push_back(7.91436);
          parameter[i].push_back(4.54097);
          parameter[i].push_back(-10.8051);
          parameter[i].push_back(4.96204);
          parameter[i].push_back(-0.724414);
        } else if(Aname=="SC2") {
          parameter[i].push_back(3.070);
          parameter[i].push_back(-0.0101448);
          parameter[i].push_back(4.67994);
          parameter[i].push_back(-0.792529);
          parameter[i].push_back(-2.09142);
          parameter[i].push_back(1.02933);
          parameter[i].push_back(-0.137787);
        } else error("Atom name not known");
      } else if(Rname=="MET") {
        if(Aname=="BB") {
          parameter[i].push_back(10.671);
          parameter[i].push_back(-0.0433724);
          parameter[i].push_back(1.13784);
          parameter[i].push_back(-0.40768);
          parameter[i].push_back(-2.00555);
          parameter[i].push_back(1.51673);
          parameter[i].push_back(-0.305547);
        } else if(Aname=="SC1") {
          parameter[i].push_back(5.85);
          parameter[i].push_back(-0.0485798);
          parameter[i].push_back(17.0391);
          parameter[i].push_back(-3.65327);
          parameter[i].push_back(-13.174);
          parameter[i].push_back(8.68286);
          parameter[i].push_back(-1.56095);
        } else error("Atom name not known");
      } else if(Rname=="PHE") {
        if(Aname=="BB") {
          parameter[i].push_back(10.741);
          parameter[i].push_back(-0.0317275);
          parameter[i].push_back(1.15599);
          parameter[i].push_back(0.0276187);
          parameter[i].push_back(-2.74757);
          parameter[i].push_back(1.88783);
          parameter[i].push_back(-0.363525);
        } else if(Aname=="SC1") {
          parameter[i].push_back(-0.636);
          parameter[i].push_back(0.527882);
          parameter[i].push_back(6.77612);
          parameter[i].push_back(3.18508);
          parameter[i].push_back(-8.92826);
          parameter[i].push_back(4.29752);
          parameter[i].push_back(-0.65187);
        } else if(Aname=="SC2") {
          parameter[i].push_back(-0.424);
          parameter[i].push_back(0.389174);
          parameter[i].push_back(4.11761);
          parameter[i].push_back(2.29527);
          parameter[i].push_back(-4.7652);
          parameter[i].push_back(1.97023);
          parameter[i].push_back(-0.262318);
        } else if(Aname=="SC3") {
          parameter[i].push_back(-0.424);
          parameter[i].push_back(0.38927);
          parameter[i].push_back(4.11708);
          parameter[i].push_back(2.29623);
          parameter[i].push_back(-4.76592);
          parameter[i].push_back(1.97055);
          parameter[i].push_back(-0.262381);
        } else error("Atom name not known");
      } else if(Rname=="PRO") {
        if(Aname=="BB") {
          parameter[i].push_back(11.434);
          parameter[i].push_back(-0.033323);
          parameter[i].push_back(0.472014);
          parameter[i].push_back(-0.290854);
          parameter[i].push_back(-1.81409);
          parameter[i].push_back(1.39751);
          parameter[i].push_back(-0.280407);
        } else if(Aname=="SC1") {
          parameter[i].push_back(-2.796);
          parameter[i].push_back(0.95668);
          parameter[i].push_back(6.84197);
          parameter[i].push_back(6.43774);
          parameter[i].push_back(-12.5068);
          parameter[i].push_back(5.64597);
          parameter[i].push_back(-0.825206);
        } else error("Atom name not known");
      } else if(Rname=="SER") {
        if(Aname=="BB") {
          parameter[i].push_back(10.699);
          parameter[i].push_back(-0.0325828);
          parameter[i].push_back(1.20329);
          parameter[i].push_back(-0.0674351);
          parameter[i].push_back(-2.60749);
          parameter[i].push_back(1.80318);
          parameter[i].push_back(-0.346803);
        } else if(Aname=="SC1") {
          parameter[i].push_back(3.298);
          parameter[i].push_back(-0.0366801);
          parameter[i].push_back(5.11077);
          parameter[i].push_back(-1.46774);
          parameter[i].push_back(-1.48421);
          parameter[i].push_back(0.800326);
          parameter[i].push_back(-0.108314);
        } else error("Atom name not known");
      } else if(Rname=="THR") {
        if(Aname=="BB") {
          parameter[i].push_back(10.697);
          parameter[i].push_back(-0.0242955);
          parameter[i].push_back(1.24671);
          parameter[i].push_back(0.146423);
          parameter[i].push_back(-2.97429);
          parameter[i].push_back(1.97513);
          parameter[i].push_back(-0.371479);
        } else if(Aname=="SC1") {
          parameter[i].push_back(2.366);
          parameter[i].push_back(0.0297604);
          parameter[i].push_back(11.9216);
          parameter[i].push_back(-9.32503);
          parameter[i].push_back(1.9396);
          parameter[i].push_back(0.0804861);
          parameter[i].push_back(-0.0302721);
        } else error("Atom name not known");
      } else if(Rname=="TRP") {
        if(Aname=="BB") {
          parameter[i].push_back(10.689);
          parameter[i].push_back(-0.0265879);
          parameter[i].push_back(1.17819);
          parameter[i].push_back(0.0386457);
          parameter[i].push_back(-2.75634);
          parameter[i].push_back(1.88065);
          parameter[i].push_back(-0.360217);
        } else if(Aname=="SC1") {
          parameter[i].push_back(0.084);
          parameter[i].push_back(0.752407);
          parameter[i].push_back(5.3802);
          parameter[i].push_back(4.09281);
          parameter[i].push_back(-9.28029);
          parameter[i].push_back(4.45923);
          parameter[i].push_back(-0.689008);
        } else if(Aname=="SC2") {
          parameter[i].push_back(5.739);
          parameter[i].push_back(0.0298492);
          parameter[i].push_back(4.60446);
          parameter[i].push_back(1.34463);
          parameter[i].push_back(-5.69968);
          parameter[i].push_back(2.84924);
          parameter[i].push_back(-0.433781);
        } else if(Aname=="SC3") {
          parameter[i].push_back(-0.424);
          parameter[i].push_back(0.388576);
          parameter[i].push_back(4.11859);
          parameter[i].push_back(2.29485);
          parameter[i].push_back(-4.76255);
          parameter[i].push_back(1.96849);
          parameter[i].push_back(-0.262015);
        } else if(Aname=="SC4") {
          parameter[i].push_back(-0.424);
          parameter[i].push_back(0.387685);
          parameter[i].push_back(4.12153);
          parameter[i].push_back(2.29144);
          parameter[i].push_back(-4.7589);
          parameter[i].push_back(1.96686);
          parameter[i].push_back(-0.261786);
        } else error("Atom name not known");
      } else if(Rname=="TYR") {
        if(Aname=="BB") {
          parameter[i].push_back(10.689);
          parameter[i].push_back(-0.0193526);
          parameter[i].push_back(1.18241);
          parameter[i].push_back(0.207318);
          parameter[i].push_back(-3.0041);
          parameter[i].push_back(1.99335);
          parameter[i].push_back(-0.376482);
        } else if(Aname=="SC1") {
          parameter[i].push_back(-0.636);
          parameter[i].push_back(0.528902);
          parameter[i].push_back(6.78168);
          parameter[i].push_back(3.17769);
          parameter[i].push_back(-8.93667);
          parameter[i].push_back(4.30692);
          parameter[i].push_back(-0.653993);
        } else if(Aname=="SC2") {
          parameter[i].push_back(-0.424);
          parameter[i].push_back(0.388811);
          parameter[i].push_back(4.11851);
          parameter[i].push_back(2.29545);
          parameter[i].push_back(-4.7668);
          parameter[i].push_back(1.97131);
          parameter[i].push_back(-0.262534);
        } else if(Aname=="SC3") {
          parameter[i].push_back(4.526);
          parameter[i].push_back(-0.00381305);
          parameter[i].push_back(5.8567);
          parameter[i].push_back(-0.214086);
          parameter[i].push_back(-4.63649);
          parameter[i].push_back(2.52869);
          parameter[i].push_back(-0.39894);
        } else error("Atom name not known");
      } else if(Rname=="VAL") {
        if(Aname=="BB") {
          parameter[i].push_back(10.691);
          parameter[i].push_back(-0.0162929);
          parameter[i].push_back(1.24446);
          parameter[i].push_back(0.307914);
          parameter[i].push_back(-3.27446);
          parameter[i].push_back(2.14788);
          parameter[i].push_back(-0.403259);
        } else if(Aname=="SC1") {
          parameter[i].push_back(-3.516);
          parameter[i].push_back(1.62307);
          parameter[i].push_back(5.43064);
          parameter[i].push_back(9.28809);
          parameter[i].push_back(-14.9927);
          parameter[i].push_back(6.6133);
          parameter[i].push_back(-0.964977);
        } else error("Atom name not known");
      } else error("Residue not known");
    }
  } else {
    error("MOLINFO DATA not found\n");
  }
}

void SAXSGPU::calculateASF(const vector<AtomNumber> &atoms, vector<vector<long double> > &FF_tmp, const double rho)
{
  enum { H, C, N, O, P, S, NTT };
  map<string, unsigned> AA_map;
  AA_map["H"] = H;
  AA_map["C"] = C;
  AA_map["N"] = N;
  AA_map["O"] = O;
  AA_map["P"] = P;
  AA_map["S"] = S;

  vector<vector<double> > param_a;
  vector<vector<double> > param_b;
  vector<double> param_c;
  vector<double> param_v;

  param_a.resize(NTT, vector<double>(5));
  param_b.resize(NTT, vector<double>(5));
  param_c.resize(NTT);
  param_v.resize(NTT);

  param_a[H][0] = 0.493002; param_b[H][0] = 10.5109; param_c[H] = 0.003038;
  param_a[H][1] = 0.322912; param_b[H][1] = 26.1257; param_v[H] = 5.15;
  param_a[H][2] = 0.140191; param_b[H][2] = 3.14236;
  param_a[H][3] = 0.040810; param_b[H][3] = 57.7997;
  param_a[H][4] = 0.0;      param_b[H][4] = 1.0;

  param_a[C][0] = 2.31000; param_b[C][0] = 20.8439; param_c[C] = 0.215600;
  param_a[C][1] = 1.02000; param_b[C][1] = 10.2075; param_v[C] = 16.44;
  param_a[C][2] = 1.58860; param_b[C][2] = 0.56870;
  param_a[C][3] = 0.86500; param_b[C][3] = 51.6512;
  param_a[C][4] = 0.0;     param_b[C][4] = 1.0;

  param_a[N][0] = 12.2126; param_b[N][0] = 0.00570; param_c[N] = -11.529;
  param_a[N][1] = 3.13220; param_b[N][1] = 9.89330; param_v[N] = 2.49;
  param_a[N][2] = 2.01250; param_b[N][2] = 28.9975;
  param_a[N][3] = 1.16630; param_b[N][3] = 0.58260;
  param_a[N][4] = 0.0;     param_b[N][4] = 1.0;

  param_a[O][0] = 3.04850; param_b[O][0] = 13.2771; param_c[O] = 0.250800 ;
  param_a[O][1] = 2.28680; param_b[O][1] = 5.70110; param_v[O] = 9.13;
  param_a[O][2] = 1.54630; param_b[O][2] = 0.32390;
  param_a[O][3] = 0.86700; param_b[O][3] = 32.9089;
  param_a[O][4] = 0.0;     param_b[O][4] = 1.0;

  param_a[P][0] = 6.43450; param_b[P][0] = 1.90670; param_c[P] = 1.11490;
  param_a[P][1] = 4.17910; param_b[P][1] = 27.1570; param_v[P] = 5.73;
  param_a[P][2] = 1.78000; param_b[P][2] = 0.52600;
  param_a[P][3] = 1.49080; param_b[P][3] = 68.1645;
  param_a[P][4] = 0.0;     param_b[P][4] = 1.0;

  param_a[S][0] = 6.90530; param_b[S][0] = 1.46790; param_c[S] = 0.866900;
  param_a[S][1] = 5.20340; param_b[S][1] = 22.2151; param_v[S] = 19.86;
  param_a[S][2] = 1.43790; param_b[S][2] = 0.25360;
  param_a[S][3] = 1.58630; param_b[S][3] = 56.1720;
  param_a[S][4] = 0.0;     param_b[S][4] = 1.0;

  vector<SetupMolInfo*> moldat=plumed.getActionSet().select<SetupMolInfo*>();

  if( moldat.size()==1 ) {
    log<<"  MOLINFO DATA found, using proper atom names\n";
    for(unsigned i=0; i<atoms.size(); ++i) {
      // get atom name
      string name = moldat[0]->getAtomName(atoms[i]);
      char type;
      // get atom type
      char first = name.at(0);
      // GOLDEN RULE: type is first letter, if not a number
      if (!isdigit(first)) {
        type = first;
        // otherwise is the second
      } else {
        type = name.at(1);
      }
      std::string type_s = std::string(1,type);
      if(AA_map.find(type_s) != AA_map.end()) {
        const unsigned index=AA_map[type_s];
        const double volr = pow(param_v[index], (2.0/3.0)) /(4. * M_PI);
        for(unsigned k=0; k<q_list.size(); ++k) {
          const double q = q_list[k];
          const double s = q / (4. * M_PI);
          FF_tmp[k][i] = param_c[index];
          // SUM [a_i * EXP( - b_i * (q/4pi)^2 )] Waasmaier and Kirfel (1995)
          for(unsigned j=0; j<4; j++) {
            FF_tmp[k][i] += param_a[index][j]*exp(-param_b[index][j]*s*s);
          }
          // subtract solvation: rho * v_i * EXP( (- v_i^(2/3) / (4pi)) * q^2  ) // since  D in Fraser 1978 is 2*s
          FF_tmp[k][i] -= rho*param_v[index]*exp(-volr*q*q);
        }
      } else {
        error("Wrong atom type "+type_s+" from atom name "+name+"\n");
      }
    }
  } else {
    error("MOLINFO DATA not found\n");
  }
}

}
}
