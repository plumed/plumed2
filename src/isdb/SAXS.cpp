/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2017-2020 The plumed team
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
 This class was originally written by Alexander Jussupow
 Extension for the middleman algorithm by Max Muehlbauer
 Martini beads strucutre factor for Nucleic Acids by Cristina Paissoni
*/

#include "MetainferenceBase.h"
#include "core/ActionRegister.h"
#include "core/ActionSet.h"
#include "core/SetupMolInfo.h"
#include "tools/Communicator.h"
#include "tools/Pbc.h"

#include <string>
#include <cmath>
#include <map>

#ifdef __PLUMED_HAS_ARRAYFIRE
#include <arrayfire.h>
#include <af/util.h>
#endif

#ifndef M_PI
#define M_PI           3.14159265358979323846
#endif

using namespace std;

namespace PLMD {
namespace isdb {

//+PLUMEDOC ISDB_COLVAR SAXS
/*
Calculates SAXS scattered intensity using either the Debye equation.

Intensities are calculated for a set of scattering length set using QVALUE keywords that are numbered starting from 0.
Structure factors can be either assigned using a polynomial expansion to any order using the PARAMETERS keywords;
automatically assigned to atoms using the ATOMISTIC flag reading a PDB file, a correction for the water density is
automatically added, with water density that by default is 0.334 but that can be set otherwise using WATERDENS;
automatically assigned to Martini pseudo atoms using the MARTINI flag.
The calculated intensities can be scaled using the SCALEINT keywords. This is applied by rescaling the structure factors.
Experimental reference intensities can be added using the EXPINT keywords.
By default SAXS is calculated using Debye on CPU, by adding the GPU flag it is possible to solve the equation on a GPU
if the ARRAYFIRE libraries are installed and correctly linked.
\ref METAINFERENCE can be activated using DOSCORE and the other relevant keywords.

\par Examples
in the following example the saxs intensities for a martini model are calculated. structure factors
are obtained from the pdb file indicated in the MOLINFO.

\plumedfile
MOLINFO STRUCTURE=template.pdb

SAXS ...
LABEL=saxs
ATOMS=1-355
SCALEINT=3920000
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
... SAXS

PRINT ARG=(saxs\.q-.*),(saxs\.exp-.*) FILE=colvar STRIDE=1

\endplumedfile

*/
//+ENDPLUMEDOC

class SAXS :
  public MetainferenceBase
{
private:
  enum { H, C, N, O, P, S, NTT };
  enum { ALA_BB, ARG_BB, ARG_SC1, ARG_SC2, ASN_BB, ASN_SC1, ASP_BB, ASP_SC1, CYS_BB, CYS_SC1,
         GLN_BB, GLN_SC1, GLU_BB, GLU_SC1, GLY_BB, HIS_BB, HIS_SC1, HIS_SC2, HIS_SC3, ILE_BB,
         ILE_SC1, LEU_BB, LEU_SC1, LYS_BB, LYS_SC1, LYS_SC2, MET_BB, MET_SC1, PHE_BB, PHE_SC1,
         PHE_SC2, PHE_SC3, PRO_BB, PRO_SC1, SER_BB, SER_SC1, THR_BB, THR_SC1, TRP_BB, TRP_SC1,
         TRP_SC2, TRP_SC3, TRP_SC4, TYR_BB, TYR_SC1, TYR_SC2, TYR_SC3, VAL_BB, VAL_SC1, A_BB1,
         A_BB2, A_BB3, A_SC1, A_SC2, A_SC3, A_SC4, A_3TE, A_5TE, A_TE3, A_TE5, C_BB1, C_BB2,
         C_BB3, C_SC1, C_SC2, C_SC3, C_3TE, C_5TE, C_TE3, C_TE5, G_BB1, G_BB2, G_BB3, G_SC1,
         G_SC2, G_SC3, G_SC4, G_3TE, G_5TE, G_TE3, G_TE5, U_BB1, U_BB2, U_BB3, U_SC1, U_SC2,
         U_SC3, U_3TE, U_5TE, U_TE3, U_TE5, DA_BB1, DA_BB2, DA_BB3, DA_SC1, DA_SC2, DA_SC3,
         DA_SC4, DA_3TE, DA_5TE, DA_TE3, DA_TE5, DC_BB1, DC_BB2, DC_BB3, DC_SC1, DC_SC2, DC_SC3,
         DC_3TE, DC_5TE, DC_TE3, DC_TE5, DG_BB1, DG_BB2, DG_BB3, DG_SC1, DG_SC2, DG_SC3, DG_SC4,
         DG_3TE, DG_5TE, DG_TE3, DG_TE5, DT_BB1, DT_BB2, DT_BB3, DT_SC1, DT_SC2, DT_SC3, DT_3TE,
         DT_5TE, DT_TE3, DT_TE5, NMARTINI
       };
  bool                       pbc;
  bool                       serial;
  bool                       gpu;
  int                        deviceid;
  vector<unsigned>           atoi;
  vector<double>             q_list;
  vector<double>             FF_rank;
  vector<vector<double> >    FF_value;
  vector<vector<float> >     FFf_value;

  void calculate_gpu(vector<Vector> &deriv);
  void calculate_cpu(vector<Vector> &deriv);
  void getMartiniSFparam(const vector<AtomNumber> &atoms, vector<vector<long double> > &parameter);
  double calculateASF(const vector<AtomNumber> &atoms, vector<vector<long double> > &FF_tmp, const double rho);

public:
  static void registerKeywords( Keywords& keys );
  explicit SAXS(const ActionOptions&);
  void calculate() override;
  void update() override;
};

PLUMED_REGISTER_ACTION(SAXS,"SAXS")

void SAXS::registerKeywords(Keywords& keys) {
  componentsAreNotOptional(keys);
  MetainferenceBase::registerKeywords(keys);
  keys.addFlag("NOPBC",false,"ignore the periodic boundary conditions when calculating distances");
  keys.addFlag("SERIAL",false,"Perform the calculation in serial - for debug purpose");
  keys.add("compulsory","DEVICEID","0","Identifier of the GPU to be used");
  keys.addFlag("GPU",false,"calculate SAXS using ARRAYFIRE on an accelerator device");
  keys.addFlag("ATOMISTIC",false,"calculate SAXS for an atomistic model");
  keys.addFlag("MARTINI",false,"calculate SAXS for a Martini model");
  keys.add("atoms","ATOMS","The atoms to be included in the calculation, e.g. the whole protein.");
  keys.add("numbered","QVALUE","Selected scattering lengths in Angstrom are given as QVALUE1, QVALUE2, ... .");
  keys.add("numbered","PARAMETERS","Used parameter Keywords like PARAMETERS1, PARAMETERS2. These are used to calculate the structure factor for the \\f$i\\f$th atom/bead.");
  keys.add("compulsory","WATERDENS","0.334","Density of the water to be used for the correction of atomistic structure factors.");
  keys.add("numbered","EXPINT","Add an experimental value for each q value.");
  keys.add("compulsory","SCALEINT","1.0","SCALING value of the calculated data. Useful to simplify the comparison.");
  keys.addOutputComponent("q","default","the # SAXS of q");
  keys.addOutputComponent("exp","EXPINT","the # experimental intensity");
}

SAXS::SAXS(const ActionOptions&ao):
  PLUMED_METAINF_INIT(ao),
  pbc(true),
  serial(false),
  gpu(false),
  deviceid(0)
{
  vector<AtomNumber> atoms;
  parseAtomList("ATOMS",atoms);
  const unsigned size = atoms.size();

  parseFlag("SERIAL",serial);

  bool nopbc=!pbc;
  parseFlag("NOPBC",nopbc);
  pbc=!nopbc;
  if(pbc)      log.printf("  using periodic boundary conditions\n");
  else         log.printf("  without periodic boundary conditions\n");

  parseFlag("GPU",gpu);
#ifndef  __PLUMED_HAS_ARRAYFIRE
  if(gpu) error("To use the GPU mode PLUMED must be compiled with ARRAYFIRE");
#endif

  parse("DEVICEID",deviceid);
#ifdef  __PLUMED_HAS_ARRAYFIRE
  if(gpu) {
    af::setDevice(deviceid);
    af::info();
  }
#endif

  unsigned ntarget=0;
  for(unsigned i=0;; ++i) {
    double t_list;
    if( !parseNumbered( "QVALUE", i+1, t_list) ) break;
    if(t_list<=0.) error("QVALUE cannot be less or equal to zero!\n");
    q_list.push_back(t_list);
    ntarget++;
  }
  const unsigned numq = ntarget;

  for(unsigned i=0; i<numq; i++) {
    if(q_list[i]==0.) error("it is not possible to set q=0\n");
    if(i>0&&q_list[i]<q_list[i-1]) error("QVALUE must be in ascending order");
    log.printf("  my q: %lf \n",q_list[i]);
  }

  bool atomistic=false;
  parseFlag("ATOMISTIC",atomistic);
  bool martini=false;
  parseFlag("MARTINI",martini);

  if(martini&&atomistic) error("You cannot use MARTINI and ATOMISTIC at the same time");

  double rho = 0.334;
  parse("WATERDENS", rho);

  double Iq0=0;
  vector<vector<long double> > FF_tmp;
  atoi.resize(size);
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
    FF_tmp.resize(numq,vector<long double>(size));
    for(unsigned i=0; i<size; ++i) {
      atoi[i]=i;
      for(unsigned k=0; k<numq; ++k) {
        for(unsigned j=0; j<parameter[i].size(); ++j) {
          FF_tmp[k][i]+= parameter[i][j]*powl(static_cast<long double>(q_list[k]),j);
        }
      }
    }
    for(unsigned i=0; i<size; ++i) Iq0+=parameter[i][0];
  } else if(martini) {
    //read in parameter vector
    FF_tmp.resize(numq,vector<long double>(NMARTINI));
    vector<vector<long double> > parameter;
    parameter.resize(NMARTINI);
    getMartiniSFparam(atoms, parameter);
    for(unsigned i=0; i<NMARTINI; ++i) {
      for(unsigned k=0; k<numq; ++k) {
        for(unsigned j=0; j<parameter[i].size(); ++j) {
          FF_tmp[k][i]+= parameter[i][j]*powl(static_cast<long double>(q_list[k]),j);
        }
      }
    }
    for(unsigned i=0; i<size; ++i) Iq0+=parameter[atoi[i]][0];
  } else if(atomistic) {
    FF_tmp.resize(numq,vector<long double>(NTT));
    Iq0=calculateASF(atoms, FF_tmp, rho);
  }
  double scale_int = Iq0*Iq0;

  vector<double> expint;
  expint.resize( numq );
  ntarget=0;
  for(unsigned i=0; i<numq; ++i) {
    if( !parseNumbered( "EXPINT", i+1, expint[i] ) ) break;
    ntarget++;
  }
  bool exp=false;
  if(ntarget!=numq && ntarget!=0) error("found wrong number of EXPINT values");
  if(ntarget==numq) exp=true;
  if(getDoScore()&&!exp) error("with DOSCORE you need to set the EXPINT values");

  double tmp_scale_int=1.;
  parse("SCALEINT",tmp_scale_int);

  if(tmp_scale_int!=1) scale_int /= tmp_scale_int;
  else {
    if(exp) scale_int /= expint[0];
  }

  if(!gpu) {
    FF_rank.resize(numq);
    unsigned n_atom_types=size;
    if(atomistic) n_atom_types=NTT;
    else if(martini) n_atom_types=NMARTINI;
    FF_value.resize(n_atom_types,vector<double>(numq));
    for(unsigned k=0; k<numq; ++k) {
      for(unsigned i=0; i<n_atom_types; i++) FF_value[i][k] = static_cast<double>(FF_tmp[k][i])/sqrt(scale_int);
      for(unsigned i=0; i<size; i++) FF_rank[k] += FF_value[atoi[i]][k]*FF_value[atoi[i]][k];
    }
  } else {
    FFf_value.resize(numq,vector<float>(size));
    for(unsigned k=0; k<numq; ++k) {
      for(unsigned i=0; i<size; i++) {
        FFf_value[k][i] = static_cast<float>(FF_tmp[k][atoi[i]])/sqrt(scale_int);
      }
    }
  }

  if(!getDoScore()) {
    for(unsigned i=0; i<numq; i++) {
      std::string num; Tools::convert(i,num);
      addComponentWithDerivatives("q-"+num);
      componentIsNotPeriodic("q-"+num);
    }
    if(exp) {
      for(unsigned i=0; i<numq; i++) {
        std::string num; Tools::convert(i,num);
        addComponent("exp-"+num);
        componentIsNotPeriodic("exp-"+num);
        Value* comp=getPntrToComponent("exp-"+num);
        comp->set(expint[i]);
      }
    }
  } else {
    for(unsigned i=0; i<numq; i++) {
      std::string num; Tools::convert(i,num);
      addComponent("q-"+num);
      componentIsNotPeriodic("q-"+num);
    }
    for(unsigned i=0; i<numq; i++) {
      std::string num; Tools::convert(i,num);
      addComponent("exp-"+num);
      componentIsNotPeriodic("exp-"+num);
      Value* comp=getPntrToComponent("exp-"+num);
      comp->set(expint[i]);
    }
  }

  // convert units to nm^-1
  for(unsigned i=0; i<numq; ++i) {
    q_list[i]=q_list[i]*10.0;    //factor 10 to convert from A^-1 to nm^-1
  }
  log<<"  Bibliography ";
  if(martini) {
    log<<plumed.cite("Niebling, Björling, Westenhoff, J Appl Crystallogr 47, 1190–1198 (2014).");
    log<<plumed.cite("Paissoni, Jussupow, Camilloni, J Appl Crystallogr 52, 394-402 (2019).");
  }
  if(atomistic) {
    log<<plumed.cite("Fraser, MacRae, Suzuki, J. Appl. Crystallogr., 11, 693–694 (1978).");
    log<<plumed.cite("Brown, Fox, Maslen, O'Keefe, Willis, International Tables for Crystallography C, 554–595 (International Union of Crystallography, 2006).");
  }
  log<< plumed.cite("Bonomi, Camilloni, Bioinformatics, 33, 3999 (2017)");
  log<<"\n";

  requestAtoms(atoms, false);
  if(getDoScore()) {
    setParameters(expint);
    Initialise(numq);
  }
  setDerivatives();
  checkRead();
}

void SAXS::calculate_gpu(vector<Vector> &deriv)
{
#ifdef __PLUMED_HAS_ARRAYFIRE
  const unsigned size = getNumberOfAtoms();
  const unsigned numq = q_list.size();

  std::vector<float> sum;
  sum.resize(numq);

  std::vector<float> dd;
  dd.resize(size*3*numq);

  // on gpu only the master rank run the calculation
  if(comm.Get_rank()==0) {
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
      // size,1,1,1
      af::array AFF_value(size, &FFf_value[k].front());
      // size,size,1,1
      af::array FFdist_mod = af::tile(AFF_value(af::span), 1, size)*af::transpose(af::tile(AFF_value(af::span), 1, size));

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
    sum_device.host(&sum.front());

    deriv_device = af::reorder(deriv_device, 2, 1, 0);
    deriv_device = af::flat(deriv_device);
    deriv_device.host(&dd.front());
  }

  comm.Bcast(dd, 0);
  comm.Bcast(sum, 0);

  for(unsigned k=0; k<numq; k++) {
    string num; Tools::convert(k,num);
    Value* val=getPntrToComponent("q-"+num);
    val->set(sum[k]);
    if(getDoScore()) setCalcData(k, sum[k]);
    for(unsigned i=0; i<size; i++) {
      const unsigned di = k*size*3+i*3;
      deriv[k*size+i] = Vector(2.*dd[di+0],2.*dd[di+1],2.*dd[di+2]);
    }
  }
#endif
}

void SAXS::calculate_cpu(vector<Vector> &deriv)
{
  const unsigned size = getNumberOfAtoms();
  const unsigned numq = q_list.size();

  unsigned stride = comm.Get_size();
  unsigned rank   = comm.Get_rank();
  if(serial) {
    stride = 1;
    rank   = 0;
  }

  vector<double> sum(numq,0);
  unsigned nt=OpenMP::getNumThreads();
  #pragma omp parallel num_threads(nt)
  {
    std::vector<Vector> omp_deriv(deriv.size());
    std::vector<double> omp_sum(numq,0);
    #pragma omp for nowait
    for (unsigned i=rank; i<size-1; i+=stride) {
      Vector posi = getPosition(i);
      for (unsigned j=i+1; j<size ; j++) {
        Vector c_distances = delta(posi,getPosition(j));
        double m_distances = c_distances.modulo();
        c_distances = c_distances/m_distances/m_distances;
        for (unsigned k=0; k<numq; k++) {
          unsigned kdx=k*size;
          double qdist = q_list[k]*m_distances;
          double FFF = 2.*FF_value[atoi[i]][k]*FF_value[atoi[j]][k];
          double tsq = sin(qdist)/qdist;
          double tcq = cos(qdist);
          double tmp = FFF*(tcq-tsq);
          Vector dd  = c_distances*tmp;
          if(nt>1) {
            omp_deriv[kdx+i] -=dd;
            omp_deriv[kdx+j] +=dd;
            omp_sum[k] += FFF*tsq;
          } else {
            deriv[kdx+i] -= dd;
            deriv[kdx+j] += dd;
            sum[k]       += FFF*tsq;
          }
        }
      }
    }
    #pragma omp critical
    if(nt>1) {
      for(unsigned i=0; i<deriv.size(); i++) deriv[i]+=omp_deriv[i];
      for(unsigned k=0; k<numq; k++) sum[k]+=omp_sum[k];
    }
  }

  if(!serial) {
    comm.Sum(&deriv[0][0], 3*deriv.size());
    comm.Sum(&sum[0], numq);
  }

  for (unsigned k=0; k<numq; k++) {
    sum[k]+=FF_rank[k];
    string num; Tools::convert(k,num);
    Value* val=getPntrToComponent("q-"+num);
    val->set(sum[k]);
    if(getDoScore()) setCalcData(k, sum[k]);
  }
}

void SAXS::calculate()
{
  if(pbc) makeWhole();

  const unsigned size = getNumberOfAtoms();
  const unsigned numq = q_list.size();

  vector<Vector> deriv(numq*size);
  if(gpu) calculate_gpu(deriv);
  else calculate_cpu(deriv);

  if(getDoScore()) {
    /* Metainference */
    double score = getScore();
    setScore(score);
  }

  for (unsigned k=0; k<numq; k++) {
    const unsigned kdx=k*size;
    Tensor deriv_box;
    Value* val;
    if(!getDoScore()) {
      string num; Tools::convert(k,num);
      val=getPntrToComponent("q-"+num);
      for(unsigned i=0; i<size; i++) {
        setAtomsDerivatives(val, i, deriv[kdx+i]);
        deriv_box += Tensor(getPosition(i),deriv[kdx+i]);
      }
    } else {
      val=getPntrToComponent("score");
      for(unsigned i=0; i<size; i++) {
        setAtomsDerivatives(val, i, deriv[kdx+i]*getMetaDer(k));
        deriv_box += Tensor(getPosition(i),deriv[kdx+i]*getMetaDer(k));
      }
    }
    setBoxDerivatives(val, -deriv_box);
  }
}

void SAXS::update() {
  // write status file
  if(getWstride()>0&& (getStep()%getWstride()==0 || getCPT()) ) writeStatus();
}

void SAXS::getMartiniSFparam(const vector<AtomNumber> &atoms, vector<vector<long double> > &parameter)
{
  parameter[ALA_BB].push_back(9.045);
  parameter[ALA_BB].push_back(-0.098114);
  parameter[ALA_BB].push_back(7.54281);
  parameter[ALA_BB].push_back(-1.97438);
  parameter[ALA_BB].push_back(-8.32689);
  parameter[ALA_BB].push_back(6.09318);
  parameter[ALA_BB].push_back(-1.18913);

  parameter[ARG_BB].push_back(10.729);
  parameter[ARG_BB].push_back(-0.0392574);
  parameter[ARG_BB].push_back(1.15382);
  parameter[ARG_BB].push_back(-0.155999);
  parameter[ARG_BB].push_back(-2.43619);
  parameter[ARG_BB].push_back(1.72922);
  parameter[ARG_BB].push_back(-0.33799);

  parameter[ARG_SC1].push_back(-2.796);
  parameter[ARG_SC1].push_back(0.472403);
  parameter[ARG_SC1].push_back(8.07424);
  parameter[ARG_SC1].push_back(4.37299);
  parameter[ARG_SC1].push_back(-10.7398);
  parameter[ARG_SC1].push_back(4.95677);
  parameter[ARG_SC1].push_back(-0.725797);

  parameter[ARG_SC2].push_back(15.396);
  parameter[ARG_SC2].push_back(0.0636736);
  parameter[ARG_SC2].push_back(-1.258);
  parameter[ARG_SC2].push_back(1.93135);
  parameter[ARG_SC2].push_back(-4.45031);
  parameter[ARG_SC2].push_back(2.49356);
  parameter[ARG_SC2].push_back(-0.410721);

  parameter[ASN_BB].push_back(10.738);
  parameter[ASN_BB].push_back(-0.0402162);
  parameter[ASN_BB].push_back(1.03007);
  parameter[ASN_BB].push_back(-0.254174);
  parameter[ASN_BB].push_back(-2.12015);
  parameter[ASN_BB].push_back(1.55535);
  parameter[ASN_BB].push_back(-0.30963);

  parameter[ASN_SC1].push_back(9.249);
  parameter[ASN_SC1].push_back(-0.0148678);
  parameter[ASN_SC1].push_back(5.52169);
  parameter[ASN_SC1].push_back(0.00853212);
  parameter[ASN_SC1].push_back(-6.71992);
  parameter[ASN_SC1].push_back(3.93622);
  parameter[ASN_SC1].push_back(-0.64973);

  parameter[ASP_BB].push_back(10.695);
  parameter[ASP_BB].push_back(-0.0410247);
  parameter[ASP_BB].push_back(1.03656);
  parameter[ASP_BB].push_back(-0.298558);
  parameter[ASP_BB].push_back(-2.06064);
  parameter[ASP_BB].push_back(1.53495);
  parameter[ASP_BB].push_back(-0.308365);

  parameter[ASP_SC1].push_back(9.476);
  parameter[ASP_SC1].push_back(-0.0254664);
  parameter[ASP_SC1].push_back(5.57899);
  parameter[ASP_SC1].push_back(-0.395027);
  parameter[ASP_SC1].push_back(-5.9407);
  parameter[ASP_SC1].push_back(3.48836);
  parameter[ASP_SC1].push_back(-0.569402);

  parameter[CYS_BB].push_back(10.698);
  parameter[CYS_BB].push_back(-0.0233493);
  parameter[CYS_BB].push_back(1.18257);
  parameter[CYS_BB].push_back(0.0684464);
  parameter[CYS_BB].push_back(-2.792);
  parameter[CYS_BB].push_back(1.88995);
  parameter[CYS_BB].push_back(-0.360229);

  parameter[CYS_SC1].push_back(8.199);
  parameter[CYS_SC1].push_back(-0.0261569);
  parameter[CYS_SC1].push_back(6.79677);
  parameter[CYS_SC1].push_back(-0.343845);
  parameter[CYS_SC1].push_back(-5.03578);
  parameter[CYS_SC1].push_back(2.7076);
  parameter[CYS_SC1].push_back(-0.420714);

  parameter[GLN_BB].push_back(10.728);
  parameter[GLN_BB].push_back(-0.0391984);
  parameter[GLN_BB].push_back(1.09264);
  parameter[GLN_BB].push_back(-0.261555);
  parameter[GLN_BB].push_back(-2.21245);
  parameter[GLN_BB].push_back(1.62071);
  parameter[GLN_BB].push_back(-0.322325);

  parameter[GLN_SC1].push_back(8.317);
  parameter[GLN_SC1].push_back(-0.229045);
  parameter[GLN_SC1].push_back(12.6338);
  parameter[GLN_SC1].push_back(-7.6719);
  parameter[GLN_SC1].push_back(-5.8376);
  parameter[GLN_SC1].push_back(5.53784);
  parameter[GLN_SC1].push_back(-1.12604);

  parameter[GLU_BB].push_back(10.694);
  parameter[GLU_BB].push_back(-0.0521961);
  parameter[GLU_BB].push_back(1.11153);
  parameter[GLU_BB].push_back(-0.491995);
  parameter[GLU_BB].push_back(-1.86236);
  parameter[GLU_BB].push_back(1.45332);
  parameter[GLU_BB].push_back(-0.29708);

  parameter[GLU_SC1].push_back(8.544);
  parameter[GLU_SC1].push_back(-0.249555);
  parameter[GLU_SC1].push_back(12.8031);
  parameter[GLU_SC1].push_back(-8.42696);
  parameter[GLU_SC1].push_back(-4.66486);
  parameter[GLU_SC1].push_back(4.90004);
  parameter[GLU_SC1].push_back(-1.01204);

  parameter[GLY_BB].push_back(9.977);
  parameter[GLY_BB].push_back(-0.0285799);
  parameter[GLY_BB].push_back(1.84236);
  parameter[GLY_BB].push_back(-0.0315192);
  parameter[GLY_BB].push_back(-2.88326);
  parameter[GLY_BB].push_back(1.87323);
  parameter[GLY_BB].push_back(-0.345773);

  parameter[HIS_BB].push_back(10.721);
  parameter[HIS_BB].push_back(-0.0379337);
  parameter[HIS_BB].push_back(1.06028);
  parameter[HIS_BB].push_back(-0.236143);
  parameter[HIS_BB].push_back(-2.17819);
  parameter[HIS_BB].push_back(1.58357);
  parameter[HIS_BB].push_back(-0.31345);

  parameter[HIS_SC1].push_back(-0.424);
  parameter[HIS_SC1].push_back(0.665176);
  parameter[HIS_SC1].push_back(3.4369);
  parameter[HIS_SC1].push_back(2.93795);
  parameter[HIS_SC1].push_back(-5.18288);
  parameter[HIS_SC1].push_back(2.12381);
  parameter[HIS_SC1].push_back(-0.284224);

  parameter[HIS_SC2].push_back(5.363);
  parameter[HIS_SC2].push_back(-0.0176945);
  parameter[HIS_SC2].push_back(2.9506);
  parameter[HIS_SC2].push_back(-0.387018);
  parameter[HIS_SC2].push_back(-1.83951);
  parameter[HIS_SC2].push_back(0.9703);
  parameter[HIS_SC2].push_back(-0.1458);

  parameter[HIS_SC3].push_back(5.784);
  parameter[HIS_SC3].push_back(-0.0293129);
  parameter[HIS_SC3].push_back(2.74167);
  parameter[HIS_SC3].push_back(-0.520875);
  parameter[HIS_SC3].push_back(-1.62949);
  parameter[HIS_SC3].push_back(0.902379);
  parameter[HIS_SC3].push_back(-0.139957);

  parameter[ILE_BB].push_back(10.699);
  parameter[ILE_BB].push_back(-0.0188962);
  parameter[ILE_BB].push_back(1.217);
  parameter[ILE_BB].push_back(0.242481);
  parameter[ILE_BB].push_back(-3.13898);
  parameter[ILE_BB].push_back(2.07916);
  parameter[ILE_BB].push_back(-0.392574);

  parameter[ILE_SC1].push_back(-4.448);
  parameter[ILE_SC1].push_back(1.20996);
  parameter[ILE_SC1].push_back(11.5141);
  parameter[ILE_SC1].push_back(6.98895);
  parameter[ILE_SC1].push_back(-19.1948);
  parameter[ILE_SC1].push_back(9.89207);
  parameter[ILE_SC1].push_back(-1.60877);

  parameter[LEU_BB].push_back(10.692);
  parameter[LEU_BB].push_back(-0.0414917);
  parameter[LEU_BB].push_back(1.1077);
  parameter[LEU_BB].push_back(-0.288062);
  parameter[LEU_BB].push_back(-2.17187);
  parameter[LEU_BB].push_back(1.59879);
  parameter[LEU_BB].push_back(-0.318545);

  parameter[LEU_SC1].push_back(-4.448);
  parameter[LEU_SC1].push_back(2.1063);
  parameter[LEU_SC1].push_back(6.72381);
  parameter[LEU_SC1].push_back(14.6954);
  parameter[LEU_SC1].push_back(-23.7197);
  parameter[LEU_SC1].push_back(10.7247);
  parameter[LEU_SC1].push_back(-1.59146);

  parameter[LYS_BB].push_back(10.706);
  parameter[LYS_BB].push_back(-0.0468629);
  parameter[LYS_BB].push_back(1.09477);
  parameter[LYS_BB].push_back(-0.432751);
  parameter[LYS_BB].push_back(-1.94335);
  parameter[LYS_BB].push_back(1.49109);
  parameter[LYS_BB].push_back(-0.302589);

  parameter[LYS_SC1].push_back(-2.796);
  parameter[LYS_SC1].push_back(0.508044);
  parameter[LYS_SC1].push_back(7.91436);
  parameter[LYS_SC1].push_back(4.54097);
  parameter[LYS_SC1].push_back(-10.8051);
  parameter[LYS_SC1].push_back(4.96204);
  parameter[LYS_SC1].push_back(-0.724414);

  parameter[LYS_SC2].push_back(3.070);
  parameter[LYS_SC2].push_back(-0.0101448);
  parameter[LYS_SC2].push_back(4.67994);
  parameter[LYS_SC2].push_back(-0.792529);
  parameter[LYS_SC2].push_back(-2.09142);
  parameter[LYS_SC2].push_back(1.02933);
  parameter[LYS_SC2].push_back(-0.137787);

  parameter[MET_BB].push_back(10.671);
  parameter[MET_BB].push_back(-0.0433724);
  parameter[MET_BB].push_back(1.13784);
  parameter[MET_BB].push_back(-0.40768);
  parameter[MET_BB].push_back(-2.00555);
  parameter[MET_BB].push_back(1.51673);
  parameter[MET_BB].push_back(-0.305547);

  parameter[MET_SC1].push_back(5.85);
  parameter[MET_SC1].push_back(-0.0485798);
  parameter[MET_SC1].push_back(17.0391);
  parameter[MET_SC1].push_back(-3.65327);
  parameter[MET_SC1].push_back(-13.174);
  parameter[MET_SC1].push_back(8.68286);
  parameter[MET_SC1].push_back(-1.56095);

  parameter[PHE_BB].push_back(10.741);
  parameter[PHE_BB].push_back(-0.0317275);
  parameter[PHE_BB].push_back(1.15599);
  parameter[PHE_BB].push_back(0.0276187);
  parameter[PHE_BB].push_back(-2.74757);
  parameter[PHE_BB].push_back(1.88783);
  parameter[PHE_BB].push_back(-0.363525);

  parameter[PHE_SC1].push_back(-0.636);
  parameter[PHE_SC1].push_back(0.527882);
  parameter[PHE_SC1].push_back(6.77612);
  parameter[PHE_SC1].push_back(3.18508);
  parameter[PHE_SC1].push_back(-8.92826);
  parameter[PHE_SC1].push_back(4.29752);
  parameter[PHE_SC1].push_back(-0.65187);

  parameter[PHE_SC2].push_back(-0.424);
  parameter[PHE_SC2].push_back(0.389174);
  parameter[PHE_SC2].push_back(4.11761);
  parameter[PHE_SC2].push_back(2.29527);
  parameter[PHE_SC2].push_back(-4.7652);
  parameter[PHE_SC2].push_back(1.97023);
  parameter[PHE_SC2].push_back(-0.262318);

  parameter[PHE_SC3].push_back(-0.424);
  parameter[PHE_SC3].push_back(0.38927);
  parameter[PHE_SC3].push_back(4.11708);
  parameter[PHE_SC3].push_back(2.29623);
  parameter[PHE_SC3].push_back(-4.76592);
  parameter[PHE_SC3].push_back(1.97055);
  parameter[PHE_SC3].push_back(-0.262381);

  parameter[PRO_BB].push_back(11.434);
  parameter[PRO_BB].push_back(-0.033323);
  parameter[PRO_BB].push_back(0.472014);
  parameter[PRO_BB].push_back(-0.290854);
  parameter[PRO_BB].push_back(-1.81409);
  parameter[PRO_BB].push_back(1.39751);
  parameter[PRO_BB].push_back(-0.280407);

  parameter[PRO_SC1].push_back(-2.796);
  parameter[PRO_SC1].push_back(0.95668);
  parameter[PRO_SC1].push_back(6.84197);
  parameter[PRO_SC1].push_back(6.43774);
  parameter[PRO_SC1].push_back(-12.5068);
  parameter[PRO_SC1].push_back(5.64597);
  parameter[PRO_SC1].push_back(-0.825206);

  parameter[SER_BB].push_back(10.699);
  parameter[SER_BB].push_back(-0.0325828);
  parameter[SER_BB].push_back(1.20329);
  parameter[SER_BB].push_back(-0.0674351);
  parameter[SER_BB].push_back(-2.60749);
  parameter[SER_BB].push_back(1.80318);
  parameter[SER_BB].push_back(-0.346803);

  parameter[SER_SC1].push_back(3.298);
  parameter[SER_SC1].push_back(-0.0366801);
  parameter[SER_SC1].push_back(5.11077);
  parameter[SER_SC1].push_back(-1.46774);
  parameter[SER_SC1].push_back(-1.48421);
  parameter[SER_SC1].push_back(0.800326);
  parameter[SER_SC1].push_back(-0.108314);

  parameter[THR_BB].push_back(10.697);
  parameter[THR_BB].push_back(-0.0242955);
  parameter[THR_BB].push_back(1.24671);
  parameter[THR_BB].push_back(0.146423);
  parameter[THR_BB].push_back(-2.97429);
  parameter[THR_BB].push_back(1.97513);
  parameter[THR_BB].push_back(-0.371479);

  parameter[THR_SC1].push_back(2.366);
  parameter[THR_SC1].push_back(0.0297604);
  parameter[THR_SC1].push_back(11.9216);
  parameter[THR_SC1].push_back(-9.32503);
  parameter[THR_SC1].push_back(1.9396);
  parameter[THR_SC1].push_back(0.0804861);
  parameter[THR_SC1].push_back(-0.0302721);

  parameter[TRP_BB].push_back(10.689);
  parameter[TRP_BB].push_back(-0.0265879);
  parameter[TRP_BB].push_back(1.17819);
  parameter[TRP_BB].push_back(0.0386457);
  parameter[TRP_BB].push_back(-2.75634);
  parameter[TRP_BB].push_back(1.88065);
  parameter[TRP_BB].push_back(-0.360217);

  parameter[TRP_SC1].push_back(0.084);
  parameter[TRP_SC1].push_back(0.752407);
  parameter[TRP_SC1].push_back(5.3802);
  parameter[TRP_SC1].push_back(4.09281);
  parameter[TRP_SC1].push_back(-9.28029);
  parameter[TRP_SC1].push_back(4.45923);
  parameter[TRP_SC1].push_back(-0.689008);

  parameter[TRP_SC2].push_back(5.739);
  parameter[TRP_SC2].push_back(0.0298492);
  parameter[TRP_SC2].push_back(4.60446);
  parameter[TRP_SC2].push_back(1.34463);
  parameter[TRP_SC2].push_back(-5.69968);
  parameter[TRP_SC2].push_back(2.84924);
  parameter[TRP_SC2].push_back(-0.433781);

  parameter[TRP_SC3].push_back(-0.424);
  parameter[TRP_SC3].push_back(0.388576);
  parameter[TRP_SC3].push_back(4.11859);
  parameter[TRP_SC3].push_back(2.29485);
  parameter[TRP_SC3].push_back(-4.76255);
  parameter[TRP_SC3].push_back(1.96849);
  parameter[TRP_SC3].push_back(-0.262015);

  parameter[TRP_SC4].push_back(-0.424);
  parameter[TRP_SC4].push_back(0.387685);
  parameter[TRP_SC4].push_back(4.12153);
  parameter[TRP_SC4].push_back(2.29144);
  parameter[TRP_SC4].push_back(-4.7589);
  parameter[TRP_SC4].push_back(1.96686);
  parameter[TRP_SC4].push_back(-0.261786);

  parameter[TYR_BB].push_back(10.689);
  parameter[TYR_BB].push_back(-0.0193526);
  parameter[TYR_BB].push_back(1.18241);
  parameter[TYR_BB].push_back(0.207318);
  parameter[TYR_BB].push_back(-3.0041);
  parameter[TYR_BB].push_back(1.99335);
  parameter[TYR_BB].push_back(-0.376482);

  parameter[TYR_SC1].push_back(-0.636);
  parameter[TYR_SC1].push_back(0.528902);
  parameter[TYR_SC1].push_back(6.78168);
  parameter[TYR_SC1].push_back(3.17769);
  parameter[TYR_SC1].push_back(-8.93667);
  parameter[TYR_SC1].push_back(4.30692);
  parameter[TYR_SC1].push_back(-0.653993);

  parameter[TYR_SC2].push_back(-0.424);
  parameter[TYR_SC2].push_back(0.388811);
  parameter[TYR_SC2].push_back(4.11851);
  parameter[TYR_SC2].push_back(2.29545);
  parameter[TYR_SC2].push_back(-4.7668);
  parameter[TYR_SC2].push_back(1.97131);
  parameter[TYR_SC2].push_back(-0.262534);

  parameter[TYR_SC3].push_back(4.526);
  parameter[TYR_SC3].push_back(-0.00381305);
  parameter[TYR_SC3].push_back(5.8567);
  parameter[TYR_SC3].push_back(-0.214086);
  parameter[TYR_SC3].push_back(-4.63649);
  parameter[TYR_SC3].push_back(2.52869);
  parameter[TYR_SC3].push_back(-0.39894);

  parameter[VAL_BB].push_back(10.691);
  parameter[VAL_BB].push_back(-0.0162929);
  parameter[VAL_BB].push_back(1.24446);
  parameter[VAL_BB].push_back(0.307914);
  parameter[VAL_BB].push_back(-3.27446);
  parameter[VAL_BB].push_back(2.14788);
  parameter[VAL_BB].push_back(-0.403259);

  parameter[VAL_SC1].push_back(-3.516);
  parameter[VAL_SC1].push_back(1.62307);
  parameter[VAL_SC1].push_back(5.43064);
  parameter[VAL_SC1].push_back(9.28809);
  parameter[VAL_SC1].push_back(-14.9927);
  parameter[VAL_SC1].push_back(6.6133);
  parameter[VAL_SC1].push_back(-0.964977);

  parameter[A_BB1].push_back(32.88500000);
  parameter[A_BB1].push_back(0.08339900);
  parameter[A_BB1].push_back(-7.36054400);
  parameter[A_BB1].push_back(2.19220300);
  parameter[A_BB1].push_back(-3.56523400);
  parameter[A_BB1].push_back(2.33326900);
  parameter[A_BB1].push_back(-0.39785500);

  parameter[A_BB2].push_back(3.80600000);
  parameter[A_BB2].push_back(-0.10727600);
  parameter[A_BB2].push_back(9.58854100);
  parameter[A_BB2].push_back(-6.23740500);
  parameter[A_BB2].push_back(-0.48267300);
  parameter[A_BB2].push_back(1.14119500);
  parameter[A_BB2].push_back(-0.21385600);

  parameter[A_BB3].push_back(3.59400000);
  parameter[A_BB3].push_back(0.04537300);
  parameter[A_BB3].push_back(9.59178900);
  parameter[A_BB3].push_back(-1.29202200);
  parameter[A_BB3].push_back(-7.10851000);
  parameter[A_BB3].push_back(4.05571200);
  parameter[A_BB3].push_back(-0.63372500);

  parameter[A_SC1].push_back(6.67100000);
  parameter[A_SC1].push_back(-0.00855300);
  parameter[A_SC1].push_back(1.63222400);
  parameter[A_SC1].push_back(-0.06466200);
  parameter[A_SC1].push_back(-1.48694200);
  parameter[A_SC1].push_back(0.78544600);
  parameter[A_SC1].push_back(-0.12083500);

  parameter[A_SC2].push_back(5.95100000);
  parameter[A_SC2].push_back(-0.02606600);
  parameter[A_SC2].push_back(2.54399900);
  parameter[A_SC2].push_back(-0.48436900);
  parameter[A_SC2].push_back(-1.55357400);
  parameter[A_SC2].push_back(0.86466900);
  parameter[A_SC2].push_back(-0.13509000);

  parameter[A_SC3].push_back(11.39400000);
  parameter[A_SC3].push_back(0.00871300);
  parameter[A_SC3].push_back(-0.23891300);
  parameter[A_SC3].push_back(0.48919400);
  parameter[A_SC3].push_back(-1.75289400);
  parameter[A_SC3].push_back(0.99267500);
  parameter[A_SC3].push_back(-0.16291300);

  parameter[A_SC4].push_back(6.45900000);
  parameter[A_SC4].push_back(0.01990600);
  parameter[A_SC4].push_back(4.17970400);
  parameter[A_SC4].push_back(0.97629900);
  parameter[A_SC4].push_back(-5.03297800);
  parameter[A_SC4].push_back(2.55576700);
  parameter[A_SC4].push_back(-0.39150500);

  parameter[A_3TE].push_back(4.23000000);
  parameter[A_3TE].push_back(0.00064800);
  parameter[A_3TE].push_back(0.92124600);
  parameter[A_3TE].push_back(0.08064300);
  parameter[A_3TE].push_back(-0.39054400);
  parameter[A_3TE].push_back(0.12429100);
  parameter[A_3TE].push_back(-0.01122700);

  parameter[A_5TE].push_back(4.23000000);
  parameter[A_5TE].push_back(0.00039300);
  parameter[A_5TE].push_back(0.92305100);
  parameter[A_5TE].push_back(0.07747500);
  parameter[A_5TE].push_back(-0.38792100);
  parameter[A_5TE].push_back(0.12323800);
  parameter[A_5TE].push_back(-0.01106600);

  parameter[A_TE3].push_back(7.82400000);
  parameter[A_TE3].push_back(-0.04881000);
  parameter[A_TE3].push_back(8.21557900);
  parameter[A_TE3].push_back(-0.89491400);
  parameter[A_TE3].push_back(-9.54293700);
  parameter[A_TE3].push_back(6.33122200);
  parameter[A_TE3].push_back(-1.16672900);

  parameter[A_TE5].push_back(8.03600000);
  parameter[A_TE5].push_back(0.01641200);
  parameter[A_TE5].push_back(5.14902200);
  parameter[A_TE5].push_back(0.83419700);
  parameter[A_TE5].push_back(-7.59068300);
  parameter[A_TE5].push_back(4.52063200);
  parameter[A_TE5].push_back(-0.78260800);

  parameter[C_BB1].push_back(32.88500000);
  parameter[C_BB1].push_back(0.08311100);
  parameter[C_BB1].push_back(-7.35432100);
  parameter[C_BB1].push_back(2.18610000);
  parameter[C_BB1].push_back(-3.55788300);
  parameter[C_BB1].push_back(2.32918700);
  parameter[C_BB1].push_back(-0.39720000);

  parameter[C_BB2].push_back(3.80600000);
  parameter[C_BB2].push_back(-0.10808100);
  parameter[C_BB2].push_back(9.61612600);
  parameter[C_BB2].push_back(-6.28595400);
  parameter[C_BB2].push_back(-0.45187000);
  parameter[C_BB2].push_back(1.13326000);
  parameter[C_BB2].push_back(-0.21320300);

  parameter[C_BB3].push_back(3.59400000);
  parameter[C_BB3].push_back(0.04484200);
  parameter[C_BB3].push_back(9.61919800);
  parameter[C_BB3].push_back(-1.33582800);
  parameter[C_BB3].push_back(-7.07200400);
  parameter[C_BB3].push_back(4.03952900);
  parameter[C_BB3].push_back(-0.63098200);

  parameter[C_SC1].push_back(5.95100000);
  parameter[C_SC1].push_back(-0.02911300);
  parameter[C_SC1].push_back(2.59700400);
  parameter[C_SC1].push_back(-0.55507700);
  parameter[C_SC1].push_back(-1.56344600);
  parameter[C_SC1].push_back(0.88956200);
  parameter[C_SC1].push_back(-0.14061300);

  parameter[C_SC2].push_back(11.62100000);
  parameter[C_SC2].push_back(0.01366100);
  parameter[C_SC2].push_back(-0.25959200);
  parameter[C_SC2].push_back(0.48918300);
  parameter[C_SC2].push_back(-1.52550500);
  parameter[C_SC2].push_back(0.83644100);
  parameter[C_SC2].push_back(-0.13407300);

  parameter[C_SC3].push_back(5.01900000);
  parameter[C_SC3].push_back(-0.03276100);
  parameter[C_SC3].push_back(5.53776900);
  parameter[C_SC3].push_back(-0.95105000);
  parameter[C_SC3].push_back(-3.71130800);
  parameter[C_SC3].push_back(2.16146000);
  parameter[C_SC3].push_back(-0.34918600);

  parameter[C_3TE].push_back(4.23000000);
  parameter[C_3TE].push_back(0.00057300);
  parameter[C_3TE].push_back(0.92174800);
  parameter[C_3TE].push_back(0.07964500);
  parameter[C_3TE].push_back(-0.38965700);
  parameter[C_3TE].push_back(0.12392500);
  parameter[C_3TE].push_back(-0.01117000);

  parameter[C_5TE].push_back(4.23000000);
  parameter[C_5TE].push_back(0.00071000);
  parameter[C_5TE].push_back(0.92082800);
  parameter[C_5TE].push_back(0.08150600);
  parameter[C_5TE].push_back(-0.39127000);
  parameter[C_5TE].push_back(0.12455900);
  parameter[C_5TE].push_back(-0.01126300);

  parameter[C_TE3].push_back(7.82400000);
  parameter[C_TE3].push_back(-0.05848300);
  parameter[C_TE3].push_back(8.29319900);
  parameter[C_TE3].push_back(-1.12563800);
  parameter[C_TE3].push_back(-9.42197600);
  parameter[C_TE3].push_back(6.35441700);
  parameter[C_TE3].push_back(-1.18356900);

  parameter[C_TE5].push_back(8.03600000);
  parameter[C_TE5].push_back(0.00493500);
  parameter[C_TE5].push_back(4.92622000);
  parameter[C_TE5].push_back(0.64810700);
  parameter[C_TE5].push_back(-7.05100000);
  parameter[C_TE5].push_back(4.26064400);
  parameter[C_TE5].push_back(-0.74819100);

  parameter[G_BB1].push_back(32.88500000);
  parameter[G_BB1].push_back(0.08325400);
  parameter[G_BB1].push_back(-7.35736000);
  parameter[G_BB1].push_back(2.18914800);
  parameter[G_BB1].push_back(-3.56154800);
  parameter[G_BB1].push_back(2.33120600);
  parameter[G_BB1].push_back(-0.39752300);

  parameter[G_BB2].push_back(3.80600000);
  parameter[G_BB2].push_back(-0.10788300);
  parameter[G_BB2].push_back(9.60930800);
  parameter[G_BB2].push_back(-6.27402500);
  parameter[G_BB2].push_back(-0.46192700);
  parameter[G_BB2].push_back(1.13737000);
  parameter[G_BB2].push_back(-0.21383100);

  parameter[G_BB3].push_back(3.59400000);
  parameter[G_BB3].push_back(0.04514500);
  parameter[G_BB3].push_back(9.61234700);
  parameter[G_BB3].push_back(-1.31542100);
  parameter[G_BB3].push_back(-7.09150500);
  parameter[G_BB3].push_back(4.04706200);
  parameter[G_BB3].push_back(-0.63201000);

  parameter[G_SC1].push_back(6.67100000);
  parameter[G_SC1].push_back(-0.00863200);
  parameter[G_SC1].push_back(1.63252300);
  parameter[G_SC1].push_back(-0.06567200);
  parameter[G_SC1].push_back(-1.48680500);
  parameter[G_SC1].push_back(0.78565600);
  parameter[G_SC1].push_back(-0.12088900);

  parameter[G_SC2].push_back(11.39400000);
  parameter[G_SC2].push_back(0.00912200);
  parameter[G_SC2].push_back(-0.22869000);
  parameter[G_SC2].push_back(0.49616400);
  parameter[G_SC2].push_back(-1.75039000);
  parameter[G_SC2].push_back(0.98649200);
  parameter[G_SC2].push_back(-0.16141600);

  parameter[G_SC3].push_back(10.90100000);
  parameter[G_SC3].push_back(0.02208700);
  parameter[G_SC3].push_back(0.17032800);
  parameter[G_SC3].push_back(0.73280800);
  parameter[G_SC3].push_back(-1.95292000);
  parameter[G_SC3].push_back(0.98357600);
  parameter[G_SC3].push_back(-0.14790900);

  parameter[G_SC4].push_back(6.45900000);
  parameter[G_SC4].push_back(0.02023700);
  parameter[G_SC4].push_back(4.17655400);
  parameter[G_SC4].push_back(0.98731800);
  parameter[G_SC4].push_back(-5.04352800);
  parameter[G_SC4].push_back(2.56059400);
  parameter[G_SC4].push_back(-0.39234300);

  parameter[G_3TE].push_back(4.23000000);
  parameter[G_3TE].push_back(0.00066300);
  parameter[G_3TE].push_back(0.92118800);
  parameter[G_3TE].push_back(0.08062700);
  parameter[G_3TE].push_back(-0.39041600);
  parameter[G_3TE].push_back(0.12419400);
  parameter[G_3TE].push_back(-0.01120500);

  parameter[G_5TE].push_back(4.23000000);
  parameter[G_5TE].push_back(0.00062800);
  parameter[G_5TE].push_back(0.92133500);
  parameter[G_5TE].push_back(0.08029900);
  parameter[G_5TE].push_back(-0.39015300);
  parameter[G_5TE].push_back(0.12411600);
  parameter[G_5TE].push_back(-0.01119900);

  parameter[G_TE3].push_back(7.82400000);
  parameter[G_TE3].push_back(-0.05177400);
  parameter[G_TE3].push_back(8.34606700);
  parameter[G_TE3].push_back(-1.02936300);
  parameter[G_TE3].push_back(-9.55211900);
  parameter[G_TE3].push_back(6.37776600);
  parameter[G_TE3].push_back(-1.17898000);

  parameter[G_TE5].push_back(8.03600000);
  parameter[G_TE5].push_back(0.00525100);
  parameter[G_TE5].push_back(4.71070600);
  parameter[G_TE5].push_back(0.66746900);
  parameter[G_TE5].push_back(-6.72538700);
  parameter[G_TE5].push_back(4.03644100);
  parameter[G_TE5].push_back(-0.70605700);

  parameter[U_BB1].push_back(32.88500000);
  parameter[U_BB1].push_back(0.08321400);
  parameter[U_BB1].push_back(-7.35634900);
  parameter[U_BB1].push_back(2.18826800);
  parameter[U_BB1].push_back(-3.56047400);
  parameter[U_BB1].push_back(2.33064700);
  parameter[U_BB1].push_back(-0.39744000);

  parameter[U_BB2].push_back(3.80600000);
  parameter[U_BB2].push_back(-0.10773100);
  parameter[U_BB2].push_back(9.60099900);
  parameter[U_BB2].push_back(-6.26131900);
  parameter[U_BB2].push_back(-0.46668300);
  parameter[U_BB2].push_back(1.13698100);
  parameter[U_BB2].push_back(-0.21351600);

  parameter[U_BB3].push_back(3.59400000);
  parameter[U_BB3].push_back(0.04544300);
  parameter[U_BB3].push_back(9.59625900);
  parameter[U_BB3].push_back(-1.29222200);
  parameter[U_BB3].push_back(-7.11143200);
  parameter[U_BB3].push_back(4.05687700);
  parameter[U_BB3].push_back(-0.63382800);

  parameter[U_SC1].push_back(5.95100000);
  parameter[U_SC1].push_back(-0.02924500);
  parameter[U_SC1].push_back(2.59668700);
  parameter[U_SC1].push_back(-0.56118700);
  parameter[U_SC1].push_back(-1.56477100);
  parameter[U_SC1].push_back(0.89265100);
  parameter[U_SC1].push_back(-0.14130800);

  parameter[U_SC2].push_back(10.90100000);
  parameter[U_SC2].push_back(0.02178900);
  parameter[U_SC2].push_back(0.18839000);
  parameter[U_SC2].push_back(0.72223100);
  parameter[U_SC2].push_back(-1.92581600);
  parameter[U_SC2].push_back(0.96654300);
  parameter[U_SC2].push_back(-0.14501300);

  parameter[U_SC3].push_back(5.24600000);
  parameter[U_SC3].push_back(-0.04586500);
  parameter[U_SC3].push_back(5.89978100);
  parameter[U_SC3].push_back(-1.50664700);
  parameter[U_SC3].push_back(-3.17054400);
  parameter[U_SC3].push_back(1.93717100);
  parameter[U_SC3].push_back(-0.31701000);

  parameter[U_3TE].push_back(4.23000000);
  parameter[U_3TE].push_back(0.00067500);
  parameter[U_3TE].push_back(0.92102300);
  parameter[U_3TE].push_back(0.08100800);
  parameter[U_3TE].push_back(-0.39084300);
  parameter[U_3TE].push_back(0.12441900);
  parameter[U_3TE].push_back(-0.01124900);

  parameter[U_5TE].push_back(4.23000000);
  parameter[U_5TE].push_back(0.00059000);
  parameter[U_5TE].push_back(0.92154600);
  parameter[U_5TE].push_back(0.07968200);
  parameter[U_5TE].push_back(-0.38950100);
  parameter[U_5TE].push_back(0.12382500);
  parameter[U_5TE].push_back(-0.01115100);

  parameter[U_TE3].push_back(7.82400000);
  parameter[U_TE3].push_back(-0.02968100);
  parameter[U_TE3].push_back(7.93783200);
  parameter[U_TE3].push_back(-0.33078100);
  parameter[U_TE3].push_back(-10.14120200);
  parameter[U_TE3].push_back(6.63334700);
  parameter[U_TE3].push_back(-1.22111200);

  parameter[U_TE5].push_back(8.03600000);
  parameter[U_TE5].push_back(-0.00909700);
  parameter[U_TE5].push_back(4.33193500);
  parameter[U_TE5].push_back(0.43416500);
  parameter[U_TE5].push_back(-5.80831400);
  parameter[U_TE5].push_back(3.52438800);
  parameter[U_TE5].push_back(-0.62382400);

  parameter[DA_BB1].push_back(32.88500000);
  parameter[DA_BB1].push_back(0.08179900);
  parameter[DA_BB1].push_back(-7.31735900);
  parameter[DA_BB1].push_back(2.15614500);
  parameter[DA_BB1].push_back(-3.52263200);
  parameter[DA_BB1].push_back(2.30604700);
  parameter[DA_BB1].push_back(-0.39270100);

  parameter[DA_BB2].push_back(3.80600000);
  parameter[DA_BB2].push_back(-0.10597700);
  parameter[DA_BB2].push_back(9.52537500);
  parameter[DA_BB2].push_back(-6.12991000);
  parameter[DA_BB2].push_back(-0.54092600);
  parameter[DA_BB2].push_back(1.15429100);
  parameter[DA_BB2].push_back(-0.21503500);

  parameter[DA_BB3].push_back(-1.35600000);
  parameter[DA_BB3].push_back(0.58928300);
  parameter[DA_BB3].push_back(6.71894100);
  parameter[DA_BB3].push_back(4.14050900);
  parameter[DA_BB3].push_back(-9.65859900);
  parameter[DA_BB3].push_back(4.43185000);
  parameter[DA_BB3].push_back(-0.64657300);

  parameter[DA_SC1].push_back(6.67100000);
  parameter[DA_SC1].push_back(-0.00871400);
  parameter[DA_SC1].push_back(1.63289100);
  parameter[DA_SC1].push_back(-0.06637700);
  parameter[DA_SC1].push_back(-1.48632900);
  parameter[DA_SC1].push_back(0.78551800);
  parameter[DA_SC1].push_back(-0.12087300);

  parameter[DA_SC2].push_back(5.95100000);
  parameter[DA_SC2].push_back(-0.02634300);
  parameter[DA_SC2].push_back(2.54864300);
  parameter[DA_SC2].push_back(-0.49015800);
  parameter[DA_SC2].push_back(-1.55386900);
  parameter[DA_SC2].push_back(0.86630200);
  parameter[DA_SC2].push_back(-0.13546200);

  parameter[DA_SC3].push_back(11.39400000);
  parameter[DA_SC3].push_back(0.00859500);
  parameter[DA_SC3].push_back(-0.25471400);
  parameter[DA_SC3].push_back(0.48718800);
  parameter[DA_SC3].push_back(-1.74520000);
  parameter[DA_SC3].push_back(0.99246200);
  parameter[DA_SC3].push_back(-0.16351900);

  parameter[DA_SC4].push_back(6.45900000);
  parameter[DA_SC4].push_back(0.01991800);
  parameter[DA_SC4].push_back(4.17962300);
  parameter[DA_SC4].push_back(0.97469100);
  parameter[DA_SC4].push_back(-5.02950400);
  parameter[DA_SC4].push_back(2.55371800);
  parameter[DA_SC4].push_back(-0.39113400);

  parameter[DA_3TE].push_back(4.23000000);
  parameter[DA_3TE].push_back(0.00062600);
  parameter[DA_3TE].push_back(0.92142000);
  parameter[DA_3TE].push_back(0.08016400);
  parameter[DA_3TE].push_back(-0.39000300);
  parameter[DA_3TE].push_back(0.12402500);
  parameter[DA_3TE].push_back(-0.01117900);

  parameter[DA_5TE].push_back(4.23000000);
  parameter[DA_5TE].push_back(0.00055500);
  parameter[DA_5TE].push_back(0.92183900);
  parameter[DA_5TE].push_back(0.07907600);
  parameter[DA_5TE].push_back(-0.38895100);
  parameter[DA_5TE].push_back(0.12359600);
  parameter[DA_5TE].push_back(-0.01111600);

  parameter[DA_TE3].push_back(2.87400000);
  parameter[DA_TE3].push_back(0.00112900);
  parameter[DA_TE3].push_back(12.51167200);
  parameter[DA_TE3].push_back(-7.67548000);
  parameter[DA_TE3].push_back(-2.02234000);
  parameter[DA_TE3].push_back(2.50837100);
  parameter[DA_TE3].push_back(-0.49458500);

  parameter[DA_TE5].push_back(8.03600000);
  parameter[DA_TE5].push_back(0.00473100);
  parameter[DA_TE5].push_back(4.65554400);
  parameter[DA_TE5].push_back(0.66424100);
  parameter[DA_TE5].push_back(-6.62131300);
  parameter[DA_TE5].push_back(3.96107400);
  parameter[DA_TE5].push_back(-0.69075800);

  parameter[DC_BB1].push_back(32.88500000);
  parameter[DC_BB1].push_back(0.08189900);
  parameter[DC_BB1].push_back(-7.32493500);
  parameter[DC_BB1].push_back(2.15976900);
  parameter[DC_BB1].push_back(-3.52612100);
  parameter[DC_BB1].push_back(2.31058600);
  parameter[DC_BB1].push_back(-0.39402700);

  parameter[DC_BB2].push_back(3.80600000);
  parameter[DC_BB2].push_back(-0.10559800);
  parameter[DC_BB2].push_back(9.52527700);
  parameter[DC_BB2].push_back(-6.12131700);
  parameter[DC_BB2].push_back(-0.54899400);
  parameter[DC_BB2].push_back(1.15592900);
  parameter[DC_BB2].push_back(-0.21494500);

  parameter[DC_BB3].push_back(-1.35600000);
  parameter[DC_BB3].push_back(0.55525700);
  parameter[DC_BB3].push_back(6.80305500);
  parameter[DC_BB3].push_back(4.05924700);
  parameter[DC_BB3].push_back(-9.61034700);
  parameter[DC_BB3].push_back(4.41253800);
  parameter[DC_BB3].push_back(-0.64315100);

  parameter[DC_SC1].push_back(5.95100000);
  parameter[DC_SC1].push_back(-0.02899900);
  parameter[DC_SC1].push_back(2.59587800);
  parameter[DC_SC1].push_back(-0.55388300);
  parameter[DC_SC1].push_back(-1.56395100);
  parameter[DC_SC1].push_back(0.88967400);
  parameter[DC_SC1].push_back(-0.14062500);

  parameter[DC_SC2].push_back(11.62100000);
  parameter[DC_SC2].push_back(0.01358100);
  parameter[DC_SC2].push_back(-0.24913000);
  parameter[DC_SC2].push_back(0.48787200);
  parameter[DC_SC2].push_back(-1.52867300);
  parameter[DC_SC2].push_back(0.83694900);
  parameter[DC_SC2].push_back(-0.13395300);

  parameter[DC_SC3].push_back(5.01900000);
  parameter[DC_SC3].push_back(-0.03298400);
  parameter[DC_SC3].push_back(5.54242800);
  parameter[DC_SC3].push_back(-0.96081500);
  parameter[DC_SC3].push_back(-3.71051600);
  parameter[DC_SC3].push_back(2.16500200);
  parameter[DC_SC3].push_back(-0.35023400);

  parameter[DC_3TE].push_back(4.23000000);
  parameter[DC_3TE].push_back(0.00055700);
  parameter[DC_3TE].push_back(0.92181400);
  parameter[DC_3TE].push_back(0.07924000);
  parameter[DC_3TE].push_back(-0.38916400);
  parameter[DC_3TE].push_back(0.12369900);
  parameter[DC_3TE].push_back(-0.01113300);

  parameter[DC_5TE].push_back(4.23000000);
  parameter[DC_5TE].push_back(0.00066500);
  parameter[DC_5TE].push_back(0.92103900);
  parameter[DC_5TE].push_back(0.08064600);
  parameter[DC_5TE].push_back(-0.39034900);
  parameter[DC_5TE].push_back(0.12417600);
  parameter[DC_5TE].push_back(-0.01120600);

  parameter[DC_TE3].push_back(2.87400000);
  parameter[DC_TE3].push_back(-0.05235500);
  parameter[DC_TE3].push_back(13.09201200);
  parameter[DC_TE3].push_back(-9.48128200);
  parameter[DC_TE3].push_back(-0.14958600);
  parameter[DC_TE3].push_back(1.75537200);
  parameter[DC_TE3].push_back(-0.39347500);

  parameter[DC_TE5].push_back(8.03600000);
  parameter[DC_TE5].push_back(-0.00513600);
  parameter[DC_TE5].push_back(4.67705700);
  parameter[DC_TE5].push_back(0.48333300);
  parameter[DC_TE5].push_back(-6.34511000);
  parameter[DC_TE5].push_back(3.83388500);
  parameter[DC_TE5].push_back(-0.67367800);

  parameter[DG_BB1].push_back(32.88500000);
  parameter[DG_BB1].push_back(0.08182900);
  parameter[DG_BB1].push_back(-7.32133900);
  parameter[DG_BB1].push_back(2.15767900);
  parameter[DG_BB1].push_back(-3.52369700);
  parameter[DG_BB1].push_back(2.30839600);
  parameter[DG_BB1].push_back(-0.39348300);

  parameter[DG_BB2].push_back(3.80600000);
  parameter[DG_BB2].push_back(-0.10618100);
  parameter[DG_BB2].push_back(9.54169000);
  parameter[DG_BB2].push_back(-6.15177600);
  parameter[DG_BB2].push_back(-0.53462400);
  parameter[DG_BB2].push_back(1.15581300);
  parameter[DG_BB2].push_back(-0.21567000);

  parameter[DG_BB3].push_back(-1.35600000);
  parameter[DG_BB3].push_back(0.57489100);
  parameter[DG_BB3].push_back(6.75164700);
  parameter[DG_BB3].push_back(4.11300900);
  parameter[DG_BB3].push_back(-9.63394600);
  parameter[DG_BB3].push_back(4.41675400);
  parameter[DG_BB3].push_back(-0.64339900);

  parameter[DG_SC1].push_back(6.67100000);
  parameter[DG_SC1].push_back(-0.00886600);
  parameter[DG_SC1].push_back(1.63333000);
  parameter[DG_SC1].push_back(-0.06892100);
  parameter[DG_SC1].push_back(-1.48683500);
  parameter[DG_SC1].push_back(0.78670800);
  parameter[DG_SC1].push_back(-0.12113900);

  parameter[DG_SC2].push_back(11.39400000);
  parameter[DG_SC2].push_back(0.00907900);
  parameter[DG_SC2].push_back(-0.22475500);
  parameter[DG_SC2].push_back(0.49535100);
  parameter[DG_SC2].push_back(-1.75324900);
  parameter[DG_SC2].push_back(0.98767400);
  parameter[DG_SC2].push_back(-0.16150800);

  parameter[DG_SC3].push_back(10.90100000);
  parameter[DG_SC3].push_back(0.02207600);
  parameter[DG_SC3].push_back(0.17932200);
  parameter[DG_SC3].push_back(0.73253200);
  parameter[DG_SC3].push_back(-1.95554900);
  parameter[DG_SC3].push_back(0.98339900);
  parameter[DG_SC3].push_back(-0.14763600);

  parameter[DG_SC4].push_back(6.45900000);
  parameter[DG_SC4].push_back(0.02018400);
  parameter[DG_SC4].push_back(4.17705400);
  parameter[DG_SC4].push_back(0.98531700);
  parameter[DG_SC4].push_back(-5.04354900);
  parameter[DG_SC4].push_back(2.56123700);
  parameter[DG_SC4].push_back(-0.39249300);

  parameter[DG_3TE].push_back(4.23000000);
  parameter[DG_3TE].push_back(0.00061700);
  parameter[DG_3TE].push_back(0.92140100);
  parameter[DG_3TE].push_back(0.08016400);
  parameter[DG_3TE].push_back(-0.39003500);
  parameter[DG_3TE].push_back(0.12406900);
  parameter[DG_3TE].push_back(-0.01119200);

  parameter[DG_5TE].push_back(4.23000000);
  parameter[DG_5TE].push_back(0.00064900);
  parameter[DG_5TE].push_back(0.92110500);
  parameter[DG_5TE].push_back(0.08031500);
  parameter[DG_5TE].push_back(-0.38997000);
  parameter[DG_5TE].push_back(0.12401200);
  parameter[DG_5TE].push_back(-0.01118100);

  parameter[DG_TE3].push_back(2.87400000);
  parameter[DG_TE3].push_back(0.00182000);
  parameter[DG_TE3].push_back(12.41507000);
  parameter[DG_TE3].push_back(-7.47384800);
  parameter[DG_TE3].push_back(-2.11864700);
  parameter[DG_TE3].push_back(2.50112600);
  parameter[DG_TE3].push_back(-0.48652200);

  parameter[DG_TE5].push_back(8.03600000);
  parameter[DG_TE5].push_back(0.00676400);
  parameter[DG_TE5].push_back(4.65989200);
  parameter[DG_TE5].push_back(0.78482500);
  parameter[DG_TE5].push_back(-6.86460600);
  parameter[DG_TE5].push_back(4.11675400);
  parameter[DG_TE5].push_back(-0.72249100);

  parameter[DT_BB1].push_back(32.88500000);
  parameter[DT_BB1].push_back(0.08220100);
  parameter[DT_BB1].push_back(-7.33006800);
  parameter[DT_BB1].push_back(2.16636500);
  parameter[DT_BB1].push_back(-3.53465700);
  parameter[DT_BB1].push_back(2.31447600);
  parameter[DT_BB1].push_back(-0.39445400);

  parameter[DT_BB2].push_back(3.80600000);
  parameter[DT_BB2].push_back(-0.10723000);
  parameter[DT_BB2].push_back(9.56675000);
  parameter[DT_BB2].push_back(-6.20236100);
  parameter[DT_BB2].push_back(-0.49550400);
  parameter[DT_BB2].push_back(1.14300600);
  parameter[DT_BB2].push_back(-0.21420000);

  parameter[DT_BB3].push_back(-1.35600000);
  parameter[DT_BB3].push_back(0.56737900);
  parameter[DT_BB3].push_back(6.76595400);
  parameter[DT_BB3].push_back(4.08976100);
  parameter[DT_BB3].push_back(-9.61512500);
  parameter[DT_BB3].push_back(4.40975100);
  parameter[DT_BB3].push_back(-0.64239800);

  parameter[DT_SC1].push_back(5.95100000);
  parameter[DT_SC1].push_back(-0.02926500);
  parameter[DT_SC1].push_back(2.59630300);
  parameter[DT_SC1].push_back(-0.56152200);
  parameter[DT_SC1].push_back(-1.56532600);
  parameter[DT_SC1].push_back(0.89322800);
  parameter[DT_SC1].push_back(-0.14142900);

  parameter[DT_SC2].push_back(10.90100000);
  parameter[DT_SC2].push_back(0.02183400);
  parameter[DT_SC2].push_back(0.19463000);
  parameter[DT_SC2].push_back(0.72393000);
  parameter[DT_SC2].push_back(-1.93199500);
  parameter[DT_SC2].push_back(0.96856300);
  parameter[DT_SC2].push_back(-0.14512600);

  parameter[DT_SC3].push_back(4.31400000);
  parameter[DT_SC3].push_back(-0.07745600);
  parameter[DT_SC3].push_back(12.49820300);
  parameter[DT_SC3].push_back(-7.64994200);
  parameter[DT_SC3].push_back(-3.00359600);
  parameter[DT_SC3].push_back(3.26263300);
  parameter[DT_SC3].push_back(-0.64498600);

  parameter[DT_3TE].push_back(4.23000000);
  parameter[DT_3TE].push_back(0.00062000);
  parameter[DT_3TE].push_back(0.92141100);
  parameter[DT_3TE].push_back(0.08030900);
  parameter[DT_3TE].push_back(-0.39021500);
  parameter[DT_3TE].push_back(0.12414000);
  parameter[DT_3TE].push_back(-0.01120100);

  parameter[DT_5TE].push_back(4.23000000);
  parameter[DT_5TE].push_back(0.00063700);
  parameter[DT_5TE].push_back(0.92130800);
  parameter[DT_5TE].push_back(0.08026900);
  parameter[DT_5TE].push_back(-0.39007500);
  parameter[DT_5TE].push_back(0.12406600);
  parameter[DT_5TE].push_back(-0.01118800);

  parameter[DT_TE3].push_back(2.87400000);
  parameter[DT_TE3].push_back(-0.00251200);
  parameter[DT_TE3].push_back(12.43576400);
  parameter[DT_TE3].push_back(-7.55343800);
  parameter[DT_TE3].push_back(-2.07363500);
  parameter[DT_TE3].push_back(2.51279300);
  parameter[DT_TE3].push_back(-0.49437100);

  parameter[DT_TE5].push_back(8.03600000);
  parameter[DT_TE5].push_back(0.00119900);
  parameter[DT_TE5].push_back(4.91762300);
  parameter[DT_TE5].push_back(0.65637000);
  parameter[DT_TE5].push_back(-7.23392500);
  parameter[DT_TE5].push_back(4.44636600);
  parameter[DT_TE5].push_back(-0.79467800);

  vector<SetupMolInfo*> moldat=plumed.getActionSet().select<SetupMolInfo*>();
  if( moldat.size()==1 ) {
    log<<"  MOLINFO DATA found, using proper atom names\n";
    for(unsigned i=0; i<atoms.size(); ++i) {
      string Aname = moldat[0]->getAtomName(atoms[i]);
      string Rname = moldat[0]->getResidueName(atoms[i]);
      if(Rname=="ALA") {
        if(Aname=="BB") {
          atoi[i]=ALA_BB;
        } else error("Atom name not known: "+Aname);
      } else if(Rname=="ARG") {
        if(Aname=="BB") {
          atoi[i]=ARG_BB;
        } else if(Aname=="SC1") {
          atoi[i]=ARG_SC1;
        } else if(Aname=="SC2") {
          atoi[i]=ARG_SC2;
        } else error("Atom name not known: "+Aname);
      } else if(Rname=="ASN") {
        if(Aname=="BB") {
          atoi[i]=ASN_BB;
        } else if(Aname=="SC1") {
          atoi[i]=ASN_SC1;
        } else error("Atom name not known: "+Aname);
      } else if(Rname=="ASP") {
        if(Aname=="BB") {
          atoi[i]=ASP_BB;
        } else if(Aname=="SC1") {
          atoi[i]=ASP_SC1;
        } else error("Atom name not known: "+Aname);
      } else if(Rname=="CYS") {
        if(Aname=="BB") {
          atoi[i]=CYS_BB;
        } else if(Aname=="SC1") {
          atoi[i]=CYS_SC1;
        } else error("Atom name not known: "+Aname);
      } else if(Rname=="GLN") {
        if(Aname=="BB") {
          atoi[i]=GLN_BB;
        } else if(Aname=="SC1") {
          atoi[i]=GLN_SC1;
        } else error("Atom name not known: "+Aname);
      } else if(Rname=="GLU") {
        if(Aname=="BB") {
          atoi[i]=GLU_BB;
        } else if(Aname=="SC1") {
          atoi[i]=GLU_SC1;
        } else error("Atom name not known: "+Aname);
      } else if(Rname=="GLY") {
        if(Aname=="BB") {
          atoi[i]=GLY_BB;
        } else error("Atom name not known: "+Aname);
      } else if(Rname=="HIS") {
        if(Aname=="BB") {
          atoi[i]=HIS_BB;
        } else if(Aname=="SC1") {
          atoi[i]=HIS_SC1;
        } else if(Aname=="SC2") {
          atoi[i]=HIS_SC2;
        } else if(Aname=="SC3") {
          atoi[i]=HIS_SC3;
        } else error("Atom name not known: "+Aname);
      } else if(Rname=="ILE") {
        if(Aname=="BB") {
          atoi[i]=ILE_BB;
        } else if(Aname=="SC1") {
          atoi[i]=ILE_SC1;
        } else error("Atom name not known: "+Aname);
      } else if(Rname=="LEU") {
        if(Aname=="BB") {
          atoi[i]=LEU_BB;
        } else if(Aname=="SC1") {
          atoi[i]=LEU_SC1;
        } else error("Atom name not known: "+Aname);
      } else if(Rname=="LYS") {
        if(Aname=="BB") {
          atoi[i]=LYS_BB;
        } else if(Aname=="SC1") {
          atoi[i]=LYS_SC1;
        } else if(Aname=="SC2") {
          atoi[i]=LYS_SC2;
        } else error("Atom name not known: "+Aname);
      } else if(Rname=="MET") {
        if(Aname=="BB") {
          atoi[i]=MET_BB;
        } else if(Aname=="SC1") {
          atoi[i]=MET_SC1;
        } else error("Atom name not known: "+Aname);
      } else if(Rname=="PHE") {
        if(Aname=="BB") {
          atoi[i]=PHE_BB;
        } else if(Aname=="SC1") {
          atoi[i]=PHE_SC1;
        } else if(Aname=="SC2") {
          atoi[i]=PHE_SC2;
        } else if(Aname=="SC3") {
          atoi[i]=PHE_SC3;
        } else error("Atom name not known: "+Aname);
      } else if(Rname=="PRO") {
        if(Aname=="BB") {
          atoi[i]=PRO_BB;
        } else if(Aname=="SC1") {
          atoi[i]=PRO_SC1;
        } else error("Atom name not known: "+Aname);
      } else if(Rname=="SER") {
        if(Aname=="BB") {
          atoi[i]=SER_BB;
        } else if(Aname=="SC1") {
          atoi[i]=SER_SC1;
        } else error("Atom name not known: "+Aname);
      } else if(Rname=="THR") {
        if(Aname=="BB") {
          atoi[i]=THR_BB;
        } else if(Aname=="SC1") {
          atoi[i]=THR_SC1;
        } else error("Atom name not known: "+Aname);
      } else if(Rname=="TRP") {
        if(Aname=="BB") {
          atoi[i]=TRP_BB;
        } else if(Aname=="SC1") {
          atoi[i]=TRP_SC1;
        } else if(Aname=="SC2") {
          atoi[i]=TRP_SC2;
        } else if(Aname=="SC3") {
          atoi[i]=TRP_SC3;
        } else if(Aname=="SC4") {
          atoi[i]=TRP_SC4;
        } else error("Atom name not known: "+Aname);
      } else if(Rname=="TYR") {
        if(Aname=="BB") {
          atoi[i]=TYR_BB;
        } else if(Aname=="SC1") {
          atoi[i]=TYR_SC1;
        } else if(Aname=="SC2") {
          atoi[i]=TYR_SC2;
        } else if(Aname=="SC3") {
          atoi[i]=TYR_SC3;
        } else error("Atom name not known: "+Aname);
      } else if(Rname=="VAL") {
        if(Aname=="BB") {
          atoi[i]=VAL_BB;
        } else if(Aname=="SC1") {
          atoi[i]=VAL_SC1;
        } else error("Atom name not known: "+Aname);
      } else if(Rname=="  A") {
        if(Aname=="BB1") {
          atoi[i]=A_BB1;
        } else if(Aname=="BB2") {
          atoi[i]=A_BB2;
        } else if(Aname=="BB3") {
          atoi[i]=A_BB3;
        } else if(Aname=="SC1") {
          atoi[i]=A_SC1;
        } else if(Aname=="SC2") {
          atoi[i]=A_SC2;
        } else if(Aname=="SC3") {
          atoi[i]=A_SC3;
        } else if(Aname=="SC4") {
          atoi[i]=A_SC4;
        } else if(Aname=="3TE") {
          atoi[i]=A_3TE;
        } else if(Aname=="5TE") {
          atoi[i]=A_5TE;
        } else if(Aname=="TE3") {
          atoi[i]=A_TE3;
        } else if(Aname=="TE5") {
          atoi[i]=A_TE5;
        } else error("Atom name not known: "+Aname);
      } else if(Rname=="  C") {
        if(Aname=="BB1") {
          atoi[i]=C_BB1;
        } else if(Aname=="BB2") {
          atoi[i]=C_BB2;
        } else if(Aname=="BB3") {
          atoi[i]=C_BB3;
        } else if(Aname=="SC1") {
          atoi[i]=C_SC1;
        } else if(Aname=="SC2") {
          atoi[i]=C_SC2;
        } else if(Aname=="SC3") {
          atoi[i]=C_SC3;
        } else if(Aname=="3TE") {
          atoi[i]=C_3TE;
        } else if(Aname=="5TE") {
          atoi[i]=C_5TE;
        } else if(Aname=="TE3") {
          atoi[i]=C_TE3;
        } else if(Aname=="TE5") {
          atoi[i]=C_TE5;
        } else error("Atom name not known: "+Aname);
      } else if(Rname=="  G") {
        if(Aname=="BB1") {
          atoi[i]=G_BB1;
        } else if(Aname=="BB2") {
          atoi[i]=G_BB2;
        } else if(Aname=="BB3") {
          atoi[i]=G_BB3;
        } else if(Aname=="SC1") {
          atoi[i]=G_SC1;
        } else if(Aname=="SC2") {
          atoi[i]=G_SC2;
        } else if(Aname=="SC3") {
          atoi[i]=G_SC3;
        } else if(Aname=="SC4") {
          atoi[i]=G_SC4;
        } else if(Aname=="3TE") {
          atoi[i]=G_3TE;
        } else if(Aname=="5TE") {
          atoi[i]=G_5TE;
        } else if(Aname=="TE3") {
          atoi[i]=G_TE3;
        } else if(Aname=="TE5") {
          atoi[i]=G_TE5;
        } else error("Atom name not known: "+Aname);
      } else if(Rname=="  U") {
        if(Aname=="BB1") {
          atoi[i]=U_BB1;
        } else if(Aname=="BB2") {
          atoi[i]=U_BB2;
        } else if(Aname=="BB3") {
          atoi[i]=U_BB3;
        } else if(Aname=="SC1") {
          atoi[i]=U_SC1;
        } else if(Aname=="SC2") {
          atoi[i]=U_SC2;
        } else if(Aname=="SC3") {
          atoi[i]=U_SC3;
        } else if(Aname=="3TE") {
          atoi[i]=U_3TE;
        } else if(Aname=="5TE") {
          atoi[i]=U_5TE;
        } else if(Aname=="TE3") {
          atoi[i]=U_TE3;
        } else if(Aname=="TE5") {
          atoi[i]=U_TE5;
        } else error("Atom name not known: "+Aname);
      } else if(Rname==" DA") {
        if(Aname=="BB1") {
          atoi[i]=DA_BB1;
        } else if(Aname=="BB2") {
          atoi[i]=DA_BB2;
        } else if(Aname=="BB3") {
          atoi[i]=DA_BB3;
        } else if(Aname=="SC1") {
          atoi[i]=DA_SC1;
        } else if(Aname=="SC2") {
          atoi[i]=DA_SC2;
        } else if(Aname=="SC3") {
          atoi[i]=DA_SC3;
        } else if(Aname=="SC4") {
          atoi[i]=DA_SC4;
        } else if(Aname=="3TE") {
          atoi[i]=DA_3TE;
        } else if(Aname=="5TE") {
          atoi[i]=DA_5TE;
        } else if(Aname=="TE3") {
          atoi[i]=DA_TE3;
        } else if(Aname=="TE5") {
          atoi[i]=DA_TE5;
        } else error("Atom name not known: "+Aname);
      } else if(Rname==" DC") {
        if(Aname=="BB1") {
          atoi[i]=DC_BB1;
        } else if(Aname=="BB2") {
          atoi[i]=DC_BB2;
        } else if(Aname=="BB3") {
          atoi[i]=DC_BB3;
        } else if(Aname=="SC1") {
          atoi[i]=DC_SC1;
        } else if(Aname=="SC2") {
          atoi[i]=DC_SC2;
        } else if(Aname=="SC3") {
          atoi[i]=DC_SC3;
        } else if(Aname=="3TE") {
          atoi[i]=DC_3TE;
        } else if(Aname=="5TE") {
          atoi[i]=DC_5TE;
        } else if(Aname=="TE3") {
          atoi[i]=DC_TE3;
        } else if(Aname=="TE5") {
          atoi[i]=DC_TE5;
        } else error("Atom name not known: "+Aname);
      } else if(Rname==" DG") {
        if(Aname=="BB1") {
          atoi[i]=DG_BB1;
        } else if(Aname=="BB2") {
          atoi[i]=DG_BB2;
        } else if(Aname=="BB3") {
          atoi[i]=DG_BB3;
        } else if(Aname=="SC1") {
          atoi[i]=DG_SC1;
        } else if(Aname=="SC2") {
          atoi[i]=DG_SC2;
        } else if(Aname=="SC3") {
          atoi[i]=DG_SC3;
        } else if(Aname=="SC4") {
          atoi[i]=DG_SC4;
        } else if(Aname=="3TE") {
          atoi[i]=DG_3TE;
        } else if(Aname=="5TE") {
          atoi[i]=DG_5TE;
        } else if(Aname=="TE3") {
          atoi[i]=DG_TE3;
        } else if(Aname=="TE5") {
          atoi[i]=DG_TE5;
        } else error("Atom name not known: "+Aname);
      } else if(Rname==" DT") {
        if(Aname=="BB1") {
          atoi[i]=DT_BB1;
        } else if(Aname=="BB2") {
          atoi[i]=DT_BB2;
        } else if(Aname=="BB3") {
          atoi[i]=DT_BB3;
        } else if(Aname=="SC1") {
          atoi[i]=DT_SC1;
        } else if(Aname=="SC2") {
          atoi[i]=DT_SC2;
        } else if(Aname=="SC3") {
          atoi[i]=DT_SC3;
        } else if(Aname=="3TE") {
          atoi[i]=DT_3TE;
        } else if(Aname=="5TE") {
          atoi[i]=DT_5TE;
        } else if(Aname=="TE3") {
          atoi[i]=DT_TE3;
        } else if(Aname=="TE5") {
          atoi[i]=DT_TE5;
        } else error("Atom name not known: "+Aname);
      } else error("Residue not known: "+Rname);
    }
  } else {
    error("MOLINFO DATA not found\n");
  }
}

double SAXS::calculateASF(const vector<AtomNumber> &atoms, vector<vector<long double> > &FF_tmp, const double rho)
{
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

  double Iq0=0.;
  if( moldat.size()==1 ) {
    log<<"  MOLINFO DATA found, using proper atom names\n";
    // cycle over the atom types
    for(unsigned i=0; i<NTT; i++) {
      const double volr = pow(param_v[i], (2.0/3.0)) /(4. * M_PI);
      // cycle over q
      for(unsigned k=0; k<q_list.size(); ++k) {
        const double q = q_list[k];
        const double s = q / (4. * M_PI);
        FF_tmp[k][i] = param_c[i];
        // SUM [a_i * EXP( - b_i * (q/4pi)^2 )] Waasmaier and Kirfel (1995)
        for(unsigned j=0; j<4; j++) {
          FF_tmp[k][i] += param_a[i][j]*exp(-param_b[i][j]*s*s);
        }
        // subtract solvation: rho * v_i * EXP( (- v_i^(2/3) / (4pi)) * q^2  ) // since  D in Fraser 1978 is 2*s
        FF_tmp[k][i] -= rho*param_v[i]*exp(-volr*q*q);
      }
    }
    // cycle over the atoms to assign the atom type and calculate I0
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
        atoi[i] = AA_map[type_s];
        for(unsigned j=0; j<4; j++) Iq0 += param_a[index][j];
        Iq0 = Iq0 -rho*param_v[index] + param_c[index];
      } else {
        error("Wrong atom type "+type_s+" from atom name "+name+"\n");
      }
    }
  } else {
    error("MOLINFO DATA not found\n");
  }

  return Iq0;
}

}
}
