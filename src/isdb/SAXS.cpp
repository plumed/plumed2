/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2017-2019 The plumed team
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

#ifdef __PLUMED_HAS_GSL
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_legendre.h>
#endif

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
Calculates SAXS scattered intensity using either the Debye equation or the harmonic sphere approximation.

Intensities are calculated for a set of scattering length set using QVALUE keywords that are numbered starting from 0.
Structure factors can be either assigned using a polynomial expansion to any order using the PARAMETERS keywords;
automatically assigned to atoms using the ATOMISTIC flag reading a PDB file, a correction for the water density is
automatically added, with water density that by default is 0.334 but that can be set otherwise using WATERDENS;
automatically assigned to Martini pseudo atoms using the MARTINI flag.
The calculated intensities can be scaled using the SCALEINT keywords. This is applied by rescaling the structure factors.
Experimental reference intensities can be added using the EXPINT keywords.
By default SAXS is calculated using Debye on CPU, by adding the GPU flag it is possible to solve the equation on a GPU
if the ARRAYFIRE libraries are installed and correctly linked (). Alternatively we an implementation based on Bessel functions,
BESSEL flag. This is very fast for small q values because a short expansion is enough.
An automatic choice is made for which q Bessel are used and for which the calculation is done by Debye. If one wants to force
all q values to be calculated using Bessel function this can be done using FORCE_BESSEL.
Irrespective of the method employed, \ref METAINFERENCE can be activated using DOSCORE and the other relevant keywords.

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
  bool                       pbc;
  bool                       serial;
  bool                       bessel;
  bool                       force_bessel;
  bool                       gpu;
  int                        deviceid;
  vector<double>             q_list;
  vector<double>             FF_rank;
  vector<vector<double> >    FF_value;
  vector<vector<float> >     FFf_value;
  vector<double>             avals;
  vector<double>             bvals;

  void calculate_gpu(vector<Vector> &deriv);
  void calculate_cpu(vector<Vector> &deriv);
  void getMartiniSFparam(const vector<AtomNumber> &atoms, vector<vector<long double> > &parameter);
  double calculateASF(const vector<AtomNumber> &atoms, vector<vector<long double> > &FF_tmp, const double rho);
  void bessel_calculate(vector<Vector> &deriv, vector<double> &sum, vector<Vector2d> &qRnm, const vector<double> &r_polar,
                        const vector<unsigned> &trunc, const int algorithm, const unsigned p2);
  void setup_midl(vector<double> &r_polar, vector<Vector2d> &qRnm, int &algorithm, unsigned &p2, vector<unsigned> &trunc);
  Vector2d dXHarmonics(const unsigned p2, const unsigned k, const unsigned i, const int n, const int m, const vector<Vector2d> &decRnm);
  Vector2d dYHarmonics(const unsigned p2, const unsigned k, const unsigned i, const int n, const int m, const vector<Vector2d> &decRnm);
  Vector2d dZHarmonics(const unsigned p2, const unsigned k, const unsigned i, const int n, const int m, const vector<Vector2d> &decRnm);
  void cal_coeff();

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
  keys.addFlag("BESSEL",false,"Perform the calculation using the adaptive spherical harmonic approximation");
  keys.addFlag("FORCE_BESSEL",false,"Perform the calculation using the adaptive spherical harmonic approximation, without adaptive algorithm, useful for debug only");
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
  bessel(false),
  force_bessel(false),
  gpu(false),
  deviceid(0)
{
  vector<AtomNumber> atoms;
  parseAtomList("ATOMS",atoms);
  const unsigned size = atoms.size();

  parseFlag("SERIAL",serial);

  parseFlag("BESSEL",bessel);
  parseFlag("FORCE_BESSEL",force_bessel);
  if(force_bessel) bessel = true;
#ifndef __PLUMED_HAS_GSL
  if(bessel) error("You CANNOT use BESSEL without GSL. Recompile PLUMED with GSL!\n");
#endif
  if(bessel) cal_coeff();

  bool nopbc=!pbc;
  parseFlag("NOPBC",nopbc);
  pbc=!nopbc;

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

  if(bessel&&gpu) error("You CANNOT use BESSEL on GPU!\n");

  unsigned ntarget=0;
  for(unsigned i=0;; ++i) {
    double t_list;
    if( !parseNumbered( "QVALUE", i+1, t_list) ) break;
    if(t_list<=0.) error("QVALUE cannot be less or equal to zero!\n");
    q_list.push_back(t_list);
    ntarget++;
  }
  const unsigned numq = ntarget;

  bool atomistic=false;
  parseFlag("ATOMISTIC",atomistic);
  bool martini=false;
  parseFlag("MARTINI",martini);

  if(martini&&atomistic) error("You cannot use MARTINI and ATOMISTIC at the same time");

  double rho = 0.334;
  parse("WATERDENS", rho);

  double Iq0=0;
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
    for(unsigned i=0; i<size; ++i) Iq0+=parameter[i][0];
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
    for(unsigned i=0; i<size; ++i) Iq0+=parameter[i][0];
  } else if(atomistic) {
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


  if(pbc)      log.printf("  using periodic boundary conditions\n");
  else         log.printf("  without periodic boundary conditions\n");
  for(unsigned i=0; i<numq; i++) {
    if(q_list[i]==0.) error("it is not possible to set q=0\n");
    if(i>0&&q_list[i]<q_list[i-1]) error("QVALUE must be in ascending order");
    log.printf("  my q: %lf \n",q_list[i]);
  }

  // Calculate Rank of FF_matrix
  if(tmp_scale_int!=1) scale_int /= tmp_scale_int;
  else {
    if(exp) scale_int /= expint[0];
  }

  if(!gpu) {
    FF_rank.resize(numq);
    FF_value.resize(numq,vector<double>(size));
    for(unsigned k=0; k<numq; ++k) {
      for(unsigned i=0; i<size; i++) {
        FF_value[k][i] = static_cast<double>(FF_tmp[k][i])/sqrt(scale_int);
        FF_rank[k]+=FF_value[k][i]*FF_value[k][i];
      }
    }
  } else {
    FFf_value.resize(numq,vector<float>(size));
    for(unsigned k=0; k<numq; ++k) {
      for(unsigned i=0; i<size; i++) {
        FFf_value[k][i] = static_cast<float>(FF_tmp[k][i])/sqrt(scale_int);
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
    if(bessel&&i>0&&q_list[i]<q_list[i-1]) plumed_merror("With BESSEL the Q values should be ordered from the smallest to the largest");
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
  if(bessel) log<<plumed.cite("Gumerov, Berlin, Fushman, Duraiswami, J. Comput. Chem. 33, 1981-1996 (2012).");
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
  vector<Vector> c_dist(size*size);
  vector<double> m_dist(size*size);

  vector<double> r_polar;
  vector<Vector2d> qRnm;
  vector<unsigned> trunc;
  int algorithm=-1;
  unsigned p2=0;
  bool direct = true;

  if(bessel) {
    r_polar.resize(size);
    trunc.resize(numq);
    setup_midl(r_polar, qRnm, algorithm, p2, trunc);
    if(algorithm>=0) bessel_calculate(deriv, sum, qRnm, r_polar, trunc, algorithm, p2);
    if(algorithm+1>=numq) direct=false;
    if(algorithm==-1) bessel=false;
  }

  if(direct) {
    #pragma omp parallel for num_threads(OpenMP::getNumThreads())
    for (unsigned i=rank; i<size-1; i+=stride) {
      const Vector posi=getPosition(i);
      for (unsigned j=i+1; j<size ; j++) {
        c_dist[i*size+j] = delta(posi,getPosition(j));
        m_dist[i*size+j] = c_dist[i*size+j].modulo();
      }
    }

    #pragma omp parallel for num_threads(OpenMP::getNumThreads())
    for (unsigned k=(algorithm+1); k<numq; k++) {
      const unsigned kdx=k*size;
      for (unsigned i=rank; i<size-1; i+=stride) {
        const double FF=2.*FF_value[k][i];
        Vector dsum;
        for (unsigned j=i+1; j<size ; j++) {
          const Vector c_distances = c_dist[i*size+j];
          const double m_distances = m_dist[i*size+j];
          const double qdist       = q_list[k]*m_distances;
          const double FFF = FF*FF_value[k][j];
          const double tsq = FFF*sin(qdist)/qdist;
          const double tcq = FFF*cos(qdist);
          const double tmp = (tcq-tsq)/(m_distances*m_distances);
          const Vector dd  = c_distances*tmp;
          dsum         += dd;
          deriv[kdx+j] += dd;
          sum[k]       += tsq;
        }
        deriv[kdx+i] -= dsum;
      }
    }
  }

  if(!serial) {
    comm.Sum(&deriv[0][0], 3*deriv.size());
    comm.Sum(&sum[0], numq);
  }

  if(bessel) {
    for(unsigned k=0; k<=algorithm; k++) {
      const unsigned kN = k*size;
      sum[k] *= 4.*M_PI;
      string num; Tools::convert(k,num);
      Value* val=getPntrToComponent("q-"+num);
      val->set(sum[k]);
      if(getDoScore()) setCalcData(k, sum[k]);
      for(unsigned i=0; i<size; i++) deriv[kN+i] *= 8.*M_PI*q_list[k];
    }
  }

  if(direct) {
    for (unsigned k=algorithm+1; k<numq; k++) {
      sum[k]+=FF_rank[k];
      string num; Tools::convert(k,num);
      Value* val=getPntrToComponent("q-"+num);
      val->set(sum[k]);
      if(getDoScore()) setCalcData(k, sum[k]);
    }
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

void SAXS::bessel_calculate(vector<Vector> &deriv, vector<double> &sum, vector<Vector2d> &qRnm, const vector<double> &r_polar,
                            const vector<unsigned> &trunc, const int algorithm, const unsigned p2)
{
#ifdef __PLUMED_HAS_GSL
  const unsigned size = getNumberOfAtoms();

  unsigned stride = comm.Get_size();
  unsigned rank   = comm.Get_rank();
  if(serial) {
    stride = 1;
    rank   = 0;
  }

  //calculation via Middleman method
  for(unsigned k=0; k<algorithm+1; k++) {
    const unsigned kN  = k * size;
    const unsigned p22 = trunc[k]*trunc[k];
    //double sum over the p^2 expansion terms
    vector<Vector2d> Bnm(p22);
    for(unsigned i=rank; i<size; i+=stride) {
      double pq = r_polar[i]*q_list[k];
      for(unsigned n=0; n<trunc[k]; n++) {
        //the spherical bessel functions do not depend on the order and are therefore precomputed here
        double besself = gsl_sf_bessel_jl(n,pq);
        //here conj(R(m,n))=R(-m,n) is used to decrease the terms in the sum over m by a factor of two
        for(unsigned m=0; m<(n+1); m++) {
          int order = m-n;
          int s = n*n + m;
          int t = s - 2*order;
          int x = p2*i + s;
          int y = p2*i + t;
          //real part of the spherical basis function of order m, degree n of atom i
          qRnm[x]  *= besself;
          //real part of the spherical basis function of order -m, degree n of atom i
          qRnm[y][0] = qRnm[x][0];
          //imaginary part of the spherical basis function of order -m, degree n of atom i
          qRnm[y][1] = -qRnm[x][1];
          //expansion coefficient of order m and degree n
          Bnm[s] += FF_value[k][i] * qRnm[y];
          //correction for expansion coefficient of order -m and degree n
          if(order!=0) Bnm[t] += FF_value[k][i] * qRnm[x];
        }
      }
    }

    //calculate expansion coefficients for the derivatives
    vector<Vector2d> a(3*p22);
    for(unsigned i=rank; i<size; i+=stride) {
      for(unsigned n=0; n<trunc[k]-1; n++) {
        for(unsigned m=0; m<(2*n)+1; m++) {
          unsigned t = 3*(n*n + m);
          a[t]   += FF_value[k][i] * dXHarmonics(p2,k,i,n,m,qRnm);
          a[t+1] += FF_value[k][i] * dYHarmonics(p2,k,i,n,m,qRnm);
          a[t+2] += FF_value[k][i] * dZHarmonics(p2,k,i,n,m,qRnm);
        }
      }
    }
    if(!serial) {
      comm.Sum(&Bnm[0][0],2*p22);
      comm.Sum(&a[0][0],  6*p22);
    }

    //calculation of the scattering profile I of the kth scattering wavenumber q
    for(int n=rank; n<trunc[k]; n+=stride) {
      for(int m=0; m<(2*n)+1; m++) {
        int s = n * n + m;
        sum[k] += Bnm[s][0]*Bnm[s][0] + Bnm[s][1]*Bnm[s][1];
      }
    }

    //calculation of (box)derivatives
    for(unsigned i=rank; i<size; i+=stride) {
      //vector of the derivatives of the expanded functions psi
      Vector dPsi;
      int s = p2 * i;
      double pq = r_polar[i]* q_list[k];
      for(int n=trunc[k]-1; n>=0; n--) {
        double besself = gsl_sf_bessel_jl(n,pq);
        for(int m=0; m<(2*n)+1; m++) {
          int y = n  * n + m  + s;
          int z = 3*(n*n+m);
          dPsi[0] += 0.5*(qRnm[y][0] * a[z][0]   + qRnm[y][1] * a[z][1]);
          dPsi[1] += 0.5*(qRnm[y][0] * a[z+1][1] - qRnm[y][1] * a[z+1][0]);
          dPsi[2] +=      qRnm[y][0] * a[z+2][0] + qRnm[y][1] * a[z+2][1];
          qRnm[y] /= besself;
        }
      }
      deriv[kN+i] += FF_value[k][i] * dPsi;
    }
  }
  //end of the k loop
#endif
}

void SAXS::setup_midl(vector<double> &r_polar, vector<Vector2d> &qRnm, int &algorithm, unsigned &p2, vector<unsigned> &trunc)
{
#ifdef __PLUMED_HAS_GSL
  const unsigned size = getNumberOfAtoms();
  const unsigned numq = q_list.size();

  unsigned stride = comm.Get_size();
  unsigned rank   = comm.Get_rank();
  if(serial) {
    stride = 1;
    rank   = 0;
  }

  Vector max = getPosition(0);
  Vector min = getPosition(0);
  vector<Vector> polar(size);

  // transform in polar and look for min and max dist
  for(unsigned i=0; i<size; i++) {
    Vector coord=getPosition(i);
    //r
    polar[i][0]=sqrt(coord[0]*coord[0]+coord[1]*coord[1]+coord[2]*coord[2]);
    r_polar[i] = polar[i][0];
    //cos(theta)
    polar[i][1]=coord[2]/polar[i][0];
    //phi
    polar[i][2]=atan2(coord[1],coord[0]);

    if(coord[0]<min[0]) min[0] = coord[0];
    if(coord[1]<min[1]) min[1] = coord[1];
    if(coord[2]<min[2]) min[2] = coord[2];
    if(coord[0]>max[0]) max[0] = coord[0];
    if(coord[1]>max[1]) max[1] = coord[1];
    if(coord[2]>max[2]) max[2] = coord[2];
  }
  max -= min;
  double maxdist = max[0];
  if(maxdist<max[1]) maxdist = max[1];
  if(maxdist<max[2]) maxdist = max[2];
  unsigned truncation=5+static_cast<unsigned>(1.2*maxdist*q_list[numq-1]+0.5*pow((12-log10(maxdist*q_list[numq-1])),2/3)*pow(maxdist*q_list[numq-1],1/3));
  if(truncation<10) truncation=10;
  if(truncation>99) truncation=99;
  p2=truncation*truncation;
  //dynamically set the truncation according to the scattering wavenumber.
  for(int k=numq-1; k>=0; k--) {
    trunc[k]=5+static_cast<unsigned>(1.2*maxdist*q_list[k]+0.5*pow((12-log10(maxdist*q_list[k])),2/3)*pow(maxdist*q_list[k],1/3));
    if(trunc[k]<10) trunc[k] = 10;
    if(4*trunc[k]<static_cast<unsigned>(sqrt(2*size)) && algorithm<0) algorithm=k;
  }

  if(algorithm==-1) log.printf("BESSEL is suboptimal for this system and is being disabled, unless FORCE_BESSEL is used\n");
  if(force_bessel) algorithm=numq-1;

  unsigned qRnm_size = p2*size;
  qRnm.resize(qRnm_size);
  //as the legndre polynomials and the exponential term in the basis set expansion are not function of the scattering wavenumber, they can be precomputed
  for(unsigned i=rank; i<size; i+=stride) {
    for(int n=0; n<truncation; n++) {
      for(int m=0; m<(n+1); m++) {
        int order  = m-n;
        int x      = p2*i + n*n + m;
        double gsl = gsl_sf_legendre_sphPlm(n,abs(order),polar[i][1]);
        //real part of the spherical basis function of order m, degree n of atom i
        qRnm[x][0] = gsl * cos(order*polar[i][2]);
        //imaginary part of the spherical basis function of order m, degree n of atom i
        qRnm[x][1] = gsl * sin(order*polar[i][2]);
      }
    }
  }
#endif
}

void SAXS::update() {
  // write status file
  if(getWstride()>0&& (getStep()%getWstride()==0 || getCPT()) ) writeStatus();
}

//partial derivatives of the spherical basis functions
Vector2d SAXS::dXHarmonics(const unsigned p2, const unsigned k, const unsigned i, const int n, const int m, const vector<Vector2d> &qRnm) {
  Vector2d                            dRdc = (bvals[n*(n+4)-m+1] * qRnm[p2*i+n*(n+2)+m+3] + bvals[n*(n+2)+m+1] * qRnm[p2*i+n*(n+2)+m+1]);
  if((abs(m-n-1)<=(n-1))&&((n-1)>=0)) dRdc-= bvals[n*(n+2)-m] * qRnm[p2*i+n*(n-2)+m-1];
  if((abs(m-n+1)<=(n-1))&&((n-1)>=0)) dRdc-= bvals[n*n+m] * qRnm[p2*i+n*n-2*n+m+1];
  return dRdc;
}


Vector2d SAXS::dYHarmonics(const unsigned p2, const unsigned k, const unsigned i, const int n, const int m, const vector<Vector2d> &qRnm) {
  Vector2d                            dRdc = (bvals[n*(n+4)-m+1] * qRnm[p2*i+n*(n+2)+m+3] - bvals[n*(n+2)+m+1] * qRnm[p2*i+n*(n+2)+m+1]);
  if((abs(m-n-1)<=(n-1))&&((n-1)>=0)) dRdc+= bvals[n*(n+2)-m] * qRnm[p2*i+n*(n-2)+m-1];
  if((abs(m-n+1)<=(n-1))&&((n-1)>=0)) dRdc-= bvals[n*n+m] * qRnm[p2*i+n*(n-2)+m+1];
  return dRdc;
}


Vector2d SAXS::dZHarmonics(const unsigned p2, const unsigned k, const unsigned i, const int n, const int m, const vector<Vector2d> &qRnm) {
  Vector2d                          dRdc = -avals[n*n+m]*qRnm[p2*i+n*(n+2)+m+2];
  if((abs(m-n)<=(n-1))&&((n-1)>=0)) dRdc+= avals[n*(n-2)+m]*qRnm[p2*i+n*(n-2)+m];
  return dRdc;
}

//coefficients for partial derivatives of the spherical basis functions
void SAXS::cal_coeff() {
  avals.resize(100*100);
  bvals.resize(100*100);
  for(int n=0; n<100; n++) {
    for(int m=0; m<(2*n)+1; m++) {
      double mval = m-n;
      double nval = n;
      avals[n*n+m] = -1 * sqrt(((nval+mval+1)*(nval+1-mval))/(((2*nval)+1)*((2*nval)+3)));
      bvals[n*n+m] = sqrt(((nval-mval-1)*(nval-mval))/(((2*nval)-1)*((2*nval)+1)));
      if((-n<=(m-n)) && ((m-n)<0)) bvals[n*n+m] *= -1;
    }
  }
}

void SAXS::getMartiniSFparam(const vector<AtomNumber> &atoms, vector<vector<long double> > &parameter)
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
        } else error("Atom name not known: "+Aname);
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
        } else error("Atom name not known: "+Aname);
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
        } else error("Atom name not known: "+Aname);
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
        } else error("Atom name not known: "+Aname);
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
        } else error("Atom name not known: "+Aname);
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
        } else error("Atom name not known: "+Aname);
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
        } else error("Atom name not known: "+Aname);
      } else if(Rname=="GLY") {
        if(Aname=="BB") {
          parameter[i].push_back(9.977);
          parameter[i].push_back(-0.0285799);
          parameter[i].push_back(1.84236);
          parameter[i].push_back(-0.0315192);
          parameter[i].push_back(-2.88326);
          parameter[i].push_back(1.87323);
          parameter[i].push_back(-0.345773);
        } else error("Atom name not known: "+Aname);
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
        } else error("Atom name not known: "+Aname);
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
        } else error("Atom name not known: "+Aname);
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
        } else error("Atom name not known: "+Aname);
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
        } else error("Atom name not known: "+Aname);
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
        } else error("Atom name not known: "+Aname);
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
        } else error("Atom name not known: "+Aname);
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
        } else error("Atom name not known: "+Aname);
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
        } else error("Atom name not known: "+Aname);
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
        } else error("Atom name not known: "+Aname);
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
        } else error("Atom name not known: "+Aname);
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
        } else error("Atom name not known: "+Aname);
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
        } else error("Atom name not known: "+Aname);
      } else if(Rname=="  A") {
        if(Aname=="BB1") {
          parameter[i].push_back(32.88500000);
          parameter[i].push_back(0.08339900);
          parameter[i].push_back(-7.36054400);
          parameter[i].push_back(2.19220300);
          parameter[i].push_back(-3.56523400);
          parameter[i].push_back(2.33326900);
          parameter[i].push_back(-0.39785500);
        } else if(Aname=="BB2") {
          parameter[i].push_back(3.80600000);
          parameter[i].push_back(-0.10727600);
          parameter[i].push_back(9.58854100);
          parameter[i].push_back(-6.23740500);
          parameter[i].push_back(-0.48267300);
          parameter[i].push_back(1.14119500);
          parameter[i].push_back(-0.21385600);
        } else if(Aname=="BB3") {
          parameter[i].push_back(3.59400000);
          parameter[i].push_back(0.04537300);
          parameter[i].push_back(9.59178900);
          parameter[i].push_back(-1.29202200);
          parameter[i].push_back(-7.10851000);
          parameter[i].push_back(4.05571200);
          parameter[i].push_back(-0.63372500);
        } else if(Aname=="SC1") {
          parameter[i].push_back(6.67100000);
          parameter[i].push_back(-0.00855300);
          parameter[i].push_back(1.63222400);
          parameter[i].push_back(-0.06466200);
          parameter[i].push_back(-1.48694200);
          parameter[i].push_back(0.78544600);
          parameter[i].push_back(-0.12083500);
        } else if(Aname=="SC2") {
          parameter[i].push_back(5.95100000);
          parameter[i].push_back(-0.02606600);
          parameter[i].push_back(2.54399900);
          parameter[i].push_back(-0.48436900);
          parameter[i].push_back(-1.55357400);
          parameter[i].push_back(0.86466900);
          parameter[i].push_back(-0.13509000);
        } else if(Aname=="SC3") {
          parameter[i].push_back(11.39400000);
          parameter[i].push_back(0.00871300);
          parameter[i].push_back(-0.23891300);
          parameter[i].push_back(0.48919400);
          parameter[i].push_back(-1.75289400);
          parameter[i].push_back(0.99267500);
          parameter[i].push_back(-0.16291300);
        } else if(Aname=="SC4") {
          parameter[i].push_back(6.45900000);
          parameter[i].push_back(0.01990600);
          parameter[i].push_back(4.17970400);
          parameter[i].push_back(0.97629900);
          parameter[i].push_back(-5.03297800);
          parameter[i].push_back(2.55576700);
          parameter[i].push_back(-0.39150500);
        } else if(Aname=="3TE") {
          parameter[i].push_back(4.23000000);
          parameter[i].push_back(0.00064800);
          parameter[i].push_back(0.92124600);
          parameter[i].push_back(0.08064300);
          parameter[i].push_back(-0.39054400);
          parameter[i].push_back(0.12429100);
          parameter[i].push_back(-0.01122700);
        } else if(Aname=="5TE") {
          parameter[i].push_back(4.23000000);
          parameter[i].push_back(0.00039300);
          parameter[i].push_back(0.92305100);
          parameter[i].push_back(0.07747500);
          parameter[i].push_back(-0.38792100);
          parameter[i].push_back(0.12323800);
          parameter[i].push_back(-0.01106600);
        } else if(Aname=="TE3") {
          parameter[i].push_back(7.82400000);
          parameter[i].push_back(-0.04881000);
          parameter[i].push_back(8.21557900);
          parameter[i].push_back(-0.89491400);
          parameter[i].push_back(-9.54293700);
          parameter[i].push_back(6.33122200);
          parameter[i].push_back(-1.16672900);
        } else if(Aname=="TE5") {
          parameter[i].push_back(8.03600000);
          parameter[i].push_back(0.01641200);
          parameter[i].push_back(5.14902200);
          parameter[i].push_back(0.83419700);
          parameter[i].push_back(-7.59068300);
          parameter[i].push_back(4.52063200);
          parameter[i].push_back(-0.78260800);
        } else error("Atom name not known: "+Aname);
      } else if(Rname=="  C") {
        if(Aname=="BB1") {
          parameter[i].push_back(32.88500000);
          parameter[i].push_back(0.08311100);
          parameter[i].push_back(-7.35432100);
          parameter[i].push_back(2.18610000);
          parameter[i].push_back(-3.55788300);
          parameter[i].push_back(2.32918700);
          parameter[i].push_back(-0.39720000);
        } else if(Aname=="BB2") {
          parameter[i].push_back(3.80600000);
          parameter[i].push_back(-0.10808100);
          parameter[i].push_back(9.61612600);
          parameter[i].push_back(-6.28595400);
          parameter[i].push_back(-0.45187000);
          parameter[i].push_back(1.13326000);
          parameter[i].push_back(-0.21320300);
        } else if(Aname=="BB3") {
          parameter[i].push_back(3.59400000);
          parameter[i].push_back(0.04484200);
          parameter[i].push_back(9.61919800);
          parameter[i].push_back(-1.33582800);
          parameter[i].push_back(-7.07200400);
          parameter[i].push_back(4.03952900);
          parameter[i].push_back(-0.63098200);
        } else if(Aname=="SC1") {
          parameter[i].push_back(5.95100000);
          parameter[i].push_back(-0.02911300);
          parameter[i].push_back(2.59700400);
          parameter[i].push_back(-0.55507700);
          parameter[i].push_back(-1.56344600);
          parameter[i].push_back(0.88956200);
          parameter[i].push_back(-0.14061300);
        } else if(Aname=="SC2") {
          parameter[i].push_back(11.62100000);
          parameter[i].push_back(0.01366100);
          parameter[i].push_back(-0.25959200);
          parameter[i].push_back(0.48918300);
          parameter[i].push_back(-1.52550500);
          parameter[i].push_back(0.83644100);
          parameter[i].push_back(-0.13407300);
        } else if(Aname=="SC3") {
          parameter[i].push_back(5.01900000);
          parameter[i].push_back(-0.03276100);
          parameter[i].push_back(5.53776900);
          parameter[i].push_back(-0.95105000);
          parameter[i].push_back(-3.71130800);
          parameter[i].push_back(2.16146000);
          parameter[i].push_back(-0.34918600);
        } else if(Aname=="3TE") {
          parameter[i].push_back(4.23000000);
          parameter[i].push_back(0.00057300);
          parameter[i].push_back(0.92174800);
          parameter[i].push_back(0.07964500);
          parameter[i].push_back(-0.38965700);
          parameter[i].push_back(0.12392500);
          parameter[i].push_back(-0.01117000);
        } else if(Aname=="5TE") {
          parameter[i].push_back(4.23000000);
          parameter[i].push_back(0.00071000);
          parameter[i].push_back(0.92082800);
          parameter[i].push_back(0.08150600);
          parameter[i].push_back(-0.39127000);
          parameter[i].push_back(0.12455900);
          parameter[i].push_back(-0.01126300);
        } else if(Aname=="TE3") {
          parameter[i].push_back(7.82400000);
          parameter[i].push_back(-0.05848300);
          parameter[i].push_back(8.29319900);
          parameter[i].push_back(-1.12563800);
          parameter[i].push_back(-9.42197600);
          parameter[i].push_back(6.35441700);
          parameter[i].push_back(-1.18356900);
        } else if(Aname=="TE5") {
          parameter[i].push_back(8.03600000);
          parameter[i].push_back(0.00493500);
          parameter[i].push_back(4.92622000);
          parameter[i].push_back(0.64810700);
          parameter[i].push_back(-7.05100000);
          parameter[i].push_back(4.26064400);
          parameter[i].push_back(-0.74819100);
        } else error("Atom name not known: "+Aname);
      } else if(Rname=="  G") {
        if(Aname=="BB1") {
          parameter[i].push_back(32.88500000);
          parameter[i].push_back(0.08325400);
          parameter[i].push_back(-7.35736000);
          parameter[i].push_back(2.18914800);
          parameter[i].push_back(-3.56154800);
          parameter[i].push_back(2.33120600);
          parameter[i].push_back(-0.39752300);
        } else if(Aname=="BB2") {
          parameter[i].push_back(3.80600000);
          parameter[i].push_back(-0.10788300);
          parameter[i].push_back(9.60930800);
          parameter[i].push_back(-6.27402500);
          parameter[i].push_back(-0.46192700);
          parameter[i].push_back(1.13737000);
          parameter[i].push_back(-0.21383100);
        } else if(Aname=="BB3") {
          parameter[i].push_back(3.59400000);
          parameter[i].push_back(0.04514500);
          parameter[i].push_back(9.61234700);
          parameter[i].push_back(-1.31542100);
          parameter[i].push_back(-7.09150500);
          parameter[i].push_back(4.04706200);
          parameter[i].push_back(-0.63201000);
        } else if(Aname=="SC1") {
          parameter[i].push_back(6.67100000);
          parameter[i].push_back(-0.00863200);
          parameter[i].push_back(1.63252300);
          parameter[i].push_back(-0.06567200);
          parameter[i].push_back(-1.48680500);
          parameter[i].push_back(0.78565600);
          parameter[i].push_back(-0.12088900);
        } else if(Aname=="SC2") {
          parameter[i].push_back(11.39400000);
          parameter[i].push_back(0.00912200);
          parameter[i].push_back(-0.22869000);
          parameter[i].push_back(0.49616400);
          parameter[i].push_back(-1.75039000);
          parameter[i].push_back(0.98649200);
          parameter[i].push_back(-0.16141600);
        } else if(Aname=="SC3") {
          parameter[i].push_back(10.90100000);
          parameter[i].push_back(0.02208700);
          parameter[i].push_back(0.17032800);
          parameter[i].push_back(0.73280800);
          parameter[i].push_back(-1.95292000);
          parameter[i].push_back(0.98357600);
          parameter[i].push_back(-0.14790900);
        } else if(Aname=="SC4") {
          parameter[i].push_back(6.45900000);
          parameter[i].push_back(0.02023700);
          parameter[i].push_back(4.17655400);
          parameter[i].push_back(0.98731800);
          parameter[i].push_back(-5.04352800);
          parameter[i].push_back(2.56059400);
          parameter[i].push_back(-0.39234300);
        } else if(Aname=="3TE") {
          parameter[i].push_back(4.23000000);
          parameter[i].push_back(0.00066300);
          parameter[i].push_back(0.92118800);
          parameter[i].push_back(0.08062700);
          parameter[i].push_back(-0.39041600);
          parameter[i].push_back(0.12419400);
          parameter[i].push_back(-0.01120500);
        } else if(Aname=="5TE") {
          parameter[i].push_back(4.23000000);
          parameter[i].push_back(0.00062800);
          parameter[i].push_back(0.92133500);
          parameter[i].push_back(0.08029900);
          parameter[i].push_back(-0.39015300);
          parameter[i].push_back(0.12411600);
          parameter[i].push_back(-0.01119900);
        } else if(Aname=="TE3") {
          parameter[i].push_back(7.82400000);
          parameter[i].push_back(-0.05177400);
          parameter[i].push_back(8.34606700);
          parameter[i].push_back(-1.02936300);
          parameter[i].push_back(-9.55211900);
          parameter[i].push_back(6.37776600);
          parameter[i].push_back(-1.17898000);
        } else if(Aname=="TE5") {
          parameter[i].push_back(8.03600000);
          parameter[i].push_back(0.00525100);
          parameter[i].push_back(4.71070600);
          parameter[i].push_back(0.66746900);
          parameter[i].push_back(-6.72538700);
          parameter[i].push_back(4.03644100);
          parameter[i].push_back(-0.70605700);
        } else error("Atom name not known: "+Aname);
      } else if(Rname=="  U") {
        if(Aname=="BB1") {
          parameter[i].push_back(32.88500000);
          parameter[i].push_back(0.08321400);
          parameter[i].push_back(-7.35634900);
          parameter[i].push_back(2.18826800);
          parameter[i].push_back(-3.56047400);
          parameter[i].push_back(2.33064700);
          parameter[i].push_back(-0.39744000);
        } else if(Aname=="BB2") {
          parameter[i].push_back(3.80600000);
          parameter[i].push_back(-0.10773100);
          parameter[i].push_back(9.60099900);
          parameter[i].push_back(-6.26131900);
          parameter[i].push_back(-0.46668300);
          parameter[i].push_back(1.13698100);
          parameter[i].push_back(-0.21351600);
        } else if(Aname=="BB3") {
          parameter[i].push_back(3.59400000);
          parameter[i].push_back(0.04544300);
          parameter[i].push_back(9.59625900);
          parameter[i].push_back(-1.29222200);
          parameter[i].push_back(-7.11143200);
          parameter[i].push_back(4.05687700);
          parameter[i].push_back(-0.63382800);
        } else if(Aname=="SC1") {
          parameter[i].push_back(5.95100000);
          parameter[i].push_back(-0.02924500);
          parameter[i].push_back(2.59668700);
          parameter[i].push_back(-0.56118700);
          parameter[i].push_back(-1.56477100);
          parameter[i].push_back(0.89265100);
          parameter[i].push_back(-0.14130800);
        } else if(Aname=="SC2") {
          parameter[i].push_back(10.90100000);
          parameter[i].push_back(0.02178900);
          parameter[i].push_back(0.18839000);
          parameter[i].push_back(0.72223100);
          parameter[i].push_back(-1.92581600);
          parameter[i].push_back(0.96654300);
          parameter[i].push_back(-0.14501300);
        } else if(Aname=="SC3") {
          parameter[i].push_back(5.24600000);
          parameter[i].push_back(-0.04586500);
          parameter[i].push_back(5.89978100);
          parameter[i].push_back(-1.50664700);
          parameter[i].push_back(-3.17054400);
          parameter[i].push_back(1.93717100);
          parameter[i].push_back(-0.31701000);
        } else if(Aname=="3TE") {
          parameter[i].push_back(4.23000000);
          parameter[i].push_back(0.00067500);
          parameter[i].push_back(0.92102300);
          parameter[i].push_back(0.08100800);
          parameter[i].push_back(-0.39084300);
          parameter[i].push_back(0.12441900);
          parameter[i].push_back(-0.01124900);
        } else if(Aname=="5TE") {
          parameter[i].push_back(4.23000000);
          parameter[i].push_back(0.00059000);
          parameter[i].push_back(0.92154600);
          parameter[i].push_back(0.07968200);
          parameter[i].push_back(-0.38950100);
          parameter[i].push_back(0.12382500);
          parameter[i].push_back(-0.01115100);
        } else if(Aname=="TE3") {
          parameter[i].push_back(7.82400000);
          parameter[i].push_back(-0.02968100);
          parameter[i].push_back(7.93783200);
          parameter[i].push_back(-0.33078100);
          parameter[i].push_back(-10.14120200);
          parameter[i].push_back(6.63334700);
          parameter[i].push_back(-1.22111200);
        } else if(Aname=="TE5") {
          parameter[i].push_back(8.03600000);
          parameter[i].push_back(-0.00909700);
          parameter[i].push_back(4.33193500);
          parameter[i].push_back(0.43416500);
          parameter[i].push_back(-5.80831400);
          parameter[i].push_back(3.52438800);
          parameter[i].push_back(-0.62382400);
        } else error("Atom name not known: "+Aname);
      } else if(Rname==" DA") {
        if(Aname=="BB1") {
          parameter[i].push_back(32.88500000);
          parameter[i].push_back(0.08179900);
          parameter[i].push_back(-7.31735900);
          parameter[i].push_back(2.15614500);
          parameter[i].push_back(-3.52263200);
          parameter[i].push_back(2.30604700);
          parameter[i].push_back(-0.39270100);
        } else if(Aname=="BB2") {
          parameter[i].push_back(3.80600000);
          parameter[i].push_back(-0.10597700);
          parameter[i].push_back(9.52537500);
          parameter[i].push_back(-6.12991000);
          parameter[i].push_back(-0.54092600);
          parameter[i].push_back(1.15429100);
          parameter[i].push_back(-0.21503500);
        } else if(Aname=="BB3") {
          parameter[i].push_back(-1.35600000);
          parameter[i].push_back(0.58928300);
          parameter[i].push_back(6.71894100);
          parameter[i].push_back(4.14050900);
          parameter[i].push_back(-9.65859900);
          parameter[i].push_back(4.43185000);
          parameter[i].push_back(-0.64657300);
        } else if(Aname=="SC1") {
          parameter[i].push_back(6.67100000);
          parameter[i].push_back(-0.00871400);
          parameter[i].push_back(1.63289100);
          parameter[i].push_back(-0.06637700);
          parameter[i].push_back(-1.48632900);
          parameter[i].push_back(0.78551800);
          parameter[i].push_back(-0.12087300);
        } else if(Aname=="SC2") {
          parameter[i].push_back(5.95100000);
          parameter[i].push_back(-0.02634300);
          parameter[i].push_back(2.54864300);
          parameter[i].push_back(-0.49015800);
          parameter[i].push_back(-1.55386900);
          parameter[i].push_back(0.86630200);
          parameter[i].push_back(-0.13546200);
        } else if(Aname=="SC3") {
          parameter[i].push_back(11.39400000);
          parameter[i].push_back(0.00859500);
          parameter[i].push_back(-0.25471400);
          parameter[i].push_back(0.48718800);
          parameter[i].push_back(-1.74520000);
          parameter[i].push_back(0.99246200);
          parameter[i].push_back(-0.16351900);
        } else if(Aname=="SC4") {
          parameter[i].push_back(6.45900000);
          parameter[i].push_back(0.01991800);
          parameter[i].push_back(4.17962300);
          parameter[i].push_back(0.97469100);
          parameter[i].push_back(-5.02950400);
          parameter[i].push_back(2.55371800);
          parameter[i].push_back(-0.39113400);
        } else if(Aname=="3TE") {
          parameter[i].push_back(4.23000000);
          parameter[i].push_back(0.00062600);
          parameter[i].push_back(0.92142000);
          parameter[i].push_back(0.08016400);
          parameter[i].push_back(-0.39000300);
          parameter[i].push_back(0.12402500);
          parameter[i].push_back(-0.01117900);
        } else if(Aname=="5TE") {
          parameter[i].push_back(4.23000000);
          parameter[i].push_back(0.00055500);
          parameter[i].push_back(0.92183900);
          parameter[i].push_back(0.07907600);
          parameter[i].push_back(-0.38895100);
          parameter[i].push_back(0.12359600);
          parameter[i].push_back(-0.01111600);
        } else if(Aname=="TE3") {
          parameter[i].push_back(2.87400000);
          parameter[i].push_back(0.00112900);
          parameter[i].push_back(12.51167200);
          parameter[i].push_back(-7.67548000);
          parameter[i].push_back(-2.02234000);
          parameter[i].push_back(2.50837100);
          parameter[i].push_back(-0.49458500);
        } else if(Aname=="TE5") {
          parameter[i].push_back(8.03600000);
          parameter[i].push_back(0.00473100);
          parameter[i].push_back(4.65554400);
          parameter[i].push_back(0.66424100);
          parameter[i].push_back(-6.62131300);
          parameter[i].push_back(3.96107400);
          parameter[i].push_back(-0.69075800);
        } else error("Atom name not known: "+Aname);
      } else if(Rname==" DC") {
        if(Aname=="BB1") {
          parameter[i].push_back(32.88500000);
          parameter[i].push_back(0.08189900);
          parameter[i].push_back(-7.32493500);
          parameter[i].push_back(2.15976900);
          parameter[i].push_back(-3.52612100);
          parameter[i].push_back(2.31058600);
          parameter[i].push_back(-0.39402700);
        } else if(Aname=="BB2") {
          parameter[i].push_back(3.80600000);
          parameter[i].push_back(-0.10559800);
          parameter[i].push_back(9.52527700);
          parameter[i].push_back(-6.12131700);
          parameter[i].push_back(-0.54899400);
          parameter[i].push_back(1.15592900);
          parameter[i].push_back(-0.21494500);
        } else if(Aname=="BB3") {
          parameter[i].push_back(-1.35600000);
          parameter[i].push_back(0.55525700);
          parameter[i].push_back(6.80305500);
          parameter[i].push_back(4.05924700);
          parameter[i].push_back(-9.61034700);
          parameter[i].push_back(4.41253800);
          parameter[i].push_back(-0.64315100);
        } else if(Aname=="SC1") {
          parameter[i].push_back(5.95100000);
          parameter[i].push_back(-0.02899900);
          parameter[i].push_back(2.59587800);
          parameter[i].push_back(-0.55388300);
          parameter[i].push_back(-1.56395100);
          parameter[i].push_back(0.88967400);
          parameter[i].push_back(-0.14062500);
        } else if(Aname=="SC2") {
          parameter[i].push_back(11.62100000);
          parameter[i].push_back(0.01358100);
          parameter[i].push_back(-0.24913000);
          parameter[i].push_back(0.48787200);
          parameter[i].push_back(-1.52867300);
          parameter[i].push_back(0.83694900);
          parameter[i].push_back(-0.13395300);
        } else if(Aname=="SC3") {
          parameter[i].push_back(5.01900000);
          parameter[i].push_back(-0.03298400);
          parameter[i].push_back(5.54242800);
          parameter[i].push_back(-0.96081500);
          parameter[i].push_back(-3.71051600);
          parameter[i].push_back(2.16500200);
          parameter[i].push_back(-0.35023400);
        } else if(Aname=="3TE") {
          parameter[i].push_back(4.23000000);
          parameter[i].push_back(0.00055700);
          parameter[i].push_back(0.92181400);
          parameter[i].push_back(0.07924000);
          parameter[i].push_back(-0.38916400);
          parameter[i].push_back(0.12369900);
          parameter[i].push_back(-0.01113300);
        } else if(Aname=="5TE") {
          parameter[i].push_back(4.23000000);
          parameter[i].push_back(0.00066500);
          parameter[i].push_back(0.92103900);
          parameter[i].push_back(0.08064600);
          parameter[i].push_back(-0.39034900);
          parameter[i].push_back(0.12417600);
          parameter[i].push_back(-0.01120600);
        } else if(Aname=="TE3") {
          parameter[i].push_back(2.87400000);
          parameter[i].push_back(-0.05235500);
          parameter[i].push_back(13.09201200);
          parameter[i].push_back(-9.48128200);
          parameter[i].push_back(-0.14958600);
          parameter[i].push_back(1.75537200);
          parameter[i].push_back(-0.39347500);
        } else if(Aname=="TE5") {
          parameter[i].push_back(8.03600000);
          parameter[i].push_back(-0.00513600);
          parameter[i].push_back(4.67705700);
          parameter[i].push_back(0.48333300);
          parameter[i].push_back(-6.34511000);
          parameter[i].push_back(3.83388500);
          parameter[i].push_back(-0.67367800);
        } else error("Atom name not known: "+Aname);
      } else if(Rname==" DG") {
        if(Aname=="BB1") {
          parameter[i].push_back(32.88500000);
          parameter[i].push_back(0.08182900);
          parameter[i].push_back(-7.32133900);
          parameter[i].push_back(2.15767900);
          parameter[i].push_back(-3.52369700);
          parameter[i].push_back(2.30839600);
          parameter[i].push_back(-0.39348300);
        } else if(Aname=="BB2") {
          parameter[i].push_back(3.80600000);
          parameter[i].push_back(-0.10618100);
          parameter[i].push_back(9.54169000);
          parameter[i].push_back(-6.15177600);
          parameter[i].push_back(-0.53462400);
          parameter[i].push_back(1.15581300);
          parameter[i].push_back(-0.21567000);
        } else if(Aname=="BB3") {
          parameter[i].push_back(-1.35600000);
          parameter[i].push_back(0.57489100);
          parameter[i].push_back(6.75164700);
          parameter[i].push_back(4.11300900);
          parameter[i].push_back(-9.63394600);
          parameter[i].push_back(4.41675400);
          parameter[i].push_back(-0.64339900);
        } else if(Aname=="SC1") {
          parameter[i].push_back(6.67100000);
          parameter[i].push_back(-0.00886600);
          parameter[i].push_back(1.63333000);
          parameter[i].push_back(-0.06892100);
          parameter[i].push_back(-1.48683500);
          parameter[i].push_back(0.78670800);
          parameter[i].push_back(-0.12113900);
        } else if(Aname=="SC2") {
          parameter[i].push_back(11.39400000);
          parameter[i].push_back(0.00907900);
          parameter[i].push_back(-0.22475500);
          parameter[i].push_back(0.49535100);
          parameter[i].push_back(-1.75324900);
          parameter[i].push_back(0.98767400);
          parameter[i].push_back(-0.16150800);
        } else if(Aname=="SC3") {
          parameter[i].push_back(10.90100000);
          parameter[i].push_back(0.02207600);
          parameter[i].push_back(0.17932200);
          parameter[i].push_back(0.73253200);
          parameter[i].push_back(-1.95554900);
          parameter[i].push_back(0.98339900);
          parameter[i].push_back(-0.14763600);
        } else if(Aname=="SC4") {
          parameter[i].push_back(6.45900000);
          parameter[i].push_back(0.02018400);
          parameter[i].push_back(4.17705400);
          parameter[i].push_back(0.98531700);
          parameter[i].push_back(-5.04354900);
          parameter[i].push_back(2.56123700);
          parameter[i].push_back(-0.39249300);
        } else if(Aname=="3TE") {
          parameter[i].push_back(4.23000000);
          parameter[i].push_back(0.00061700);
          parameter[i].push_back(0.92140100);
          parameter[i].push_back(0.08016400);
          parameter[i].push_back(-0.39003500);
          parameter[i].push_back(0.12406900);
          parameter[i].push_back(-0.01119200);
        } else if(Aname=="5TE") {
          parameter[i].push_back(4.23000000);
          parameter[i].push_back(0.00064900);
          parameter[i].push_back(0.92110500);
          parameter[i].push_back(0.08031500);
          parameter[i].push_back(-0.38997000);
          parameter[i].push_back(0.12401200);
          parameter[i].push_back(-0.01118100);
        } else if(Aname=="TE3") {
          parameter[i].push_back(2.87400000);
          parameter[i].push_back(0.00182000);
          parameter[i].push_back(12.41507000);
          parameter[i].push_back(-7.47384800);
          parameter[i].push_back(-2.11864700);
          parameter[i].push_back(2.50112600);
          parameter[i].push_back(-0.48652200);
        } else if(Aname=="TE5") {
          parameter[i].push_back(8.03600000);
          parameter[i].push_back(0.00676400);
          parameter[i].push_back(4.65989200);
          parameter[i].push_back(0.78482500);
          parameter[i].push_back(-6.86460600);
          parameter[i].push_back(4.11675400);
          parameter[i].push_back(-0.72249100);
        } else error("Atom name not known: "+Aname);
      } else if(Rname==" DT") {
        if(Aname=="BB1") {
          parameter[i].push_back(32.88500000);
          parameter[i].push_back(0.08220100);
          parameter[i].push_back(-7.33006800);
          parameter[i].push_back(2.16636500);
          parameter[i].push_back(-3.53465700);
          parameter[i].push_back(2.31447600);
          parameter[i].push_back(-0.39445400);
        } else if(Aname=="BB2") {
          parameter[i].push_back(3.80600000);
          parameter[i].push_back(-0.10723000);
          parameter[i].push_back(9.56675000);
          parameter[i].push_back(-6.20236100);
          parameter[i].push_back(-0.49550400);
          parameter[i].push_back(1.14300600);
          parameter[i].push_back(-0.21420000);
        } else if(Aname=="BB3") {
          parameter[i].push_back(-1.35600000);
          parameter[i].push_back(0.56737900);
          parameter[i].push_back(6.76595400);
          parameter[i].push_back(4.08976100);
          parameter[i].push_back(-9.61512500);
          parameter[i].push_back(4.40975100);
          parameter[i].push_back(-0.64239800);
        } else if(Aname=="SC1") {
          parameter[i].push_back(5.95100000);
          parameter[i].push_back(-0.02926500);
          parameter[i].push_back(2.59630300);
          parameter[i].push_back(-0.56152200);
          parameter[i].push_back(-1.56532600);
          parameter[i].push_back(0.89322800);
          parameter[i].push_back(-0.14142900);
        } else if(Aname=="SC2") {
          parameter[i].push_back(10.90100000);
          parameter[i].push_back(0.02183400);
          parameter[i].push_back(0.19463000);
          parameter[i].push_back(0.72393000);
          parameter[i].push_back(-1.93199500);
          parameter[i].push_back(0.96856300);
          parameter[i].push_back(-0.14512600);
        } else if(Aname=="SC3") {
          parameter[i].push_back(4.31400000);
          parameter[i].push_back(-0.07745600);
          parameter[i].push_back(12.49820300);
          parameter[i].push_back(-7.64994200);
          parameter[i].push_back(-3.00359600);
          parameter[i].push_back(3.26263300);
          parameter[i].push_back(-0.64498600);
        } else if(Aname=="3TE") {
          parameter[i].push_back(4.23000000);
          parameter[i].push_back(0.00062000);
          parameter[i].push_back(0.92141100);
          parameter[i].push_back(0.08030900);
          parameter[i].push_back(-0.39021500);
          parameter[i].push_back(0.12414000);
          parameter[i].push_back(-0.01120100);
        } else if(Aname=="5TE") {
          parameter[i].push_back(4.23000000);
          parameter[i].push_back(0.00063700);
          parameter[i].push_back(0.92130800);
          parameter[i].push_back(0.08026900);
          parameter[i].push_back(-0.39007500);
          parameter[i].push_back(0.12406600);
          parameter[i].push_back(-0.01118800);
        } else if(Aname=="TE3") {
          parameter[i].push_back(2.87400000);
          parameter[i].push_back(-0.00251200);
          parameter[i].push_back(12.43576400);
          parameter[i].push_back(-7.55343800);
          parameter[i].push_back(-2.07363500);
          parameter[i].push_back(2.51279300);
          parameter[i].push_back(-0.49437100);
        } else if(Aname=="TE5") {
          parameter[i].push_back(8.03600000);
          parameter[i].push_back(0.00119900);
          parameter[i].push_back(4.91762300);
          parameter[i].push_back(0.65637000);
          parameter[i].push_back(-7.23392500);
          parameter[i].push_back(4.44636600);
          parameter[i].push_back(-0.79467800);
        } else error("Atom name not known: "+Aname);
      } else error("Residue not known: "+Rname);
    }
  } else {
    error("MOLINFO DATA not found\n");
  }
}

double SAXS::calculateASF(const vector<AtomNumber> &atoms, vector<vector<long double> > &FF_tmp, const double rho)
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

  double Iq0=0.;
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
